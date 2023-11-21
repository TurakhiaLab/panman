#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <unordered_set>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <ctime>
#include <iomanip>
#include <mutex>
#include "kseq.hpp"
#include "PangenomeMATV2.hpp"

KSEQ_INIT(int, read)


std::string PangenomeMAT2::getDate(){
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);
    std::string date;
    date += std::to_string(now->tm_year + 1900)
         + std::to_string(now->tm_mon + 1)
         +  std::to_string(now->tm_mday);
    return date;
}

PangenomeMAT2::Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

PangenomeMAT2::Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

PangenomeMAT2::Block::Block(MATNew::block b){
    primaryBlockId = (b.blockid() >> 32);
    if(b.blockgapexist()){
        secondaryBlockId = (b.blockid() & 0xFFFFFFFF);
    } else {
        secondaryBlockId = -1;
    }
    
    chromosomeName = b.chromosomename();
    for(int i = 0; i < b.consensusseq_size(); i++){
        consensusSeq.push_back(b.consensusseq(i));
    }
}

void PangenomeMAT2::stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}
void stringSplitStr(std::string s, std::string delimiter, std::vector<std::string>& words) {
    int32_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        words.push_back (token);
    }

    words.push_back (s.substr (pos_start));
}

PangenomeMAT2::Node* PangenomeMAT2::Tree::createTreeFromNewickString(std::string newickString) {

    newickString = PangenomeMAT2::stripString(newickString);

    PangenomeMAT2::Node* newTreeRoot;

    std::vector<std::string> leaves;
    std::vector<size_t> numOpen;
    std::vector<size_t> numClose;
    std::vector<std::queue<float>> branchLen (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    stringSplit(newickString, ',', s1);

    numOpen.reserve(s1.size());
    numClose.reserve(s1.size());

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        size_t leafDepth = 0;

        bool stop = false;
        bool branchStart = false;
        std::string leaf = "";
        std::string branch = "";

        for (auto c: s) {
            if (c == ':') {
                stop = true;
                branch = "";
                branchStart = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branchLen.size() <= level) {
                    branchLen.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                branchLen[level].push(len);
                level--;
                branchStart = false;
            } else if (!stop) {
                leaf += c;
                branchStart = false;
                leafDepth = level;

            } else if (branchStart) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        leaves.push_back(std::move(leaf));
        numOpen.push_back(no);
        numClose.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        branchLen[level].push(len);

        // Adjusting max and mean depths
        m_maxDepth = std::max(m_maxDepth, leafDepth);
        m_meanDepth += leafDepth;

    }

    m_meanDepth /= leaves.size();

    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    m_numLeaves = leaves.size();

    std::stack<Node*> parentStack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = numOpen[i];
        auto nc = numClose[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = newInternalNodeId();
            Node* newNode = nullptr;
            if (parentStack.size() == 0) {
                newNode = new Node(nid, branchLen[level].front());
                newTreeRoot = newNode;
            } else {
                newNode = new Node(nid, parentStack.top(), branchLen[level].front());
            }
            branchLen[level].pop();
            level++;

            allNodes[nid] = newNode;
            parentStack.push(newNode);
        }
        Node* leafNode = new Node(leaf, parentStack.top(), branchLen[level].front());
        allLeaves.push_back(leafNode);

        allNodes[leaf] = leafNode;

        branchLen[level].pop();
        for (size_t j=0; j<nc; j++) {
            parentStack.pop();
            level--;
        }
    }

    if (newTreeRoot == nullptr) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    return newTreeRoot;
}

void PangenomeMAT2::Tree::assignMutationsToNodes(Node* root, size_t& currentIndex, std::vector< MATNew::node >& nodes){
    std::vector< PangenomeMAT2::NucMut > storedNucMutation;

    for(int i = 0; i < nodes[currentIndex].nucmutation_size(); i++){
        storedNucMutation.push_back( PangenomeMAT2::NucMut(nodes[currentIndex].nucmutation(i)) );
    }

    std::vector< PangenomeMAT2::BlockMut > storedBlockMutation;
    for(int i = 0; i < nodes[currentIndex].blockmutation_size(); i++){
        PangenomeMAT2::BlockMut tempBlockMut;
        tempBlockMut.loadFromProtobuf(nodes[currentIndex].blockmutation(i));
        storedBlockMutation.push_back(tempBlockMut);
    }

    for(int i = 0; i < nodes[currentIndex].annotations_size(); i++){
        root->annotations.push_back(nodes[currentIndex].annotations(i));
        annotationsToNodes[nodes[currentIndex].annotations(i)].push_back(root->identifier);
    }

    root->nucMutation = storedNucMutation;
    root->blockMutation = storedBlockMutation;

    for(auto child: root->children){
        currentIndex++;
        assignMutationsToNodes(child, currentIndex, nodes);
    }

}

void PangenomeMAT2::Tree::invertTree(PangenomeMAT2::Node* root){
    for(auto child: root->children){
        invertTree(child);
    }
    std::reverse(root->children.begin(), root->children.end());
}

PangenomeMAT2::Tree::Tree(std::ifstream& fin){

    MATNew::tree mainTree;

    if(!mainTree.ParseFromIstream(&fin)){
        throw std::invalid_argument("Could not read tree from input file.");
    }

    // Create tree
    root = createTreeFromNewickString(mainTree.newick());
    invertTree(root);

    std::vector< MATNew::node > storedNodes;
    for(int i = 0; i < mainTree.nodes_size(); i++){
        storedNodes.push_back(mainTree.nodes(i));
    }

    size_t initialIndex = 0;

    assignMutationsToNodes(root, initialIndex, storedNodes);

    // Block sequence
    for(int i = 0; i < mainTree.blocks_size(); i++){
        blocks.emplace_back(mainTree.blocks(i));
    }

    // Gap List
    for(int i = 0; i < mainTree.gaps_size(); i++){
        PangenomeMAT2::GapList tempGaps;
        tempGaps.primaryBlockId = (mainTree.gaps(i).blockid() >> 32);
        tempGaps.secondaryBlockId = (mainTree.gaps(i).blockgapexist() ? (mainTree.gaps(i).blockid() & 0xFFFF): -1);
        for(int j = 0; j < mainTree.gaps(i).nucposition_size(); j++){
            tempGaps.nucPosition.push_back(mainTree.gaps(i).nucposition(j));
            tempGaps.nucGapLength.push_back(mainTree.gaps(i).nucgaplength(j));
        }
        gaps.push_back(tempGaps);
    }

    // Block gap list
    for(int i = 0; i < mainTree.blockgaps().blockposition_size(); i++){
        blockGaps.blockPosition.push_back(mainTree.blockgaps().blockposition(i));
        blockGaps.blockGapLength.push_back(mainTree.blockgaps().blockgaplength(i));
    }

    // Setting up global coordinates for fast retrieval
    setupGlobalCoordinates();

}

int getTotalParsimonyParallelHelper(PangenomeMAT2::Node* root, PangenomeMAT2::NucMutationType nucMutType, PangenomeMAT2::BlockMutationType blockMutType){
    int totalMutations = 0;

    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->nucMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++){
            
            if(((root->nucMutation[i].mutInfo) & 0x7) == nucMutType){
                if(nucMutType == PangenomeMAT2::NucMutationType::NS){
                    init += ((root->nucMutation[i].mutInfo) >> 4); // Length of contiguous mutation in case of substitution
                } else {
                    init++;
                }
            }
        }
        return init;
    }, [&](int x, int y){
        return x + y;
    });

    if(blockMutType != PangenomeMAT2::BlockMutationType::NONE){
        totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
            for(int i = r.begin(); i != r.end(); i++){
                if(root->blockMutation[i].blockMutInfo == blockMutType){
                    init++;
                }
            }
            return init;
        }, [&](int x, int y){
            return x + y;
        });
    }


    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->children.size()), 0, [&](tbb::blocked_range<int>& r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++){
            init += getTotalParsimonyParallelHelper(root->children[i], nucMutType, blockMutType);
        }
        return init;
    },
    [](int x, int y) -> int {
        return x+y;
    });

    return totalMutations;
}

int PangenomeMAT2::Tree::getTotalParsimonyParallel(NucMutationType nucMutType, BlockMutationType blockMutType){

    return getTotalParsimonyParallelHelper(root, nucMutType, blockMutType);

}

void PangenomeMAT2::Tree::printSummary(){

    std::cout << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    std::cout << "Total Samples in Tree: " << m_numLeaves << std::endl;
    std::cout << "Total Substitutions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::NI, PangenomeMAT2::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::ND, PangenomeMAT2::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP Substitutions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::NSNPS) << std::endl;
    std::cout << "Total SNP Insertions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::NSNPI) << std::endl;
    std::cout << "Total SNP Deletions: " << getTotalParsimonyParallel(PangenomeMAT2::NucMutationType::NSNPD) << std::endl;
    std::cout << "Max Tree Depth: " << m_maxDepth << std::endl;
    std::cout << "Mean Tree Depth: " << m_meanDepth << std::endl;

}

void PangenomeMAT2::Tree::printBfs(Node* node){
    if(node == nullptr){
        node = root;
    }

    // Traversal test
    std::queue<Node *> bfsQueue;
    size_t prevLev = 0;
    
    bfsQueue.push(node);

    while(!bfsQueue.empty()){
        Node* current = bfsQueue.front();
        bfsQueue.pop();

        if(current->level != prevLev){
            std::cout << '\n';
            prevLev = current->level;
        }
        std::cout << '(' << current->identifier << "," << current->branchLength << ") ";

        for(auto child: current->children){
            bfsQueue.push(child);
        }
    }
    std::cout << '\n';
}

void PangenomeMAT2::printSequenceLines(const std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
    const std::vector< std::pair< bool, std::vector< bool > > >& blockExists, size_t lineSize, bool aligned, std::ofstream& fout){

    std::string line;

    for(size_t i = 0; i < blockExists.size(); i++){
        for(size_t j = 0; j < blockExists[i].second.size(); j++){
            if(blockExists[i].second[j]){

                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        if(sequence[i].second[j][k].second[w] != '-'){
                            line += sequence[i].second[j][k].second[w];
                        } else if(aligned){
                            line += '-';
                        }
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                    if(sequence[i].second[j][k].first != 'x'){
                        if(sequence[i].second[j][k].first != '-'){
                            line += sequence[i].second[j][k].first;
                        } else if(aligned){
                            line += '-';
                        }
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                }
            } else if(aligned) {
                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        line += '-';
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                    if(sequence[i].second[j][k].first != 'x'){
                        line += '-';
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                }
            }
        }

        if(blockExists[i].first){
            for(size_t j = 0; j < sequence[i].first.size(); j++){
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++){
                    if(sequence[i].first[j].second[k] != '-'){
                        line += sequence[i].first[j].second[k];
                    } else if(aligned){
                        line += '-';
                    }
                    if(line.length() == lineSize){
                        fout << line << '\n';
                        line = "";
                    }
                }
                if(sequence[i].first[j].first != 'x'){
                    if(sequence[i].first[j].first != '-'){
                        line += sequence[i].first[j].first;
                    } else if(aligned){
                        line += '-';
                    }
                    if(line.length() == lineSize){
                        fout << line << '\n';
                        line = "";
                    }
                }
            }
        }

    }


    if(line.length()){
        fout << line << '\n';
        line = "";
    }
}


std::string getSequence(const std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
    const std::vector< std::pair< bool, std::vector< bool > > >& blockExists, bool aligned){

    std::string fastaSequence = "";

    for(size_t i = 0; i < blockExists.size(); i++){
        for(size_t j = 0; j < blockExists[i].second.size(); j++){
            if(blockExists[i].second[j]){

                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        if(sequence[i].second[j][k].second[w] != '-'){
                            fastaSequence += sequence[i].second[j][k].second[w];
                        } else if(aligned){
                            fastaSequence += '-';
                        }
                    }
                    if(sequence[i].second[j][k].first != 'x'){
                        if(sequence[i].second[j][k].first != '-'){
                            fastaSequence += sequence[i].second[j][k].first;
                        } else if(aligned){
                            fastaSequence += '-';
                        }
                    }
                }
            } else if(aligned) {
                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        fastaSequence += '-';
                    }
                    if(sequence[i].second[j][k].first != 'x'){
                        fastaSequence += '-';
                    }
                }
            }
        }

        if(blockExists[i].first){
            for(size_t j = 0; j < sequence[i].first.size(); j++){
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++){
                    if(sequence[i].first[j].second[k] != '-'){
                        fastaSequence += sequence[i].first[j].second[k];
                    } else if(aligned){
                        fastaSequence += '-';
                    }
                }
                if(sequence[i].first[j].first != 'x'){
                    if(sequence[i].first[j].first != '-'){
                        fastaSequence += sequence[i].first[j].first;
                    } else if(aligned){
                        fastaSequence += '-';
                    }
                }
            }
        }

    }
    

    return fastaSequence;

}

char PangenomeMAT2::getNucleotideFromCode(int code){
    switch(code){
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 5:
            return 'R';
        case 10:
            return 'Y';
        case 6:
            return 'S';
        case 9:
            return 'W';
        case 12:
            return 'K';
        case 3:
            return 'M';
        case 14:
            return 'B';
        case 13:
            return 'D';
        case 11:
            return 'H';
        case 7:
            return 'V';
        default:
            return 'N';
    }
}

int reverseNucs(int nucs){
    int res = 0;
    for(int i = 0; i < 6; i++){
        res  = (res ^ (((nucs >> (4*i)) & 0xF) << 4*(5-i)));
    }
    return res;
}

void printFASTAHelper(PangenomeMAT2::Node* root,\
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
    std::vector< std::pair< bool, std::vector< bool > > >& blockExists,\
    std::ofstream& fout, bool aligned = false){

    // Apply mutations
    // Block mutations - ignored for now since the block IDs don't seem right in the files

    std::vector< std::tuple< int32_t, int32_t, bool, bool > > blockMutationInfo;

    for(auto mutation: root->blockMutation){
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = (mutation.blockMutInfo);

        if(type == 1){
            bool oldVal;
            if(secondaryBlockId != -1){
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = true;
            } else {
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = true;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, true) );
        } else {
            bool oldVal;
            if(secondaryBlockId != -1){
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = false;
            } else {
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = false;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, false) );
        }

    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++){
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3){
            // Either S, I or D

            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(type == PangenomeMAT2::NucMutationType::NS){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NI){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::ND){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }
                    }
                }
            }
        } 
        else {
            if(type == PangenomeMAT2::NucMutationType::NSNPS){
                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPI){
                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPD){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    }
                }
            }
        }
    }
    
        // Print sequence
    fout << '>' << root->identifier << std::endl;

    PangenomeMAT2::printSequenceLines(sequence, blockExists, 70, aligned, fout);

    if(root->children.size() > 0){
        // DFS on children
        for(PangenomeMAT2::Node* child: root->children){
            printFASTAHelper(child, sequence, blockExists, fout, aligned);
        }
    }

    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
        }
    }

    // Undo nuc mutations
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }
}

void PangenomeMAT2::Tree::printFASTA(std::ofstream& fout, bool aligned){
    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(blocks.size() + 1, {false, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    }

    int32_t maxBlockId = 0;

    for(size_t i = 0; i < blocks.size(); i++){
        
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                
                if(secondaryBlockId != -1){
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }

            if(endFlag){
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1){
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }

    printFASTAHelper(root, sequence, blockExists, fout, aligned);

}

// Merge parent node and child node into parent node
void PangenomeMAT2::Tree::mergeNodes(PangenomeMAT2::Node* par, PangenomeMAT2::Node* chi){
    
    par->identifier = chi->identifier;
    par->branchLength += chi->branchLength;
    par->children = chi->children;

    // For block mutations, we cancel out irrelevant mutations
    std::map< std::pair<int, int>, PangenomeMAT2::BlockMutationType > bidMutations;

    for(auto mutation: par->blockMutation){
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;

        bool type = (mutation.blockMutInfo);
        if(type == PangenomeMAT2::BlockMutationType::BI){
            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = PangenomeMAT2::BlockMutationType::BI;
        } else {
            if(bidMutations.find(std::make_pair(primaryBlockId, secondaryBlockId)) != bidMutations.end()){
                if(bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] == PangenomeMAT2::BlockMutationType::BI){
                    // If it was insertion earlier, cancel out
                    bidMutations.erase(std::make_pair(primaryBlockId, secondaryBlockId));
                }
                // Otherwise, it remains deletion
            } else {
                bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = PangenomeMAT2::BlockMutationType::BD;
            }
        }
    }

    for(auto mutation: chi->blockMutation){
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;

        int type = (mutation.blockMutInfo);
        if(type == PangenomeMAT2::BlockMutationType::BI){
            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = PangenomeMAT2::BlockMutationType::BI;
        } else {
            if(bidMutations.find(std::make_pair(primaryBlockId, secondaryBlockId)) != bidMutations.end()){
                if(bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] == PangenomeMAT2::BlockMutationType::BI){
                    // If it was insertion earlier, cancel out
                    bidMutations.erase(std::make_pair(primaryBlockId, secondaryBlockId));
                }
                // Otherwise, it remains deletion
            } else {
                bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = PangenomeMAT2::BlockMutationType::BD;
            }
        }
    }

    std::vector< PangenomeMAT2::BlockMut > newBlockMutation;
    for(auto mutation: bidMutations){
        if(mutation.second == PangenomeMAT2::BlockMutationType::BI){
            PangenomeMAT2::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = PangenomeMAT2::BlockMutationType::BI;
            newBlockMutation.push_back( tempBlockMut );
        } else {
            PangenomeMAT2::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = PangenomeMAT2::BlockMutationType::BD;
            newBlockMutation.push_back( tempBlockMut );
        }
    }

    par->blockMutation = newBlockMutation;

    for(auto mutation: chi->nucMutation){
        par->nucMutation.push_back(mutation);
    }

    delete chi;
}

// Replace old type, char pair with new type char pair
std::pair< int, int > PangenomeMAT2::replaceMutation(std::pair<int,int> oldMutation, std::pair<int, int> newMutation){
    std::pair<int, int> ans = newMutation;
    if(oldMutation.first == newMutation.first){
        ans = newMutation;
    } else if(oldMutation.first == PangenomeMAT2::NucMutationType::NSNPS){
        // Insertion after substitution (doesn't make sense but just in case)
        if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPI){
            ans.first = PangenomeMAT2::NucMutationType::NSNPS;
        } else if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPD){
            ans = newMutation;
        }
    } else if(oldMutation.first == PangenomeMAT2::NucMutationType::NSNPI){
        if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPS){
            ans.first = PangenomeMAT2::NucMutationType::NSNPI;
        } else if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPD){
            // Cancel out the two mutations if deletion after insertion
            ans = std::make_pair(404, 404);
        }
    } else if(oldMutation.first == PangenomeMAT2::NucMutationType::NSNPD){
        if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPI){
            ans.first = PangenomeMAT2::NucMutationType::NSNPS;
        } else if(newMutation.first == PangenomeMAT2::NucMutationType::NSNPS){
            // Substitution after deletion. Doesn't make sense but still
            ans.first = PangenomeMAT2::NucMutationType::NSNPI;
        }
    }
    return ans;
}

bool PangenomeMAT2::Tree::debugSimilarity(const std::vector< PangenomeMAT2::NucMut > array1, const std::vector< PangenomeMAT2::NucMut > array2){
    std::map< std::tuple< int, int, int, int >, std::pair< int, int > > mutationRecords1, mutationRecords2;

    for(auto mutation: array1){
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = ((mutation.mutInfo) >> 4);

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT2::NucMutationType::NS:
                newType = PangenomeMAT2::NucMutationType::NSNPS;
                break;
            case PangenomeMAT2::NucMutationType::ND:
                newType = PangenomeMAT2::NucMutationType::NSNPD;
                break;
            case PangenomeMAT2::NucMutationType::NI:
                newType = PangenomeMAT2::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar = (((mutation.nucs) >> (4*(5 - i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords1.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords1.end()){
                    mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords1.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords1.end()){
                    mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    for(auto mutation: array2){
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = ((mutation.mutInfo) >> 4);

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT2::NucMutationType::NS:
                newType = PangenomeMAT2::NucMutationType::NSNPS;
                break;
            case PangenomeMAT2::NucMutationType::ND:
                newType = PangenomeMAT2::NucMutationType::NSNPD;
                break;
            case PangenomeMAT2::NucMutationType::NI:
                newType = PangenomeMAT2::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar = (((mutation.nucs) >> (4*(5 - i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords2.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords2.end()){
                    mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords2.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords2.end()){
                    mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    std::vector< std::tuple< int, int, int, int, int, int > > mutationArray1, mutationArray2;
    for(auto u: mutationRecords1){
        mutationArray1.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }
    for(auto u: mutationRecords2){
        mutationArray2.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }

    if(mutationArray1.size() != mutationArray2.size()){
        std::cout << "sizes don't match " << mutationArray1.size() << " " << mutationArray2.size() << std::endl;
        return false;
    }

    for(size_t i = 0; i < mutationArray1.size(); i++){
        if(mutationArray1[i] != mutationArray2[i]){
            std::cout << i << "th index doesn't match" << std::endl;
            return false;
        }
    }

    return true;
}

std::vector< PangenomeMAT2::NucMut > consolidateNucMutations(const std::vector< PangenomeMAT2::NucMut >& nucMutation){
    // primaryBid, secondaryBid, pos, gap_pos -> type, nuc
    std::map< std::tuple< int32_t, int32_t, int32_t, int32_t >, std::pair< int, int > > mutationRecords;
    for(auto mutation: nucMutation){
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = (((mutation.mutInfo) >> 4));

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT2::NucMutationType::NS:
                newType = PangenomeMAT2::NucMutationType::NSNPS;
                break;
            case PangenomeMAT2::NucMutationType::ND:
                newType = PangenomeMAT2::NucMutationType::NSNPD;
                break;
            case PangenomeMAT2::NucMutationType::NI:
                newType = PangenomeMAT2::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar = (((mutation.nucs) >> (4*(5-i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords.end()){
                    mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = PangenomeMAT2::replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords.end()){
                    mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = PangenomeMAT2::replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
    std::vector< std::tuple< int, int, int, int, int, int > > mutationArray;
    for(auto u: mutationRecords){
        mutationArray.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }
    
    // mutation array is already sorted since mutationRecord was sorted
    std::vector< PangenomeMAT2::NucMut > consolidatedMutationArray;

    for(size_t i = 0; i < mutationArray.size(); i++){
        size_t j = i + 1;
        for(; j < std::min(i + 6, mutationArray.size()); j++){
            if(std::get<3>(mutationArray[i]) != -1){
                // gapPos exists
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && std::get<1>(mutationArray[i]) == std::get<1>(mutationArray[j]) && std::get<2>(mutationArray[i]) == std::get<2>(mutationArray[j])
                    && std::get<4>(mutationArray[i]) == std::get<4>(mutationArray[j]) && (size_t)(std::get<3>(mutationArray[j]) - std::get<3>(mutationArray[i])) == j - i)){
                    break;
                }
            } else {
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && std::get<1>(mutationArray[i]) == std::get<1>(mutationArray[j]) && (size_t)(std::get<2>(mutationArray[j]) - std::get<2>(mutationArray[i])) == j - i
                    && std::get<4>(mutationArray[i]) == std::get<4>(mutationArray[j]) && std::get<3>(mutationArray[j]) == std::get<3>(mutationArray[i]))){
                    break;
                }
            }
        }

        if(j - i <= 1){
            consolidatedMutationArray.push_back(PangenomeMAT2::NucMut(mutationArray[i]));
            continue;
        }
        // combine mutations from i to j
        auto newMutation = PangenomeMAT2::NucMut(mutationArray, i, j);

        consolidatedMutationArray.push_back(newMutation);

        i = j - 1;
    }

    return consolidatedMutationArray;

}

void PangenomeMAT2::Tree::dfsExpansion(PangenomeMAT2::Node* node, std::vector< PangenomeMAT2::Node* >& vec){
    vec.push_back(node);
    for(auto child: node->children){
        dfsExpansion(child, vec);
    }
}

std::string PangenomeMAT2::Tree::getNewickString(Node* node){
    invertTree(node);

    std::vector< PangenomeMAT2::Node* > traversal;
    dfsExpansion(node, traversal);

    std::string newick;

    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->branchLength;
        
        if(curr_level < level){
            if (!prev_open) {
                newick += ',';
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                newick += '(';
                prev_open = true;
            }
            if (n->children.size() == 0) {

                newick += n->identifier;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += branch_length;
                }
                prev_open = false;
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        } else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                newick += ')';

                newick += node_stack.top();

                if (branch_length_stack.top() >= 0) {
                    newick += ':';
                    newick += branch_length_stack.top();
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->children.size() == 0) {
                
                newick += ',';
                newick += n->identifier;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += branch_length;
                }
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        } else {
            prev_open = false;
            if (n->children.size() == 0) {
                
                newick += ',';
                newick += n->identifier;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += branch_length;
                }
            } else {
                node_stack.push(n->identifier);
                branch_length_stack.push(branch_length);
            }
        }
        curr_level = level;
    }
    size_t remaining = node_stack.size();
    for (size_t i = 0; i < remaining; i++) {
        newick += ')';
        newick += node_stack.top();
        
        if (branch_length_stack.top() >= 0) {
            newick += ':';
            newick += branch_length_stack.top();
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    newick += ';';

    invertTree(node);

    return newick;

}

void PangenomeMAT2::Tree::compressTreeParallel(PangenomeMAT2::Node* node, size_t level){
    node->level = level;

    if(node->children.size() == 0){
        return;
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, (int)node->children.size()), [&](tbb::blocked_range<int> r){
        for(int i = r.begin(); i < r.end(); i++){
            while(node->children[i]->children.size() == 1){
                mergeNodes(node->children[i], node->children[i]->children[0]);
            }
            // consolidate mutations of parent
            auto oldVector = node->children[i]->nucMutation;

            node->children[i]->nucMutation = consolidateNucMutations(node->children[i]->nucMutation);

            if(!debugSimilarity(oldVector, node->children[i]->nucMutation)){
                std::cout << "Inaccuracy observed in subtree extract." << std::endl;
            }

            compressTreeParallel(node->children[i], level + 1);
        }
    });
}

PangenomeMAT2::Node* subtreeExtractParallelHelper(PangenomeMAT2::Node* node, const tbb::concurrent_unordered_map< PangenomeMAT2::Node*, size_t >& ticks){
    if(ticks.find(node) == ticks.end()){
        return nullptr;
    }

    PangenomeMAT2::Node* newNode = new PangenomeMAT2::Node(node->identifier, node->branchLength);

    for(auto mutation: node->nucMutation){
        newNode->nucMutation.push_back(mutation);
    }

    for(auto mutation: node->blockMutation){
        newNode->blockMutation.push_back(mutation);
    }

    newNode->children.resize(node->children.size(), nullptr);

    tbb::parallel_for(tbb::blocked_range(0, (int)node->children.size()), [&](tbb::blocked_range<int> r){
        for(int i = r.begin(); i < r.end(); i++){
            PangenomeMAT2::Node* child = node->children[i];
            if(ticks.find(child) != ticks.end()){

                PangenomeMAT2::Node* newChild = subtreeExtractParallelHelper(child, ticks);

                newChild->parent = newNode;
                newNode->children[i] = newChild;
            }
        }
    });

    size_t i = 0, j = 0;
    while(j < newNode->children.size()){
        if(newNode->children[j] != nullptr){
            std::swap(newNode->children[i], newNode->children[j]);
            i++;
        }
        j++;
    }
    newNode->children.resize(i);
    
    return newNode;

}

PangenomeMAT2::Node* PangenomeMAT2::Tree::subtreeExtractParallel(std::vector< std::string > nodeIds){
    tbb::concurrent_vector< PangenomeMAT2::Node* > requiredNodes;

    std::atomic<bool> idDoesntExist = false;

    tbb::parallel_for_each(nodeIds.begin(), nodeIds.end(), [&]( std::string& id ) {
        if(allNodes.find(id) != allNodes.end()){
            requiredNodes.push_back(allNodes[id]);
        } else {
            idDoesntExist = true;
        }
    });

    if(idDoesntExist){
        std::cout << "Error: Some of the specified node identifiers don't exist!!!" << std::endl;
        return nullptr;
    }

    tbb::concurrent_unordered_map< PangenomeMAT2::Node*, size_t > ticks;

    tbb::parallel_for_each(requiredNodes.begin(), requiredNodes.end(), [&](PangenomeMAT2::Node*& node){
        Node* current = node;

        while(current != nullptr){
            ticks[current]++;
            current = current->parent;
        }
    });

    PangenomeMAT2::Node* newTreeRoot = subtreeExtractParallelHelper(root, ticks);

    compressTreeParallel(newTreeRoot, 1);

    return newTreeRoot;
}

void PangenomeMAT2::Tree::getNodesPreorder(PangenomeMAT2::Node* root, MATNew::tree& treeToWrite){
    
    MATNew::node n;

    for(size_t i = 0; i < root->nucMutation.size(); i++){
        const PangenomeMAT2::NucMut& mutation = root->nucMutation[i];

        MATNew::nucMut nm;
        nm.set_nucposition(mutation.nucPosition);
        if(mutation.nucGapPosition != -1){
            nm.set_nucgapposition(mutation.nucGapPosition);
            nm.set_nucgapexist(true);
        } else {
            nm.set_nucgapexist(false);
        }
        if(mutation.secondaryBlockId != -1){
            nm.set_blockid(((int64_t)mutation.primaryBlockId << 32) + mutation.secondaryBlockId);
            nm.set_blockgapexist(true);
        } else {
            nm.set_blockid(((int64_t)mutation.primaryBlockId << 32));
            nm.set_blockgapexist(false);
        }
        nm.set_mutinfo((((mutation.nucs) >> (24 - (mutation.mutInfo >> 4)*4)) << 8) + mutation.mutInfo);

        n.add_nucmutation();
        *n.mutable_nucmutation(i) = nm;
    }

    for(size_t i = 0; i < root->blockMutation.size(); i++){
        const PangenomeMAT2::BlockMut& mutation = root->blockMutation[i];

        MATNew::blockMut bm;
        if(mutation.secondaryBlockId != -1){
            bm.set_blockid(((int64_t)mutation.primaryBlockId << 32) + mutation.secondaryBlockId);
            bm.set_blockgapexist(true);
        } else {
            bm.set_blockid(((int64_t)mutation.primaryBlockId << 32));
            bm.set_blockgapexist(false);
        }
        bm.set_blockmutinfo(mutation.blockMutInfo);

        n.add_blockmutation();
        *n.mutable_blockmutation(i) = bm;
    }

    for(size_t i = 0; i < root->annotations.size(); i++){
        n.add_annotations(root->annotations[i]);
    }

    treeToWrite.add_nodes();
    *treeToWrite.mutable_nodes( treeToWrite.nodes_size() - 1 ) = n;

    for(auto child: root->children){
        getNodesPreorder(child, treeToWrite);
    }
}

void PangenomeMAT2::Tree::writeToFile(std::ofstream& fout, PangenomeMAT2::Node* node){
    if(node == nullptr){
        node = root;
    }

    MATNew::tree treeToWrite;
    getNodesPreorder(node, treeToWrite);

    std::string newick = getNewickString(node);

    treeToWrite.set_newick(newick);

    for(auto block: blocks){
        MATNew::block b;
        if(block.secondaryBlockId != -1){
            b.set_blockid(((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId);
            b.set_blockgapexist(true);
        } else {
            b.set_blockid(((int64_t)block.primaryBlockId << 32));
            b.set_blockgapexist(false);
        }
        b.set_chromosomename(block.chromosomeName);
        for(auto n: block.consensusSeq){
            b.add_consensusseq(n);
        }
        treeToWrite.add_blocks();
        *treeToWrite.mutable_blocks( treeToWrite.blocks_size() - 1 ) = b;
    }

    for(size_t i = 0; i < gaps.size(); i++){
        MATNew::gapList gl;
        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            gl.add_nucposition(gaps[i].nucPosition[j]);
            gl.add_nucgaplength(gaps[i].nucGapLength[j]);
        }
        if(gaps[i].secondaryBlockId != -1){
            gl.set_blockid(((int64_t)gaps[i].primaryBlockId << 32) + gaps[i].secondaryBlockId);
            gl.set_blockgapexist(true);
        } else {
            gl.set_blockid(((int64_t)gaps[i].primaryBlockId << 32));
            gl.set_blockgapexist(false);
        }
        treeToWrite.add_gaps();
        *treeToWrite.mutable_gaps( treeToWrite.gaps_size() - 1 ) = gl;
    }

    if (!treeToWrite.SerializeToOstream(&fout)) {
		std::cerr << "Failed to write output file." << std::endl;
    }

}

std::string PangenomeMAT2::Tree::getStringFromReference(std::string reference, bool aligned){

    Node* referenceNode = nullptr;
    
    for(auto u: allNodes){
     //   if(u.second->children.size() == 0 && u.first == reference){
        if(u.first == reference){
            referenceNode = u.second;
            break;
        }
    }

    if(referenceNode == nullptr){
        return "Error: Reference sequence with matching name not found!";
    }

    std::vector< PangenomeMAT2::Node* > path;
    Node* it = referenceNode;

    while(it != root){
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(root);

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(blocks.size() + 1, {false, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    }

    int32_t maxBlockId = 0;

    for(size_t i = 0; i < blocks.size(); i++){
        
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                
                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);

                if(secondaryBlockId != -1){
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }
            if(endFlag){
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1){
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }

    // Get all blocks on the path
    for(auto node = path.rbegin(); node != path.rend(); node++){
        for(auto mutation: (*node)->blockMutation){
            int primaryBlockId = mutation.primaryBlockId;
            int secondaryBlockId = mutation.secondaryBlockId;
            int type = (mutation.blockMutInfo);

            if(type == PangenomeMAT2::BlockMutationType::BI){
                if(secondaryBlockId != -1){
                    blockExists[primaryBlockId].second[secondaryBlockId] = true;
                } else {
                    blockExists[primaryBlockId].first = true;
                }
            } else {
                if(secondaryBlockId != -1){
                    blockExists[primaryBlockId].second[secondaryBlockId] = false;
                } else {
                    blockExists[primaryBlockId].first = false;
                }
            }
        }
    }

    // Apply nucleotide mutations
    for(auto node = path.rbegin(); node != path.rend(); node++){

        for(size_t i = 0; i < (*node)->nucMutation.size(); i++){

            int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = (*node)->nucMutation[i].secondaryBlockId;

            if(secondaryBlockId != -1){
                if(!blockExists[primaryBlockId].second[secondaryBlockId]){
                    continue;
                }
            } else {
                if(!blockExists[primaryBlockId].first){
                    continue;
                }
            }

            int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
            int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
            uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3){

                int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                if(type == PangenomeMAT2::NucMutationType::NS){
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                }
                else if(type == PangenomeMAT2::NucMutationType::NI){
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                }
                else if(type == PangenomeMAT2::NucMutationType::ND){
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            }

                        }
                    } else {
                        if(nucGapPosition != -1){
                            for(int j = 0; j < len; j++){
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            }
                        }
                    }
                }
            } 
            else {
                if(type == PangenomeMAT2::NucMutationType::NSNPS){
                    newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                }
                else if(type == PangenomeMAT2::NucMutationType::NSNPI){
                    newVal = PangenomeMAT2::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                }
                else if(type == PangenomeMAT2::NucMutationType::NSNPD){
                    if(secondaryBlockId != -1){
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        }
                    } else {
                        if(nucGapPosition != -1){
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    std::string sequenceString;
    for(size_t i = 0; i < sequence.size(); i++){
        for(size_t j = 0; j < sequence[i].second.size(); j++){
            if(blockExists[i].second[j]){
                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        if(sequence[i].second[j][k].second[w] == 'x' || sequence[i].second[j][k].second[w] == '-' ){
                            if(aligned){
                                sequenceString+='-';
                            }
                        } else {
                            sequenceString += sequence[i].second[j][k].second[w];
                        }
                    }
                    if(sequence[i].second[j][k].first == 'x' || sequence[i].second[j][k].first == '-'){
                        if(aligned){
                            sequenceString+='-';
                        }
                    } else {
                        sequenceString += sequence[i].second[j][k].first;
                    }
                }
            } else {
                if(aligned){
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                            sequenceString+='-';
                        }
                        sequenceString+='-';
                    }
                }
            }
        }
        if(blockExists[i].first){
            for(size_t j = 0; j < sequence[i].first.size(); j++){
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++){
                    if(sequence[i].first[j].second[k] == 'x' || sequence[i].first[j].second[k] == '-' ){
                        // This shouldn't be possible but I'm still keeping it since it doesn't hurt
                        if(aligned){
                            sequenceString += '-';
                        }
                    } else {
                        sequenceString += sequence[i].first[j].second[k];
                    }
                }
                if(sequence[i].first[j].first == 'x' || sequence[i].first[j].first == '-'){
                    if(aligned){
                        sequenceString += '-';
                    }
                } else {
                    sequenceString += sequence[i].first[j].first;
                }
            }
        } else {
            if(aligned){
                for(size_t j = 0; j < sequence[i].first.size(); j++){
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++){
                        sequenceString+='-';
                    }
                    sequenceString+='-';
                }
            }
        }
    }

    return sequenceString;

}

std::string PangenomeMAT2::stripString(std::string s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            return s.substr(i);
        }
    }
    return s;
}

std::string PangenomeMAT2::stripGaps(std::string sequenceString){
    std::string result;
    for(auto u: sequenceString){
        if(u != '-'){
            result+=u;
        }
    }
    return result;
}

bool PangenomeMAT2::Tree::verifyVCFFile(std::ifstream& fin){

    for(auto u: allNodes){
        if(u.second->children.size() == 0){
            fin.clear();
            fin.seekg(0);
            if(getSequenceFromVCF(u.first, fin) != stripGaps(getStringFromReference(u.first))){
                return false;
            }
        }
    }

    return true;

}

void PangenomeMAT2::Tree::setupGlobalCoordinates(){
    globalCoordinates.resize(blocks.size()+1);
    
    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        globalCoordinates[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
    }

    int32_t maxBlockId = 0;
    for(size_t i = 0; i < blocks.size(); i++){
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);
        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                
                if(secondaryBlockId != -1){
                    globalCoordinates[primaryBlockId].second[secondaryBlockId].push_back({0, {}});
                } else {
                    globalCoordinates[primaryBlockId].first.push_back({0, {}});
                }
            }

            if(endFlag){
                break;
            }
        }

        if(secondaryBlockId != -1){
            globalCoordinates[primaryBlockId].second[secondaryBlockId].push_back({0, {}});
        } else {
            globalCoordinates[primaryBlockId].first.push_back({0, {}});
        }
    }

    globalCoordinates.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                globalCoordinates[primaryBId].second[secondaryBId][pos].second.resize(len, 0);
            } else {
                globalCoordinates[primaryBId].first[pos].second.resize(len, 0);
            }
        }
    }

    // Assigning coordinates
    int ctr = 0;
    for(size_t i = 0; i < globalCoordinates.size(); i++){
        for(size_t j = 0; j < globalCoordinates[i].second.size(); j++){
            for(size_t k = 0; k < globalCoordinates[i].second[j].size(); k++){
                for(size_t w = 0; w < globalCoordinates[i].second[j][k].second.size(); w++){
                    globalCoordinates[i].second[j][k].second[w] = ctr;
                    ctr++;
                }
                globalCoordinates[i].second[j][k].first = ctr;
                ctr++;
            }
        }
        for(size_t j = 0; j < globalCoordinates[i].first.size(); j++){
            for(size_t k = 0; k < globalCoordinates[i].first[j].second.size(); k++){
                globalCoordinates[i].first[j].second[k] = ctr;
                ctr++;
            }
            globalCoordinates[i].first[j].first = ctr;
            ctr++;
        }
    }
}

size_t PangenomeMAT2::Tree::getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition, int nucGapPosition){
    if(secondaryBlockId == -1){
        if(nucGapPosition == -1){
            return globalCoordinates[primaryBlockId].first[nucPosition].first;
        }
        return globalCoordinates[primaryBlockId].first[nucPosition].second[nucGapPosition];
    } else {
        if(nucGapPosition == -1){
            return globalCoordinates[primaryBlockId].second[secondaryBlockId][nucPosition].first;
        }
        return globalCoordinates[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
    }
}

void PangenomeMAT2::Tree::vcfToFASTA(std::ifstream& fin, std::ofstream& fout){
    for(auto u: allNodes){
        if(u.second->children.size() == 0){
            fin.clear();
            fin.seekg(0);
            std::string sequenceString = getSequenceFromVCF(u.first, fin);
            fout << '>' << u.first << '\n';
            for(size_t i = 0; i < sequenceString.size(); i+=70){
                fout << sequenceString.substr(i,std::min(70, (int)sequenceString.length() - (int)i)) << '\n';
            }
        }
    }
}

AuxilaryMAT::Node* PangenomeMAT2::Tree::convertToAuxMatHelper(PangenomeMAT2::Node* root, std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
    std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > >& coordinates,\
    std::vector< std::pair< bool, std::vector< bool > > >& blockExists
){

    std::map< std::tuple< int, int, int, int >, char > mutations;
    std::vector< std::tuple< int32_t, int32_t, bool, bool > > blockMutationInfo;

    AuxilaryMAT::Node* newNode = new AuxilaryMAT::Node;
    newNode->identifier = root->identifier;

    for(auto mutation: root->blockMutation){
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = (mutation.blockMutInfo);

        if(type == 1){
            bool oldVal;
            if(secondaryBlockId != -1){
                if(!blockExists[primaryBlockId].second[secondaryBlockId]){
                    for(size_t i = 0; i < sequence[primaryBlockId].second[secondaryBlockId].size(); i++){
                        for(size_t j = 0; j < sequence[primaryBlockId].second[secondaryBlockId][i].second.size(); j++){
                            if(sequence[primaryBlockId].second[secondaryBlockId][i].second[j] != '-' && sequence[primaryBlockId].second[secondaryBlockId][i].second[j] != 'x'){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, j)] = sequence[primaryBlockId].second[secondaryBlockId][i].second[j];
                            }
                        }
                        if(sequence[primaryBlockId].second[secondaryBlockId][i].first != '-' && sequence[primaryBlockId].second[secondaryBlockId][i].first != 'x'){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, -1)] = sequence[primaryBlockId].second[secondaryBlockId][i].first;
                        }
                    }
                }
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = true;
            } else {
                if(!blockExists[primaryBlockId].first){
                    for(size_t i = 0; i < sequence[primaryBlockId].first.size(); i++){
                        for(size_t j = 0; j < sequence[primaryBlockId].first[i].second.size(); j++){
                            if(sequence[primaryBlockId].first[i].second[j] != '-' && sequence[primaryBlockId].first[i].second[j] != 'x'){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, j)] = sequence[primaryBlockId].first[i].second[j];
                            }
                        }
                        if(sequence[primaryBlockId].first[i].first != '-' && sequence[primaryBlockId].first[i].first != 'x'){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, -1)] = sequence[primaryBlockId].first[i].first;
                        }
                    }
                }
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = true;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, true) );
        } else {
            bool oldVal;
            if(secondaryBlockId != -1){
                if(blockExists[primaryBlockId].second[secondaryBlockId]){
                    for(size_t i = 0; i < sequence[primaryBlockId].second[secondaryBlockId].size(); i++){
                        for(size_t j = 0; j < sequence[primaryBlockId].second[secondaryBlockId][i].second.size(); j++){
                            if(sequence[primaryBlockId].second[secondaryBlockId][i].second[j] != '-' && sequence[primaryBlockId].second[secondaryBlockId][i].second[j] != 'x'){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, j)] = '-';
                            }
                        }
                        if(sequence[primaryBlockId].second[secondaryBlockId][i].first != '-' && sequence[primaryBlockId].second[secondaryBlockId][i].first != 'x'){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, -1)] = '-';
                        }
                    }
                }
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = false;
            } else {
                if(blockExists[primaryBlockId].first){
                    for(size_t i = 0; i < sequence[primaryBlockId].first.size(); i++){
                        for(size_t j = 0; j < sequence[primaryBlockId].first[i].second.size(); j++){
                            if(sequence[primaryBlockId].first[i].second[j] != '-' && sequence[primaryBlockId].first[i].second[j] != 'x'){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, j)] = '-';
                            }
                        }
                        if(sequence[primaryBlockId].first[i].first != '-' && sequence[primaryBlockId].first[i].first != 'x'){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, i, -1)] = '-';
                        }
                    }
                }
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = false;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, false) );
        }
    }

    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++){
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3){
            // Either S, I or D

            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(type == PangenomeMAT2::NucMutationType::NS){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition)] = newVal;
                            }
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition + j, nucGapPosition)] = newVal;
                            }
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NI){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition)] = newVal;
                            }
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition + j, nucGapPosition)] = newVal;
                            }
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::ND){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j)] = '-';
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                            if(blockExists[primaryBlockId].second[secondaryBlockId]){
                                mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition)] = '-';
                            }
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition, nucGapPosition+j)] = '-';
                            }
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                            if(blockExists[primaryBlockId].first){
                                mutations[std::make_tuple(primaryBlockId, -1, nucPosition + j, nucGapPosition)] = '-';
                            }
                        }
                    }
                }
            }
        } 
        else {
            if(type == PangenomeMAT2::NucMutationType::NSNPS){
                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPI){
                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = newVal;
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPD){
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = '-';
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                        if(blockExists[primaryBlockId].second[secondaryBlockId]){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = '-';
                        }
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = '-';
                        }
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                        if(blockExists[primaryBlockId].first){
                            mutations[std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition)] = '-';
                        }
                    }
                }
            }
        }
    }

    for(auto mut: mutations){
        uint32_t coor;
        if(std::get<1>(mut.first) != -1){
            if(std::get<3>(mut.first) != -1){
                coor = coordinates[std::get<0>(mut.first)].second[std::get<1>(mut.first)][std::get<2>(mut.first)].second[std::get<3>(mut.first)];
            } else {
                coor = coordinates[std::get<0>(mut.first)].second[std::get<1>(mut.first)][std::get<2>(mut.first)].first;
            }
        } else {
            if(std::get<3>(mut.first) != -1){
                coor = coordinates[std::get<0>(mut.first)].first[std::get<2>(mut.first)].second[std::get<3>(mut.first)];
            } else {
                coor = coordinates[std::get<0>(mut.first)].first[std::get<2>(mut.first)].first;
            }
        }
        newNode->substitutions.push_back({coor, mut.second});
    }

    for(auto child: root->children){
        newNode->children.push_back(convertToAuxMatHelper(child, sequence, coordinates, blockExists));
    }

    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
        }
    }

    // Undo nuc mutations
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }

    return newNode;
}

AuxilaryMAT::Tree* PangenomeMAT2::Tree::convertToAuxMat(){
    AuxilaryMAT::Tree* auxTree = new AuxilaryMAT::Tree;

    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(blocks.size() + 1);
    std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > > coordinates(blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(blocks.size() + 1, {false, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        coordinates[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    }

    // Creating blocks
    int32_t maxBlockId = 0;

    for(size_t i = 0; i < blocks.size(); i++){
        
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                
                if(secondaryBlockId != -1){
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                    coordinates[primaryBlockId].second[secondaryBlockId].push_back({0, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                    coordinates[primaryBlockId].first.push_back({0, {}});
                }
            }

            if(endFlag){
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1){
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
            coordinates[primaryBlockId].second[secondaryBlockId].push_back({0, {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
            coordinates[primaryBlockId].first.push_back({0, {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    coordinates.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
                coordinates[primaryBId].second[secondaryBId][pos].second.resize(len, 0);
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
                coordinates[primaryBId].first[pos].second.resize(len, 0);
            }
        }
    }

    // Assigning actual coordinates
    int ctr = 0;
    for(size_t i = 0; i < coordinates.size(); i++){
        for(size_t j = 0; j < coordinates[i].second.size(); j++){
            for(size_t k = 0; k < coordinates[i].second[j].size(); k++){
                for(size_t w = 0; w < coordinates[i].second[j][k].second.size(); w++){
                    coordinates[i].second[j][k].second[w] = ctr;
                    ctr++;
                }
                coordinates[i].second[j][k].first = ctr;
                ctr++;
            }
        }
        for(size_t j = 0; j < coordinates[i].first.size(); j++){
            for(size_t k = 0; k < coordinates[i].first[j].second.size(); k++){
                coordinates[i].first[j].second[k] = ctr;
                ctr++;
            }
            coordinates[i].first[j].first = ctr;
            ctr++;
        }
    }
    // Update length of consensus sequence
    auxTree->consensusSeqLength = ctr;
    auxTree->root = convertToAuxMatHelper(root, sequence, coordinates, blockExists);

    return auxTree;

}

std::string PangenomeMAT2::Tree::getSequenceFromVCF(std::string sequenceId, std::ifstream& fin){
    std::string line;

    // get reference line
    for(int i = 0; i < 4; i++){
        std::getline(fin, line);
    }

    if(line.substr(0,12) != "##reference="){
        std::cout << "Incorrect line format: " << line << std::endl;
        return "";
    }

    std::string referenceSequenceId = line.substr(12);

    std::string referenceSequence = stripGaps(getStringFromReference(referenceSequenceId));

    if(sequenceId == referenceSequenceId){
        return referenceSequence;
    }

    // column headers
    std::getline(fin, line);

    std::vector< std::string > columnWords;
    std::string word;

    for(size_t i = 0; i < line.size(); i++){
        if(line[i] != ' ' && line[i]!='\t'){
            word += line[i];
        } else {
            if(word.length()){
                columnWords.push_back(word);
                word = "";
            }
        }
    }
    if(word.length()){
        columnWords.push_back(word);
    }

    int sequenceIndex = -1;

    for(size_t i = 9; i < columnWords.size(); i++){
        if(columnWords[i] == sequenceId){
            sequenceIndex = i;
            break;
        }
    }

    if(sequenceIndex == -1){
        std::cout << "sequence not found! "<< sequenceId << std::endl;
        return "";
    }

    // To account for insertions
    std::vector< std::pair< char, std::vector< char > > > alteredSequence;
    for(auto u: referenceSequence){
        alteredSequence.push_back({u, {}});
    }
    alteredSequence.push_back({'-', {}});

    while(getline(fin, line)){
        std::vector< std::string > words;
        std::string word;

        for(size_t i = 0; i < line.size(); i++){
            if(line[i] != ' '){
                word += line[i];
            } else {
                if(word.length()){
                    words.push_back(word);
                    word = "";
                }
            }
        }
        if(word.length()){
            words.push_back(word);
            word="";
        }

        int choice = std::stoll(words[sequenceIndex]);
        if(choice == 0){
            continue;
        }

        choice--;

        int position = std::stoll(words[1]);

        std::string ref = words[3];
        std::string altStrings = words[4];

        std::string currentAlt;
        std::vector< std::string > altChoices;

        for(auto u: altStrings){
            if(u != ','){
                currentAlt += u;
            } else {
                if(currentAlt.length()){
                    altChoices.push_back(currentAlt);
                    currentAlt = "";
                }
            }
        }

        if(currentAlt.length()){
            altChoices.push_back(currentAlt);
            currentAlt = "";
        }

        std::string alt = altChoices[choice];

        if(ref != "."){
            int len = ref.length();
            for(int i = position; i < position + len; i++){
                alteredSequence[i].first = '-';
            }
        }

        if(alt != "."){
            if(alt.length() && alteredSequence[position].second.size()){
                std::cout << "VCF Error: alternate sequence already exists at position " << position <<"!" << std::endl;
                std::cout << sequenceId << " " << referenceSequenceId << std::endl;
            }
            for(size_t i = 0; i < alt.length(); i++){
                alteredSequence[position].second.push_back(alt[i]);
            }
        }

    }

    std::string finalSequence;
    for(size_t i = 0; i < alteredSequence.size(); i++){
        for(size_t j = 0; j < alteredSequence[i].second.size();j++){
            if(alteredSequence[i].second[j] != '-'){
                finalSequence += alteredSequence[i].second[j];
            }
        }
        if(alteredSequence[i].first != '-'){
            finalSequence += alteredSequence[i].first;
        }

    }

    std::string alteredSequenceOriginal = stripGaps(getStringFromReference(sequenceId));

    return finalSequence;

}

void PangenomeMAT2::Tree::printFASTAParallel(std::ofstream& fout, bool aligned){
    std::mutex fastaMutex;
    size_t lineSize = 70;

    tbb::parallel_for_each(allNodes, [&](auto n){
        std::string sequence;
        sequence = getStringFromReference(n.first, aligned);
        
        fastaMutex.lock();
        fout << '>' << n.first << '\n';
        for(size_t i = 0; i < sequence.size(); i+=lineSize){
            fout << sequence.substr(i, std::min(lineSize, sequence.size() - i)) << '\n';
        }
        fastaMutex.unlock();
    });
}

void processNode(std::string referenceSequence, std::string altSequence, std::string nid, std::string parentId, std::ofstream &fout, std::unordered_map<std::string, PangenomeMAT2::Node*> &allNodes) {

    size_t recordID = 0;

    std::mutex vcfMapMutex;
    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    if(altSequence.length() != referenceSequence.length()){
        //std::cout << altSequence << "\n=====\n" << referenceSequence << "\n";
        std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
        return;
    }
  // std::cerr << ">ref " << parentId < "\n" << referenceSequence << "\n>alt " << nid << "\n" << altSequence << std::endl;


    std::string currentRefString, currentAltString;
    int currentCoordinate = 0;

    int diffStart = 0;

    for(size_t i = 0; i < referenceSequence.length(); i++){

        if(referenceSequence[i] == '-' && altSequence[i] == '-'){
            continue;
        } else if(referenceSequence[i] != '-' && altSequence[i] == '-'){
            if(currentRefString == "" && currentAltString == ""){
                diffStart = currentCoordinate;
            }

            currentRefString += referenceSequence[i];
        } else if(referenceSequence[i] == '-' && altSequence[i] != '-'){
            if(currentRefString == "" && currentAltString == ""){
                diffStart = currentCoordinate;
            }

            currentAltString += altSequence[i];
        } else if(referenceSequence[i] != altSequence[i]){
            if(currentRefString == "" && currentAltString == ""){
                diffStart = currentCoordinate;
            }
            if(currentRefString == currentAltString){
                currentRefString = "";
                currentAltString = "";
                diffStart = currentCoordinate;
            }
            currentRefString += referenceSequence[i];
            currentAltString += altSequence[i];
        } else if(referenceSequence[i] == altSequence[i]){
            if(currentRefString == currentAltString){
                // Reset
                diffStart = currentCoordinate;
                currentRefString = "";
                currentRefString += referenceSequence[i];
                currentAltString = currentRefString;
            } else {
                // Create VCF record at position i
                if(currentRefString == ""){
                    currentRefString += referenceSequence[i];
                    currentAltString += altSequence[i];
                    diffStart = currentCoordinate;
                    vcfMapMutex.lock();
                    vcfMap[diffStart][currentRefString][currentAltString].push_back(nid);
                    vcfMapMutex.unlock();
                    diffStart = currentCoordinate+1;
                    currentRefString = "";
                    currentAltString = "";
                } else {
                    vcfMapMutex.lock();
                    vcfMap[diffStart][currentRefString][currentAltString].push_back(nid);
                    vcfMapMutex.unlock();

                    // Reset
                    diffStart = currentCoordinate;
                    currentRefString = "";
                    currentRefString += referenceSequence[i];
                    currentAltString = currentRefString;
                }
            }
        }

        if(referenceSequence[i] != '-'){
            currentCoordinate++;
        }
    }

    if(currentRefString != currentAltString){
        vcfMapMutex.lock();
        vcfMap[diffStart][currentRefString][currentAltString].push_back(nid);
        vcfMapMutex.unlock();

        // Reset
        diffStart = referenceSequence.size();
        currentRefString = "";
        currentAltString = currentRefString;
    }
    



    fout << "##fileformat=VCFv" << VCF_VERSION << '\n';
    fout << "##fileDate=" << PangenomeMAT2::getDate() << '\n';
    fout << "##source=PanMATv" << PMAT_VERSION << '\n';
    fout << "##reference=" << parentId << '\n';
    fout << "#CHROM\t" << "POS\t" << "ID\t" << "REF\t" << "ALT\t" << "QUAL\t" << "FILTER\t" << "INFO\t" << "FORMAT\t";
    
    fout << nid + "\t";
    fout << '\n';

    std::map< std::string, size_t > sequenceIds;
    sequenceIds[nid] = 0;


    for(auto u: vcfMap){
        for(auto v: u.second){
            if(v.first == ""){
                fout << ".\t" << u.first << "\t" << recordID++ << "\t" << ".\t";
            } else {
                fout << ".\t" << u.first << "\t" << recordID++ << "\t" << v.first << "\t";
            }
            
            std::map< std::string, size_t > tempSequenceIds = sequenceIds;

            int ctr = 1;
            std::string altStrings;

            for(auto w: v.second){
                altStrings += (w.first == "" ? ".": w.first);
                altStrings += ",";
                for(auto uu: w.second){
                    tempSequenceIds[uu] = ctr;
                }
                ctr++;
            }

            altStrings.pop_back();

            fout << altStrings << "\t.\t.\t.\t.\t";

            for(auto w: tempSequenceIds){
                fout << w.second << "\t";
            }

            fout << '\n';
        }
    }
}

std::pair<int32_t, int32_t> getRecomputePositions(std::pair<int32_t, int32_t> p, std::string &gappedSequence, int32_t k) {
    int32_t mutPos = p.first;
    int32_t mutLen = p.second;

    int32_t i = mutPos;
    int32_t curr = i;
    while (i > std::max((int32_t) 0, (int32_t) (mutPos - k + 1))) {

        if (gappedSequence[curr] != '-') {
            i--;
        }
        curr--;
    }
    int32_t start = curr;

    i = mutPos + mutLen;
    curr = i;
    while (i < (mutPos+mutLen) + k + 1) {
        if (gappedSequence[curr] != '-') {
            i++;
        }
        curr++;
    }
    int32_t stop = curr;
    
    return std::make_pair(start, stop);

    
}

int32_t alignedDist(int32_t pos, int32_t k, std::string &gappedSequence) {
    int32_t i = pos;
    int32_t curr = i;
    while (i < pos + k - 1 && curr < gappedSequence.size()) {

        if (gappedSequence[curr] != '-') {
            i++;
        }
        curr++;
    }
    return curr - pos;
}


// bool compareTuples(const std::pair<int32_t, int32_t>& a, const std::pair<int32_t, int32_t>& b) {
//     return a.first < b.first;
// }

// std::vector<std::pair<int32_t, int32_t>> mergeOverlappingTuples(std::vector<std::pair<int32_t, int32_t>>& tuples) {
//     std::sort(tuples.begin(), tuples.end(), compareTuples);

//     std::vector<std::pair<int32_t, int32_t>> mergedTuples;
//     mergedTuples.push_back(tuples[0]);

//     for (int32_t i = 1; i < tuples.size(); i++) {
//         int32_t currentStart = tuples[i].first;
//         int32_t currentEnd = tuples[i].second;
//         int32_t mergedEnd = mergedTuples.back().second;

//         if (currentStart > mergedEnd) {
//             mergedTuples.push_back(tuples[i]);
//         }
//         else {
//             mergedTuples.back().second = std::max(currentEnd, mergedEnd);
//         }
//     }

//     return mergedTuples;
// }

void discardSyncmers(std::vector<kmer_t> &inSyncmers, const std::vector<std::pair<int32_t, int32_t>>& B, std::string &gappedSequence, std::unordered_map<std::string, kmer_t> &to_insert, std::unordered_map<std::string, bool> &variable_syncmers, seedIndex &index, std::string nid) {

    for (int32_t i = inSyncmers.size() - 1; i >= 0; i--) {
        const kmer_t s = inSyncmers[i];
        for (const auto& b : B) {

            if (s.pos >= b.first && s.pos+alignedDist(s.pos, 15, gappedSequence) < b.second) {

                variable_syncmers[s.seq] = true; // track variant sites
                auto it = to_insert.find(s.seq);
          
                if (it == to_insert.end()) {
                    index.deletions[nid].push_back(kmer_t{s.seq, s.pos, i}); 
                    inSyncmers.erase(inSyncmers.begin() + i);
                } else {
                    to_insert.erase(it->first);
                }
                break;
            }
        }
    }
}
void PangenomeMAT2::Tree::indexSyncmersHelper(PangenomeMAT2::Node* root,\
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
    std::vector< std::pair< bool, std::vector< bool > > >& blockExists, seedIndex &index, std::unordered_map<std::string, PangenomeMAT2::Node*> &allNodes, std::vector<kmer_t> &syncmers,\
    std::unordered_map<std::string, int32_t> &counts, std::unordered_map<std::string, bool> &variable_syncmers){


    //std::cout << root->identifier << "\n";
    
    std::vector< std::tuple< int32_t, int32_t, bool, bool > > blockMutationInfo;

    for(const auto &mutation: root->blockMutation){
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = (mutation.blockMutInfo);

        if(type == 1){
            bool oldVal;
            if(secondaryBlockId != -1){
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = true;
            } else {
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = true;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, true) );
        } else {
            bool oldVal;
            if(secondaryBlockId != -1){
                oldVal = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = false;
            } else {
                oldVal = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = false;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldVal, false) );
        }

    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Apply mutations to "sequence"
    std::vector<std::pair<size_t, size_t>> mutPositions;

    for(size_t i = 0; i < root->nucMutation.size(); i++){
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';
        size_t globalCoord = getGlobalCoordinate(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition);
        if(type < 3){
            // Either S, I or D

            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(type == PangenomeMAT2::NucMutationType::NS){
                mutPositions.push_back(std::make_pair(globalCoord, len));
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
        //                    std::cout << "MUT Substitution " << oldVal << globalCoord << newVal << "\n";
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            //                std::cout << "MUT Substitution " << oldVal << globalCoord << newVal << "\n";

                        }

                    }
                    
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
             //               std::cout << "MUT Substitution " << oldVal << globalCoord << newVal << "\n";

                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
               //             std::cout << "MUT Substitution " << oldVal << globalCoord << newVal << "\n";

                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NI){
                mutPositions.push_back(std::make_pair(globalCoord, len));

                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                 //           std::cout << "MUT Insertion " << oldVal << globalCoord << newVal << "\n";

                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                //            std::cout << "MUT Insertion " << oldVal << globalCoord << newVal << "\n";

                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                //            std::cout << "MUT Insertion " << oldVal << globalCoord << newVal << "\n";
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                //            std::cout << "MUT Insertion " << oldVal << globalCoord << newVal << "\n";
                        }
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::ND){
                mutPositions.push_back(std::make_pair(globalCoord, 0));
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                       //     std::cout << "MUT deletion " << oldVal << globalCoord << "\n";
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                      //      std::cout << "MUT deletion " << oldVal << globalCoord << "\n";
                        }

                    }
                } else {
                    if(nucGapPosition != -1){
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                      //      std::cout << "MUT deletion " << oldVal << globalCoord << "\n";
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                      //      std::cout << "MUT deletion " << oldVal << globalCoord << "\n";
                        }
                    }
                }
            }
        } 
        else {
            if(type == PangenomeMAT2::NucMutationType::NSNPS){
                mutPositions.push_back(std::make_pair(globalCoord, 0));

                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                  //      std::cout << "MUT SNP S " << oldVal << globalCoord << newVal << "\n";

                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                  //      std::cout << "MUT SNP S " << oldVal << globalCoord << newVal << "\n";

                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                  //      std::cout << "MUT SNP S " << oldVal << globalCoord << newVal << "\n";
                    } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                   //     std::cout << "MUT SNP S " << oldVal << globalCoord << newVal << "\n";
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPI){
            mutPositions.push_back(std::make_pair(globalCoord, 1));

                newVal = PangenomeMAT2::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    //    std::cout << "MUT SNP I " << oldVal << globalCoord << newVal << "\n";
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    //    std::cout << "MUT SNP I " << oldVal << globalCoord << newVal << "\n";
                    }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                    //     std::cout << "MUT SNP I " << oldVal << globalCoord << newVal << "\n";
                   } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));   
                    //    std::cout << "MUT SNP I " << oldVal << globalCoord << newVal << "\n";
                    }
                }
            }
            else if(type == PangenomeMAT2::NucMutationType::NSNPD){
            mutPositions.push_back(std::make_pair(globalCoord, -1));

                if(secondaryBlockId != -1){
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    //      std::cout << "MUT SNP D " << oldVal << globalCoord << newVal << "\n";
                  } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                     //      std::cout << "MUT SNP D " << oldVal << globalCoord << newVal << "\n";
                   }
                } else {
                    if(nucGapPosition != -1){
                        char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                      //      std::cout << "MUT SNP D " << oldVal << globalCoord << newVal << "\n";
                  } else {
                        char oldVal = sequence[primaryBlockId].first[nucPosition].first;
                        sequence[primaryBlockId].first[nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                     //        std::cout << "MUT SNP D " << oldVal << globalCoord << newVal << "\n";
                 }
                }
            }
        }
    }

    // Aligned sequence of the current node
    std::string currNodeSequence = getSequence(sequence, blockExists, true);
//    std::cout << currNodeSequence << "\n";

    std::vector<std::pair<int32_t, int32_t>> recompute;
    for (const std::pair<int32_t, int32_t> &p : mutPositions) {
        auto r = getRecomputePositions(p, currNodeSequence, 15);
        recompute.push_back(r);
    }

    std::unordered_map<std::string, kmer_t> to_insert;
    //std::cout << "\n***** " << root->identifier << "\n";

    for (auto &range : recompute) {
        //std::cout << "redo range: " << range.first << ", " << range.second << "\n";
        int32_t seqLen = currNodeSequence.size();
        std::string redo = currNodeSequence.substr(std::max(0, range.first), std::min(seqLen - range.first, range.second - range.first)); 
        auto redone = syncmerize(redo, 15, 8, false, true, std::max(0,range.first));
        for (const kmer_t &syncmer : redone) {
            //std::cout << " => redone: " << syncmer.seq << " " << syncmer.pos << " idx: " << syncmer.idx << "\n";
            to_insert[syncmer.seq] = syncmer;
        }
    }
    discardSyncmers(syncmers, recompute, currNodeSequence, to_insert, variable_syncmers, index, root->identifier); //modifies mutated_syncmers
    for (auto &s : to_insert) {
        syncmers.push_back(s.second);
        index.insertions[root->identifier].push_back(s.second);
    }
    for (auto d : index.deletions[root->identifier]) {
        //std::cout << "(x) " << d.seq << " " << d.pos << " idx: " << d.idx << "\n";
    }
    for (auto d : index.insertions[root->identifier]) {
        //std::cout << "(+) " << d.seq << " " << d.pos << " idx: " << d.idx << "\n";
    }
    std::cout << "\n:" << root->identifier << "\n";
    for (kmer_t s : syncmers) {
        std::cout << ">" << s.seq << "\n" << s.seq << "\n";
    }
    //std::cout << "\n";
   for(PangenomeMAT2::Node* child: root->children){
        indexSyncmersHelper(child, sequence, blockExists, index, allNodes, syncmers, counts, variable_syncmers);
    }
    
    //std::cout << "undo by deleting last " << index.insertions[root->identifier].size() << "\n";
    //std::cout << "size before " << syncmers.size() << "\n";
    syncmers.erase(syncmers.end() - index.insertions[root->identifier].size(), syncmers.end());
    //std::cout << "size after " << syncmers.size() << "\n";

    for (int32_t i = index.deletions[root->identifier].size() - 1; i >= 0; i--) {
        //std::cout << "undo by inserting " << index.deletions[root->identifier][i].seq << " at " << index.deletions[root->identifier][i].idx << "\n";
        syncmers.insert(syncmers.begin() + index.deletions[root->identifier][i].idx, index.deletions[root->identifier][i]);
    }

    //std::cout << "after undo:\n";
    for (kmer_t s : syncmers) {
        //std::cout << s.seq << ",";
    }
    //std::cout << "\n";

    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
        }
    }

    // Undo nuc mutations
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++){
        auto mutation = *it;
        if(std::get<1>(mutation) != -1){
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1){
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }
}

std::string reverse_complement(std::string dna_sequence) {
    std::string complement = "";
    for (char c : dna_sequence) {
        switch (c) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            default: complement += c; break;
        }
    }
    std::reverse(complement.begin(), complement.end());
    return complement;
}

void PangenomeMAT2::Tree::indexSyncmers(std::ofstream& fout){

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(blocks.size() + 1, {false, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    }
    
    int32_t maxBlockId = 0;

    for(size_t i = 0; i < blocks.size(); i++){
        
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                
                if(secondaryBlockId != -1){
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                    blockExists[primaryBlockId].second[secondaryBlockId] = true;
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                    blockExists[primaryBlockId].first = true;

                }
            }

            if(endFlag){
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1){
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }
    
    std::string consensusSequence = getSequence(sequence, blockExists, true);

    std::vector<kmer_t> initialSyncmers;
    std::unordered_map<std::string, int> counts;


    initialSyncmers = syncmerize(consensusSequence, 15, 8, false, true, 0);


    seedIndex index;
    std::unordered_map<std::string, bool> variable_syncmers;
    indexSyncmersHelper(root, sequence, blockExists, index, allNodes, initialSyncmers, counts, variable_syncmers);
    std::unordered_map<std::string, bool> invariants;
    std::vector<kmer_t> initialShaved;
    for (const kmer_t &syncmer : initialSyncmers) {
        if (variable_syncmers.find(syncmer.seq) == variable_syncmers.end()) {
            invariants[syncmer.seq] = true;
        } else {
            initialShaved.push_back(kmer_t{syncmer.seq, syncmer.pos, -1});
        }
    }
    std::cout << "total consensus seeds: " << initialSyncmers.size() << "\n";
    std::cout << "# invariant seeds: " << invariants.size() << "\n";
    index.rootSeeds = initialSyncmers;
    seedIndex shaved = shaveIndex(index, invariants, initialShaved);
    shaved.rootSeeds = initialShaved;
    writeIndex(fout, shaved);
}

template <typename INT, typename T>
void removeIndices(std::vector<T>& v, std::stack<INT>& rm )
{
    if (rm.size() < 1) {
        return;
    }
  int32_t rmVal = rm.top();
  rm.pop();
  v.erase(
    std::remove_if(std::begin(v), std::end(v), [&](T& elem)
    {
        if (!rm.empty() && &elem - &v[0] == rmVal)
        {
           rmVal = rm.top();
           rm.pop();
           return true;
        }
        return false;
    }),
    std::end(v)
  );
}
std::vector<kmer_t> syncmersFromFastq(std::string fastqPath,  std::vector<read_t> &reads) {
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    seq = kseq_init(fileno(fp));
    std::vector<std::string> input;
    std::vector<std::string> input_names;
    
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        std::string this_seq  = seq->seq.s;
        std::string this_name = seq->name.s;

        input.push_back(this_seq);
        input_names.push_back(this_name);
    }
    float est_coverage = 0; //TODO change this to 1
    int k = 15;
    int s = 8;
    bool open = false;
    
    std::unordered_set<kmer_t, KHash> syncmers;
    std::unordered_map<std::string, int> counts;
    std::unordered_map<std::string, int> counts_rc;

    //std::cerr << "length: " << input.size() << "\n";
    reads.resize(input.size());

    for (int i = 0; i < input.size(); i++) {        
        read_t this_read;
        std::string seq = input[i];
        std::string name = input_names[i];
        
        this_read.seq = seq;
        this_read.name = name;


        std::string rc = reverse_complement(seq);
        std::vector<kmer_t> these = syncmerize(seq, k, s, false, false, 0);
        std::vector<kmer_t> these_rc = syncmerize(rc, k, s, false, false, 0);
        
        
        for (const auto &m : these) {
            if (counts.find(m.seq) == counts.end()) {
                counts[m.seq] = 1;
            } else {
                counts[m.seq] += 1;
            }
            if (counts[m.seq] > est_coverage) {
                syncmers.insert(m);

                //std::cerr << syncmers.begin() << "\n";
 //               m.pos = m.pos + k - 1;
//                m.reversed = false;
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + k - 1, -1, 0, false});
                //this_read.read_coord.push_back(m.pos + k - 1);
                //this_read.reversed.push_back(false);
            }
        }
        
        for (const auto &m : these_rc) {
            if (counts_rc.find(m.seq) == counts_rc.end()) {
                counts_rc[m.seq] = 1;
            } else {
                counts_rc[m.seq] += 1;
            }
            if (counts_rc[m.seq] > est_coverage) {
                syncmers.insert(m);
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + k - 1, -1, 0, true});
            }
        }
        reads[i] = this_read;
    }
    std::vector<kmer_t> v;
    v.insert(v.end(), syncmers.begin(), syncmers.end());

    return v;
}



void updateJaccard(dynamicJaccard &dj, std::unordered_map<std::string, bool> &readSyncmers, std::vector<kmer_t> &deletedSyncmers, std::vector<kmer_t> &insertedSyncmers) {
    for (const kmer_t &syncmer : deletedSyncmers) {
        if (readSyncmers.find(syncmer.seq) != readSyncmers.end()) {
            dj.intersectionSize -= 1;
        } else {
            dj.unionSize -= 1;
        }
    }
    for (const kmer_t &syncmer : insertedSyncmers) {
        if (readSyncmers.find(syncmer.seq) != readSyncmers.end()) {
            dj.intersectionSize += 1;
        } else {
            dj.unionSize += 1;
        }
    }

    dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
}

void PangenomeMAT2::Tree::placeDFS(Node *currNode, std::vector<kmer_t> &currNodeSyncmers, std::unordered_map<std::string, bool> &querySyncmers, seedIndex &index, dynamicJaccard dj, std::unordered_map<std::string, float> &scores) {
    
    std::stack<int32_t> delIndices;
    for (const kmer_t &s : index.deletions[currNode->identifier]) {
        delIndices.push(s.idx);
    }

    removeIndices(currNodeSyncmers, delIndices);

    for (const kmer_t &s : index.insertions[currNode->identifier]) {
        currNodeSyncmers.push_back(s);
    }

    updateJaccard(dj, querySyncmers, index.deletions[currNode->identifier], index.insertions[currNode->identifier]);

    scores[currNode->identifier] = dj.jaccardIndex;

    for (Node *child : currNode->children) {
        placeDFS(child, currNodeSyncmers, querySyncmers, index, dj, scores);
    }

    currNodeSyncmers.erase(currNodeSyncmers.end() - index.insertions[currNode->identifier].size(), currNodeSyncmers.end());

    for (int32_t i = index.deletions[currNode->identifier].size() - 1; i >= 0; i--) {
        currNodeSyncmers.insert(currNodeSyncmers.begin() + index.deletions[currNode->identifier][i].idx, index.deletions[currNode->identifier][i]);
    }

}

struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

void PangenomeMAT2::Tree::shaveDFS(Node *currNode, std::vector<kmer_t> &currNodeSyncmers, seedIndex &index_orig, std::unordered_map<std::string, std::vector<kmer_t>> &deletions, std::unordered_map<std::string, bool> &invariants) {
    
    std::vector<int32_t> di;
    for (kmer_t d : index_orig.deletions[currNode->identifier]) {
        auto r = std::find(currNodeSyncmers.begin(), currNodeSyncmers.end(), d);
        if (r != currNodeSyncmers.end()) {
            deletions[currNode->identifier].push_back(kmer_t{r->seq, r->pos, static_cast<int32_t>(r - currNodeSyncmers.begin())});
            di.push_back(static_cast<int32_t>(r - currNodeSyncmers.begin()));
        }
    }
    std::sort(deletions[currNode->identifier].begin(), deletions[currNode->identifier].end(),
         [](const kmer_t &x, const kmer_t &y){ return (x.idx > y.idx); });

    std::sort(di.begin(), di.end());

    std::stack<int32_t> delIndices;
    for (auto i = di.rbegin(); i != di.rend(); i++) {
        delIndices.push(*i);
    }
    
    removeIndices(currNodeSyncmers, delIndices);

    for (const kmer_t &s : index_orig.insertions[currNode->identifier]) {
        currNodeSyncmers.push_back(s);
    }
    
    for (Node *child : currNode->children) {
        shaveDFS(child, currNodeSyncmers, index_orig, deletions, invariants);
    }

    currNodeSyncmers.erase(currNodeSyncmers.end() - index_orig.insertions[currNode->identifier].size(), currNodeSyncmers.end());
    
    for (int32_t i = deletions[currNode->identifier].size() - 1; i >= 0; i--) {
        currNodeSyncmers.insert(currNodeSyncmers.begin() + deletions[currNode->identifier][i].idx, deletions[currNode->identifier][i]);
    }
    
}



seedIndex PangenomeMAT2::Tree::shaveIndex(seedIndex &index, std::unordered_map<std::string, bool> &invariants, std::vector<kmer_t> &initialSyncmers) {
    seedIndex shaved;
    shaved.insertions = index.insertions;
    std::unordered_map<std::string, std::vector<kmer_t>> deletions;
    shaveDFS(root, initialSyncmers, index, deletions, invariants);
    shaved.deletions = deletions;
    return shaved;
}



void PangenomeMAT2::Tree::writeIndexDFS(Node *currNode, seedIndex &index, std::stringstream &ss, std::vector<kmer_t> &seeds) {
    ss << currNode->identifier << "\t";
  

    std::stack<int32_t> delIndices;
    for (const kmer_t &s : index.deletions[currNode->identifier]) {
        ss << s.idx << " ";
        delIndices.push(s.idx);
    }

    removeIndices(seeds, delIndices);

    ss << "^";
    for (const kmer_t &s : index.insertions[currNode->identifier]) {
        seeds.push_back(s);
        ss << s.seq << " ";
    }

    ss << "\n";

    for (Node *child : currNode->children) {
        writeIndexDFS(child, index, ss, seeds);
    }

    seeds.erase(seeds.end() - index.insertions[currNode->identifier].size(), seeds.end());
    
    for (int32_t i = index.deletions[currNode->identifier].size() - 1; i >= 0; i--) {
        seeds.insert(seeds.begin() + index.deletions[currNode->identifier][i].idx, index.deletions[currNode->identifier][i]);
    }
    
}

void PangenomeMAT2::Tree::writeIndex(std::ofstream &fout, seedIndex &index) {
    for (const kmer_t &s : index.rootSeeds) {
        fout << s.seq << " ";
    }
    fout << "\n";
    std::stringstream ss;
    writeIndexDFS(root, index, ss, index.rootSeeds);
    fout << ss.str();
}


void PangenomeMAT2::Tree::loadIndex(std::ifstream &indexFile, seedIndex &index) {
    std::string rootSyncmersString;
    std::vector<kmer_t> rootSyncmers;
    std::getline(indexFile, rootSyncmersString);
    std::vector< std::string > rootSplt;
    PangenomeMAT2::stringSplit(rootSyncmersString, ' ', rootSplt);
    for (const std::string &s : rootSplt) {
        if (s.size() > 1) {
            rootSyncmers.push_back(kmer_t{s, 0});
        }
    }

    std::string line;
    while (getline(indexFile, line)) {
        if (line.size() < 2) {
            continue;
        }
        std::vector<std::string> spltNid;
        std::vector<std::string> spltParts;

        std::vector<std::string> spltDel = {};
        std::vector<std::string> spltIns = {};


        stringSplit(line, '\t', spltNid);
        std::string nodeId = spltNid[0];

        stringSplit(spltNid[1], '^', spltParts);
        stringSplit(spltParts[0], ' ', spltDel);
        if (spltParts.size() > 1) {
            stringSplit(spltParts[1], ' ', spltIns);
        }

        std::vector<kmer_t> deletions;
        for (const std::string &s : spltDel) {
            deletions.push_back(kmer_t{"", 0, std::atoi(s.c_str())});
        }
        std::vector<kmer_t> insertions;
        for (const std::string &s : spltIns) {
            if (s.size() > 1) {
                insertions.push_back(kmer_t{s, 0});
            }
        }
        
        index.rootSeeds = rootSyncmers;
        index.insertions[nodeId] = insertions;
        index.deletions[nodeId] = deletions;
  
    }
}

struct Counter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type&) { ++count; }
  int32_t count = 0;
};

template<typename T1, typename T2>
int32_t intersection_size(const T1& s1, const T2& s2)
{
  Counter c;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
  return c.count;
}

extern "C" {
    void align_reads(const char *reference, int n_reads, const char **reads, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions);
}

void PangenomeMAT2::Tree::placeSample(std::string fastqPath, seedIndex &index){
    
    std::vector<read_t> reads;

    std::vector<kmer_t> readSyncmers = syncmersFromFastq(fastqPath, reads);

    int k = 15;
    int s = 8;

    std::cerr << "\n";
    std::cerr << "Placing sample...\n";


    struct dynamicJaccard dj;
    dj.intersectionSize = intersection_size(index.rootSeeds, readSyncmers);
    dj.unionSize = index.rootSeeds.size() + readSyncmers.size() - dj.intersectionSize;
    dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
    
    std::unordered_map<std::string, float> scores;
    std::unordered_map<std::string, bool> readSyncmersMap;
    for (const auto &k : readSyncmers) {
        readSyncmersMap[k.seq] = true;
    }

    placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);
    std::vector<std::pair<std::string, float>> v;
    for ( const auto &p : scores ) {
        v.push_back(std::make_pair(p.first, p.second));
    } 
    std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
        return left.second > right.second;
    });

    std::string best_match = v[0].first;
    for (const auto &s : v) {
        std::cerr << s.first << ": " << s.second << "\n";
    }
    //First in this s vector is the highest scoring node
    //   that will be reference

    std::string ref_seq = getStringFromReference(best_match, false);


    std::vector<kmer_t> ref_syncmers = syncmerize(ref_seq, k, s, false, false, 0);
    

    
    //Finding syncmer matches
    for(int r = 0; r < reads.size() ; r++){

        auto readSyncmer = reads[r].kmers.begin();
        auto refSyncmer = ref_syncmers.begin();

        std::vector<kmer_t> matchingSyncmers;

        while(readSyncmer != reads[r].kmers.end() && refSyncmer != ref_syncmers.end()) {
            if(*readSyncmer < *refSyncmer){
                readSyncmer++;
            }else if (*refSyncmer < *readSyncmer){
                refSyncmer++;
            }else{

                kmer_t matched;
                matched.seq  = (*readSyncmer).seq;
                matched.pos  = (*readSyncmer).pos;
                matched.pos2 = (*refSyncmer).pos + k - 1;
                matched.reversed  = (*readSyncmer).reversed;
                matchingSyncmers.push_back(matched);
                
                readSyncmer++;
                refSyncmer++;
            }
        }

        reads[r].kmers = matchingSyncmers;
    }

    std::cout << "\n" << ref_seq << std::endl;


    //Figure out a better way to do this
    //C++ to C interface will have to be cleaned up
    const char *reference = ref_seq.c_str();
    int n_reads = reads.size();
    const char **read_strings = (const char **)malloc(n_reads*sizeof(char *));
    int *r_lens         = (int *)malloc(n_reads*sizeof(int));
    int *seed_counts    = (int *)malloc(n_reads*sizeof(int));

    uint8_t **reversed  = (uint8_t **)malloc(n_reads*sizeof(uint8_t *));
    int **ref_positions = (int **)malloc(n_reads*sizeof(int *));
    int **qry_positions = (int **)malloc(n_reads*sizeof(int *));

    for(int i = 0; i < n_reads; i++) {
        int n_seeds = reads[i].kmers.size();
        seed_counts[i] = n_seeds;
        read_strings[i] = reads[i].seq.c_str();
        r_lens[i] = reads[i].seq.length();

        uint8_t *reversed_array = (uint8_t *)malloc(n_seeds*sizeof(uint8_t));
        int *ref_pos_array = (int *)malloc(n_seeds*sizeof(int));
        int *qry_pos_array = (int *)malloc(n_seeds*sizeof(int));

        int j = 0;
        for (auto it = reads[i].kmers.begin(); it != reads[i].kmers.end(); it++) {
            reversed_array[j] = (*it).reversed;
            qry_pos_array[j] = (*it).pos;
            ref_pos_array[j] = (*it).pos2;
            j++;
        }

        reversed[i]      = reversed_array;
        ref_positions[i] = ref_pos_array;
        qry_positions[i] = qry_pos_array;
    }
    
    std::cout << "\n@SQ	SN:reference	LN:" << ref_seq.length() << std::endl;
    align_reads(reference, n_reads, read_strings, r_lens, seed_counts, reversed, ref_positions, qry_positions);


    for(int i = 0; i < n_reads; i++) {
        free(reversed[i]);
        free(ref_positions[i]);
        free(qry_positions[i]);
    }
    free(qry_positions);
    free(ref_positions);
    free(reversed);
    free(seed_counts);
    free(r_lens);
    free(read_strings);
    
}

void PangenomeMAT2::Tree::printVCFParallel(std::string reference, std::ofstream& fout) {

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "Error: Reference sequence with matching name not found!"){
        std::cerr << referenceSequence << std::endl;
        return;
    }

    size_t recordID = 0;

    std::mutex vcfMapMutex;
    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    tbb::parallel_for_each(allNodes, [&](auto& n){
        if(n.first != reference){
            std::string altSequence = getStringFromReference(n.first);
            if(altSequence.length() != referenceSequence.length()){
                std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
                return;
            }

            std::string currentRefString, currentAltString;
            int currentCoordinate = 0;

            int diffStart = 0;

            for(size_t i = 0; i < referenceSequence.length(); i++){

                if(referenceSequence[i] == '-' && altSequence[i] == '-'){
                    continue;
                } else if(referenceSequence[i] != '-' && altSequence[i] == '-'){
                    if(currentRefString == "" && currentAltString == ""){
                        diffStart = currentCoordinate;
                    }

                    currentRefString += referenceSequence[i];
                } else if(referenceSequence[i] == '-' && altSequence[i] != '-'){
                    if(currentRefString == "" && currentAltString == ""){
                        diffStart = currentCoordinate;
                    }

                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] != altSequence[i]){
                    if(currentRefString == "" && currentAltString == ""){
                        diffStart = currentCoordinate;
                    }
                    if(currentRefString == currentAltString){
                        currentRefString = "";
                        currentAltString = "";
                        diffStart = currentCoordinate;
                    }
                    currentRefString += referenceSequence[i];
                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] == altSequence[i]){
                    if(currentRefString == currentAltString){
                        // Reset
                        diffStart = currentCoordinate;
                        currentRefString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString = currentRefString;
                    } else {
                        // Create VCF record at position i
                        if(currentRefString == ""){
                            currentRefString += referenceSequence[i];
                            currentAltString += altSequence[i];
                            diffStart = currentCoordinate;
                            vcfMapMutex.lock();
                            vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                            vcfMapMutex.unlock();
                            diffStart = currentCoordinate+1;
                            currentRefString = "";
                            currentAltString = "";
                        } else {
                            vcfMapMutex.lock();
                            vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                            vcfMapMutex.unlock();

                            // Reset
                            diffStart = currentCoordinate;
                            currentRefString = "";
                            currentRefString += referenceSequence[i];
                            currentAltString = currentRefString;
                        }
                    }
                }

                if(referenceSequence[i] != '-'){
                    currentCoordinate++;
                }
            }

            if(currentRefString != currentAltString){
                vcfMapMutex.lock();
                vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                vcfMapMutex.unlock();

                // Reset
                diffStart = referenceSequence.size();
                currentRefString = "";
                currentAltString = currentRefString;
            }
        }
    });

    std::mutex sequenceIdsMutex;
    std::map< std::string, size_t > sequenceIds;
    tbb::parallel_for_each(allNodes, [&](auto& u){
        //if(u.second->children.size() == 0 && u.first != reference){
        if(u.first != reference){

            sequenceIdsMutex.lock();
            sequenceIds[u.first] = 0;
            sequenceIdsMutex.unlock();
        }
    });


    fout << "##fileformat=VCFv" << VCF_VERSION << '\n';
    fout << "##fileDate=" << PangenomeMAT2::getDate() << '\n';
    fout << "##source=PanMATv" << PMAT_VERSION << '\n';
    fout << "##reference=" << reference << '\n';
    fout << "#CHROM\t" << "POS\t" << "ID\t" << "REF\t" << "ALT\t" << "QUAL\t" << "FILTER\t" << "INFO\t" << "FORMAT\t";
    
    // fout << std::left << std::setw(20) << "#CHROM " << std::setw(20) << "POS " << std::setw(20) << "ID " << std::setw(20) << "REF " << std::setw(20) << "ALT " << std::setw(20) << "QUAL " << std::setw(20) << "FILTER " << std::setw(20) << "INFO " << std::setw(20) << "FORMAT ";
    for(auto u: sequenceIds){
        fout << u.first + "\t";
    }
    fout << '\n';

    for(auto u: vcfMap){
        for(auto v: u.second){
            if(v.first == ""){
                fout << ".\t" << u.first << "\t" << recordID++ << "\t" << ".\t";
            } else {
                fout << ".\t" << u.first << "\t" << recordID++ << "\t" << v.first << "\t";
            }
            
            std::map< std::string, size_t > tempSequenceIds = sequenceIds;

            int ctr = 1;
            std::string altStrings;

            for(auto w: v.second){
                altStrings += (w.first == "" ? ".": w.first);
                altStrings += ",";
                for(auto uu: w.second){
                    tempSequenceIds[uu] = ctr;
                }
                ctr++;
            }

            altStrings.pop_back();

            fout << altStrings << "\t.\t.\t.\t.\t";

            for(auto w: tempSequenceIds){
                fout << w.second << "\t";
            }

            fout << '\n';
        }
    }
}
std::vector< std::string > PangenomeMAT2::Tree::searchByAnnotation(std::string annotation){
    if(annotationsToNodes.find(annotation) != annotationsToNodes.end()){
        return annotationsToNodes[annotation];
    }
    return {};
}

void PangenomeMAT2::Tree::annotate(std::ifstream& fin){
    std::string line;
    while(getline(fin, line)){
        std::string word;
        std::string nodeId;

        // Extract node ID
        size_t i = 0;
        for(;i < line.length() && line[i]!=','; i++){
            word+=line[i];
        }

        word = stripString(word);

        if(word.length()){
            nodeId = word;
            word = "";
        } else {
            std::cout << "File in incorrect format. Line: " << line << std::endl;
            return;
        }

        if(i >= line.length()){
            // comma not found
            std::cout << "File in incorrect format. Line: " << line << std::endl;
            return;
        }

        if(allNodes.find(nodeId) == allNodes.end()){
            std::cout << "Node ID not found. Line: " << line << std::endl;
            return;
        }

        Node* nodeToAnnotate = allNodes[nodeId];

        // Extract annotations
        for(;i < line.length(); i++){
            if(line[i] != ','){
                word += line[i];
            } else {
                word = stripString(word);
                if(word.length()){
                    std::string annotation = word;
                    nodeToAnnotate->annotations.push_back(annotation);
                    annotationsToNodes[annotation].push_back(nodeId);
                    word = "";
                }
            }
        }

        word = stripString(word);
        if(word.length()){
            std::string annotation = word;
            nodeToAnnotate->annotations.push_back(annotation);
            annotationsToNodes[annotation].push_back(nodeId);
            word = "";
        }

    }
}
