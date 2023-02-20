#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <ctime>
#include <iomanip>
#include <mutex>

#include "PangenomeMATV2.hpp"

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

    for(int i = 0; i < nodes[currentIndex].mutations_size(); i++){
        if(nodes[currentIndex].mutations(i).mutationtype() == false){
            for(auto nucMut: nodes[currentIndex].mutations(i).nucmutation()){
                storedNucMutation.push_back( PangenomeMAT2::NucMut(nucMut, nodes[currentIndex].mutations(i).blockid(), nodes[currentIndex].mutations(i).blockgapexist()) );
            }
        }
    }

    std::vector< PangenomeMAT2::BlockMut > storedBlockMutation;
    for(int i = 0; i < nodes[currentIndex].mutations_size(); i++){
        PangenomeMAT2::BlockMut tempBlockMut;
        if(nodes[currentIndex].mutations(i).mutationtype()){
            tempBlockMut.loadFromProtobuf(nodes[currentIndex].mutations(i));
            storedBlockMutation.push_back(tempBlockMut);
        }
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
                            // if(primaryBlockId >= sequence.size()){
                            //     std::cout << "A " << root->identifier << " " << primaryBlockId << " " << sequence.size() << std::endl;
                            // }
                            // if(nucPosition >= sequence[primaryBlockId].first.size()){
                            //     std::cout << "B " << root->identifier << " " << nucPosition << " " << sequence[primaryBlockId].first.size() << std::endl;
                            // }
                            // if(nucGapPosition+j >= sequence[primaryBlockId].first[nucPosition].second.size()){
                            //     std::cout << "C " << primaryBlockId << " " << nucPosition << " " << i << " " << root->identifier << " " << nucGapPosition << " " << nucGapPosition+j << " " << sequence[primaryBlockId].first[nucPosition].second.size() << std::endl;
                            // }
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            // if(nucGapPosition+j >= sequence[primaryBlockId].first[nucPosition].second.size()){
                            //     std::cout << "C " << i << " " << root->identifier << " " << nucGapPosition << " " << nucGapPosition+j << " " << sequence[primaryBlockId].first[nucPosition].second.size() << std::endl;
                            // }
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
                            // if(nucGapPosition+j >= sequence[primaryBlockId].first[nucPosition].second.size()){
                            //     std::cout << "C " << primaryBlockId << " " << nucPosition << " " << i << " " << root->identifier << " " << nucGapPosition << " " << nucGapPosition+j << " " << sequence[primaryBlockId].first[nucPosition].second.size() << std::endl;
                            // }
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            // if(nucGapPosition+j >= sequence[primaryBlockId].first[nucPosition].second.size()){
                            //     std::cout << "C " << primaryBlockId << " " << nucPosition << " " << i << " " << root->identifier << " " << nucGapPosition << " " << nucGapPosition+j << " " << sequence[primaryBlockId].first[nucPosition].second.size() << std::endl;
                            // }
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
    
    if(root->children.size() == 0){
        // Print sequence
        fout << '>' << root->identifier << std::endl;

        PangenomeMAT2::printSequenceLines(sequence, blockExists, 70, aligned, fout);

    } else {
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
    std::map< std::pair< int32_t, int32_t >, std::vector< MATNew::nucMut > > blockToNucMutations;

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

        nm.set_mutinfo((((mutation.nucs) >> (24 - (mutation.mutInfo >> 4)*4)) << 8) + mutation.mutInfo);
        blockToNucMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].push_back(nm);
    }

    for(auto u: blockToNucMutations){
        MATNew::mutation mutation;
        mutation.set_mutationtype(false);
        int32_t primaryBlockId = u.first.first;
        int32_t secondaryBlockId = u.first.second;
        if(secondaryBlockId != -1){
            mutation.set_blockid(((int64_t)primaryBlockId << 32) + secondaryBlockId);
            mutation.set_blockgapexist(true);
        } else {
            mutation.set_blockid(((int64_t)primaryBlockId << 32));
            mutation.set_blockgapexist(false);
        }
        for(auto v: u.second){
            mutation.add_nucmutation();
            *mutation.mutable_nucmutation(mutation.nucmutation_size() - 1) = v;
        }
        n.add_mutations();
        *n.mutable_mutations(n.mutations_size() - 1) = mutation;
    }

    for(size_t i = 0; i < root->blockMutation.size(); i++){
        const PangenomeMAT2::BlockMut& mutation = root->blockMutation[i];

        MATNew::mutation bm;
        bm.set_mutationtype(true);

        if(mutation.secondaryBlockId != -1){
            bm.set_blockid(((int64_t)mutation.primaryBlockId << 32) + mutation.secondaryBlockId);
            bm.set_blockgapexist(true);
        } else {
            bm.set_blockid(((int64_t)mutation.primaryBlockId << 32));
            bm.set_blockgapexist(false);
        }
        bm.set_blockmutinfo(mutation.blockMutInfo);

        n.add_mutations();
        *n.mutable_mutations(i) = bm;
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
        if(u.second->children.size() == 0 && u.first == reference){
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

    // if(reference == "ON650465.1"){
    //         // getStringFromReference("ON650465.1", false);
    //         for(int i = sequence[0].first.size()-10; i < sequence[0].first.size(); i++){
    //             for(int j = 0; j < sequence[0].first[i].second.size(); j++){
    //                 if(sequence[0].first[i].second[j]!='-')std::cout << sequence[0].first[i].second[j];
    //             }
    //             std::cout << sequence[0].first[i].first;
    //         }
    //         std::cout << std::endl;
    //     }

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

    // if(reference == "ON650465.1"){
    //     for(auto p: path){
    //             std::cout << p->identifier << std::endl;
    //         }
    //     for(int i = sequence[0].first.size()-10; i < sequence[0].first.size(); i++){
    //         for(int j = 0; j < sequence[0].first[i].second.size(); j++){
    //             if(sequence[0].first[i].second[j]!='-')std::cout << sequence[0].first[i].second[j];
    //         }
    //         std::cout << sequence[0].first[i].first;
    //     }
    //     std::cout << std::endl;
    // }

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

std::string PangenomeMAT2::stripGaps(const std::string sequenceString){
    std::string result;
    for(auto u: sequenceString){
        if(u != '-' && u != 'x'){
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
        if(n.second->children.size() == 0){
            std::string sequence;
            sequence = getStringFromReference(n.first, aligned);
            
            fastaMutex.lock();
            fout << '>' << n.first << '\n';
            for(size_t i = 0; i < sequence.size(); i+=lineSize){
                fout << sequence.substr(i, std::min(lineSize, sequence.size() - i)) << '\n';
            }
            fastaMutex.unlock();
        }

    });
}

void PangenomeMAT2::Tree::printVCFParallel(std::string reference, std::ofstream& fout){

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "Error: Reference sequence with matching name not found!"){
        std::cerr << referenceSequence << std::endl;
        return;
    }

    size_t recordID = 0;

    std::mutex vcfMapMutex;
    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    tbb::parallel_for_each(allNodes, [&](auto& n){
        if(n.second->children.size() == 0 && n.first != reference){
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
        if(u.second->children.size() == 0 && u.first != reference){
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
        if(u.first != sequenceIds.rbegin()->first){
            fout << u.first + "\t";
        } else {
            fout << u.first;
        }
    }
    fout << '\n';

    for(auto u: vcfMap){
        for(auto v: u.second){
            if(v.first == ""){
                fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << ".\t";
            } else {
                fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << v.first << "\t";
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
                if(w.first != sequenceIds.rbegin()->first){
                    fout << w.second << "\t";
                } else {
                    fout << w.second;
                }
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

void PangenomeMAT2::Tree::convertToGFA(std::ofstream& fout){

    size_t autoIncrId = 0;
    // std::map< std::tuple< int32_t, int32_t, int32_t, int32_t, std::string >, uint64_t > sequenceNodes;
    std::map< std::pair< size_t, std::string >, size_t > allSequenceNodes;
    std::map< std::string, std::vector< size_t > > paths;

    for(auto u: allNodes){
        if(u.second->children.size() != 0){
            continue;
        }

        std::string alignedSequence = getStringFromReference(u.first, true);

        std::string currentSequence;
        size_t currentStart = 0;
        std::vector< size_t > sequenceNodeIds;

        for(size_t i = 0; i < alignedSequence.length(); i++){
            if(currentSequence.length() == 0){
                currentStart = i;
            }
            currentSequence += alignedSequence[i];
            if(currentSequence.length() == 32){
                currentSequence = stripGaps(currentSequence);
                if(currentSequence.length()){
                    if(allSequenceNodes.find(std::make_pair(currentStart, currentSequence)) == allSequenceNodes.end()){
                        allSequenceNodes[std::make_pair(currentStart, currentSequence)] = autoIncrId;
                        sequenceNodeIds.push_back(autoIncrId);
                        autoIncrId++;
                    } else {
                        sequenceNodeIds.push_back(allSequenceNodes[std::make_pair(currentStart, currentSequence)]);
                    }
                }
                currentSequence = "";
            }
        }
        if(currentSequence.length()){
            currentSequence = stripGaps(currentSequence);
            if(currentSequence.length()){
                if(allSequenceNodes.find(std::make_pair(currentStart, currentSequence)) == allSequenceNodes.end()){
                    allSequenceNodes[std::make_pair(currentStart, currentSequence)] = autoIncrId;
                    sequenceNodeIds.push_back(autoIncrId);
                    autoIncrId++;
                } else {
                    sequenceNodeIds.push_back(allSequenceNodes[std::make_pair(currentStart, currentSequence)]);
                }
            }
            currentSequence = "";
        }
        paths[u.first] = sequenceNodeIds;
    }

    // for(auto u: paths){
    //     for(size_t i = 1; i < u.second.size(); i++){
    //         edges.insert(std::make_pair(u.second[i-1], u.second[i]));
    //     }
    // }

    std::map< size_t, std::string > finalNodes;
    for(auto u: allSequenceNodes){
        finalNodes[u.second] = u.first.second;
    }


    // Graph and its transpose
    std::vector< std::pair< std::string, size_t > > G[autoIncrId];
    std::vector< std::pair< std::string, size_t > > GT[autoIncrId];

    for(auto u: paths){
        for(size_t i = 1; i < u.second.size(); i++){
            G[u.second[i-1]].push_back(std::make_pair(u.first, u.second[i]));
            GT[u.second[i]].push_back(std::make_pair(u.first, u.second[i-1]));
        }
    }

    for(size_t i = 0; i < autoIncrId; i++){
        // Sort so we can compare the sequence IDs of incoming and outgoing edges for equality
        sort(G[i].begin(), G[i].end());
        sort(GT[i].begin(), GT[i].end());
    }

    for(size_t i = 0; i < autoIncrId; i++){
        if(finalNodes.find(i) == finalNodes.end()){
            continue;
        }

        while(true){
            // check if a node's edges only go to one next next node and there is an outgoing edge for every incoming edge
            if(G[i].size() != GT[i].size() || G[i].size() == 0){
                break;
            }
            bool check = true;

            for(size_t j = 0; j < G[i].size(); j++){
                if(G[i][j].first != GT[i][j].first){
                    check = false;
                    break;
                } else if(j > 0 && G[i][j].second != G[i][j-1].second){
                    check = false;
                    break;
                }
            }

            if(!check){
                break;
            }

            size_t dest = G[i][0].second;
            if(G[i].size() != GT[dest].size()){
                break;
            }

            for(size_t j = 0; j < G[i].size(); j++){
                if(G[i][j].first != GT[dest][j].first && GT[dest][j].second != i){
                    check = false;
                    break;
                }
            }
            if(!check){
                break;
            }
            // combine i and dest
            finalNodes[i] += finalNodes[dest];
            finalNodes.erase(dest);
            G[i].clear();
            G[i] = G[dest];
        }
    }

    std::set< std::pair< size_t, size_t > > edges;
    
    for(size_t i = 0; i < autoIncrId; i++){
        if(finalNodes.find(i) == finalNodes.end()){
            continue;
        }
        for(auto edge: G[i]){
            edges.insert(std::make_pair(i, edge.second));
        }
    }

    for(auto p: paths){
        std::vector< size_t > newPath;
        for(auto n: p.second){
            if(finalNodes.find(n) != finalNodes.end()){
                newPath.push_back(n);
            }
        }
        paths[p.first] = newPath;
    }

    fout << "H\tVN:Z:1.1\n";
    for(auto u: finalNodes){
        fout << "S\t" << u.first << "\t" << u.second << "\n";
    }

    for(auto u: edges){
        fout << "L\t" << u.first << "\t+\t" << u.second << "\t+\t0M\n";
    }

    for(auto u: paths){
        fout << "P\t" << u.first << "\t";
        for(size_t i = 0; i < u.second.size(); i++){
            fout << u.second[i] << "+";
            if(i != u.second.size() - 1){
                fout << ",";
            }
        }
        fout << "\t*\n";
    }
}

void PangenomeMAT2::Tree::printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout){
    std::map< std::string, std::string > nodes;
    std::map< std::string, std::vector< std::string > > paths;
    std::string line;
    while(getline(fin, line, '\n')){
        std::vector< std::string > separatedLine;
        stringSplit(line, '\t', separatedLine);
        if(separatedLine[0] == "S"){
            nodes[separatedLine[1]] = separatedLine[2];
        } else if(separatedLine[0] == "P"){
            std::vector< std::string > v;
            stringSplit(separatedLine[2], ',', v);
            for(size_t i = 0; i < v.size(); i++){
                v[i].pop_back();
            }
            paths[separatedLine[1]] = v;
        }
    }
    for(auto p: paths){
        fout << ">" << p.first << "\n";
        std::string sequence;
        for(auto s: p.second){
            sequence += nodes[s];
        }
        for(size_t i = 0; i < sequence.size(); i+=70){
            fout << sequence.substr(i,std::min((size_t)70, sequence.size() - i)) << '\n';
        }
    }
}

void PangenomeMAT2::Tree::convertToVG(std::ofstream& fout){
    vg::Graph g;
    int32_t maxBlockId = 0;
    
    std::vector< std::pair< bool, std::vector< bool > > > globalBlockExists(blocks.size() + 1, {false, {}});
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > globalSequence(blocks.size() + 1);

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        globalSequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        globalBlockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    }

    for(size_t i = 0; i < blocks.size(); i++){
        int64_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int64_t secondaryBlockId = ((int64_t)blocks[i].secondaryBlockId & 0xFFFFFFFF);

        maxBlockId = std::max((int64_t)maxBlockId, primaryBlockId);
        uint64_t fullBlockId = ((primaryBlockId << 32)+secondaryBlockId);

        secondaryBlockId = (secondaryBlockId == 0xFFFFFFFF? -1: secondaryBlockId);

        std::string blockString;

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                
                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);

                blockString += nucleotide;
                if(secondaryBlockId != -1){
                    globalSequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    globalSequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }
            if(endFlag){
                break;
            }
        }
        if(secondaryBlockId != -1){
            globalSequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            globalSequence[primaryBlockId].first.push_back({'x', {}});
        }
        vg::Node* n = g.add_node();
        n->set_sequence(blockString);
        n->set_id(fullBlockId);
        n->set_name("");
    }
    globalSequence.resize(maxBlockId + 1);
    globalBlockExists.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId & 0xFFFFFFFF);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            if(secondaryBId != -1){
                globalSequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                globalSequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }

    for(auto u: allNodes){
        if(u.second->children.size() != 0){
            continue; 
        }
        Node* it = u.second;
        std::vector< PangenomeMAT2::Node* > path;
        while(it != root){
            path.push_back(it);
            it = it->parent;
        }
        path.push_back(root);

        std::vector< std::pair< bool, std::vector< bool > > > blockExists = globalBlockExists;
        std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence = globalSequence;

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

        std::vector< uint64_t > blockIdSequence;
        for(size_t i = 0; i < blockExists.size(); i++){
            for(size_t j = 0; j < blockExists[i].second.size(); i++){
                if(blockExists[i].second[j]){
                    int64_t primaryBlockId = i;
                    int64_t secondaryBlockId = j;
                    uint64_t fullBlockId = (primaryBlockId << 32)+secondaryBlockId;
                    blockIdSequence.push_back(fullBlockId);
                }
            }
            if(blockExists[i].first){
                int64_t primaryBlockId = i;
                int64_t secondaryBlockId = ((int64_t)(-1) & 0xFFFFFFFF);
                uint64_t fullBlockId = (primaryBlockId << 32)+secondaryBlockId;
                blockIdSequence.push_back(fullBlockId);
            }
        }

        for(size_t i = 0; i < blockIdSequence.size() - 1; i++){
            vg::Edge* e = g.add_edge();
            e->set_from(blockIdSequence[i]);
            e->set_to(blockIdSequence[i+1]);
            e->set_from_start(false);
            e->set_to_end(false);
            e->set_overlap(0);
        }

        vg::Path* p = g.add_path();
        p->set_name(u.first);
        p->set_is_circular(false);

        for(size_t i = 0; i < blockIdSequence.size(); i++){

            // compute edits
            // First compute aligned sequence strings
            std::string referenceSequence;
            std::string altSequence;
            int32_t primaryBlockId = (blockIdSequence[i] >> 32);
            int32_t secondaryBlockId = ((blockIdSequence[i] & 0xFFFFFFFF) == 0xFFFFFFFF? -1: (blockIdSequence[i] & 0xFFFFFFFF));
            std::vector< std::pair< char, std::vector< char > > > globalBlockVector;
            std::vector< std::pair< char, std::vector< char > > > localBlockVector;

            if(secondaryBlockId != -1){
                globalBlockVector = globalSequence[primaryBlockId].second[secondaryBlockId];
                localBlockVector = sequence[primaryBlockId].second[secondaryBlockId];
            } else {
                globalBlockVector = globalSequence[primaryBlockId].first;
                localBlockVector = sequence[primaryBlockId].first;
            }

            for(size_t i = 0; i < globalBlockVector.size(); i++){
                for(size_t j = 0; j < globalBlockVector[i].second.size(); j++){
                    if(globalBlockVector[i].second[j] != 'x'){
                        referenceSequence += globalBlockVector[i].second[j];
                    } else {
                        referenceSequence += '-';
                    }
                    if(localBlockVector[i].second[j] != 'x'){
                        altSequence += localBlockVector[i].second[j];
                    } else {
                        altSequence += '-';
                    }
                }
                if(globalBlockVector[i].first != 'x'){
                    referenceSequence += globalBlockVector[i].first;
                } else {
                    referenceSequence += '-';
                }
                if(localBlockVector[i].first != 'x'){
                    altSequence += localBlockVector[i].first;
                } else {
                    altSequence += '-';
                }
            }

            std::vector< std::tuple< int, int, int, std::string > > vcfMap;

            std::string currentRefString, currentAltString;
            int currentCoordinate = 0;

            int diffStart = 0;

            for(size_t i = 0; i < referenceSequence.length(); i++){

                if(referenceSequence[i] == '-' && altSequence[i] == '-'){
                    continue;
                } else if(referenceSequence[i] != '-' && altSequence[i] == '-'){
                    if(currentAltString == "" && currentRefString != ""){
                        // Deletion extends
                        currentRefString += referenceSequence[i];
                    } else if(currentAltString != "" || currentRefString != ""){
                        // Save old. Start new
                        if(currentRefString != currentAltString){
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), currentAltString));
                        } else {
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), ""));
                        }
                        currentAltString = "";
                        currentRefString = "";
                        currentRefString += referenceSequence[i];
                        diffStart = currentCoordinate;
                    } else {
                        currentRefString += referenceSequence[i];
                        diffStart = currentCoordinate;
                    }
                } else if(referenceSequence[i] == '-' && altSequence[i] != '-'){
                    if(currentRefString == "" && currentAltString != ""){
                        // Insertion extends
                        currentAltString += altSequence[i];
                    } else if(currentRefString != "" || currentAltString != ""){
                        // Save old. Start new
                        if(currentRefString != currentAltString){
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), currentAltString));
                        } else {
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), ""));
                        }
                        currentRefString = "";
                        currentAltString = "";
                        currentAltString += altSequence[i];
                        diffStart = currentCoordinate;
                    } else {
                        currentAltString += altSequence[i];
                        diffStart = currentCoordinate;
                    }
                } else if(referenceSequence[i] != altSequence[i]){
                    if(currentRefString != currentAltString && currentRefString.length() == currentAltString.length()){
                        // Substitution extends
                        currentRefString += referenceSequence[i];
                        currentAltString += altSequence[i];
                    } else if(currentRefString != "" || currentAltString != "") {
                        if(currentRefString != currentAltString){
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), currentAltString));
                        } else {
                            vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), ""));
                        }
                        currentRefString = "";
                        currentAltString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString += altSequence[i];
                        diffStart = currentCoordinate;
                    } else {
                        currentRefString = "";
                        currentAltString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString += altSequence[i];
                        diffStart = currentCoordinate;
                    }
                } else if(referenceSequence[i] == altSequence[i]){
                    
                    if(currentRefString == currentAltString){
                        if(currentRefString == ""){
                            diffStart = currentCoordinate;
                        }
                        // Extend equality
                        currentRefString += referenceSequence[i];
                        currentAltString += altSequence[i];
                    } else {
                        vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), currentAltString));
                        currentRefString = "";
                        currentAltString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString += altSequence[i];
                        diffStart = currentCoordinate;
                    }
                }

                if(referenceSequence[i] != '-'){
                    currentCoordinate++;
                }
            }

            if(currentRefString != currentAltString){
                vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), currentAltString));
            } else if(currentRefString.length() != 0){
                vcfMap.push_back(std::make_tuple(diffStart, currentRefString.length(), currentAltString.length(), ""));
            }

            // Now, add to the path
            vg::Mapping* mp = p->add_mapping();
            mp->set_rank(i);
            vg::Position* pos = mp->mutable_position();
            pos->set_node_id(blockIdSequence[i]);
            pos->set_offset(0);
            pos->set_is_reverse(false);

            for(auto v: vcfMap){
                vg::Edit* ed = mp->add_edit();
                ed->set_from_length(std::get<1>(v));
                ed->set_to_length(std::get<2>(v));
                ed->set_sequence(std::get<3>(v));
            }
        }
    }
    if(!g.SerializeToOstream(&fout)){
        std::cerr << "Error in writing to file" << std::endl;
        return;
    }
}

void PangenomeMAT2::Tree::printFASTAFromVG(std::ifstream& fin, std::ofstream& fout){
    vg::Graph g;
    if(!g.ParseFromIstream(&fin)){
        throw std::invalid_argument("Could not read VG from input file.");
        return;
    }
    for(int i = 0; i < g.path_size(); i++){
        fout << '>' << g.path(i).name() << '\n';
        bool found = false;
        for(int j = 0; j < g.path(i).mapping_size(); j++){
            for(int k = 0; k < g.node_size(); k++){
                if(g.node(k).id() == g.path(i).mapping(j).position().node_id()){
                    // node found
                    found = true;
                    std::string consensusSeq = g.node(k).sequence();
                    // if(g.path(i).name() == "ON650465.1"){
                    //     std::cout << consensusSeq.substr(0,10) << std::endl;
                    // }
                    std::vector< std::pair< char, std::vector< char > > > v;
                    for(auto c: consensusSeq){
                        v.push_back({c,{}});
                    }
                    v.push_back({'-',{}});
                    int currentptr = 0;
                    for(int w = 0; w < g.path(i).mapping(j).edit_size(); w++){
                        int fl = g.path(i).mapping(j).edit(w).from_length();
                        int tl = g.path(i).mapping(j).edit(w).to_length();
                        std::string seq = g.path(i).mapping(j).edit(w).sequence();
                        if(fl == tl && seq.length() == 0){
                            currentptr+=fl;
                            continue;
                        }
                        for(int ii = 0; ii < fl; ii++){
                            v[currentptr+ii].first = '-';
                        }
                        for(int ii = 0; ii < tl; ii++){
                            v[currentptr].second.push_back(seq[ii]);
                        }
                        currentptr+=fl;
                    }
                    std::string newSeq;
                    for(size_t w = 0; w < v.size(); w++){
                        for(size_t ii = 0; ii < v[w].second.size(); ii++){
                            if(v[w].second[ii] != '-'){
                                newSeq += (v[w].second[ii]);
                            }
                        }
                        if(v[w].first != '-'){
                            newSeq += (v[w].first);
                        }
                    }
                    
                    for(size_t w = 0; w < newSeq.size(); w+=70){
                        fout << newSeq.substr(w, std::min(newSeq.size() - w, (size_t)70)) << '\n';
                    }
                    break;
                }
            }
            if(found){
                break;
            }
        }
    }
}