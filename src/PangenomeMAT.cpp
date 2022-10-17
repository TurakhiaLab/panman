#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>

#include "PangenomeMAT.hpp"

// #define DEBUG

PangenomeMAT::Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

PangenomeMAT::Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

PangenomeMAT::Block::Block(MAT::block b){
    blockId = b.block_id();
    chromosomeName = b.chromosome_name();
    for(int i = 0; i < b.consensus_seq_size(); i++){
        consensusSeq.push_back(b.consensus_seq(i));
    }
}

void PangenomeMAT::stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    // TIMEIT();
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        // if ((end_pos == start_pos) || end_pos >= s.length()) {
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

PangenomeMAT::Node* PangenomeMAT::Tree::createTreeFromNewickString(std::string newickString) {

    PangenomeMAT::Node* newTreeRoot;

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

void PangenomeMAT::Tree::assignMutationsToNodes(Node* root, size_t currentIndex, std::vector< MAT::node >& nodes){
    std::vector< PangenomeMAT::NucMut > storedNucMutation;
    for(int i = 0; i < nodes[currentIndex].nuc_mutation_size(); i++){
        storedNucMutation.push_back( PangenomeMAT::NucMut(nodes[currentIndex].nuc_mutation(i)) );
    }

    PangenomeMAT::BlockMut storedBlockMutation;
    storedBlockMutation.loadFromProtobuf(nodes[currentIndex].block_mutation());

    root->nucMutation = storedNucMutation;
    root->blockMutation = storedBlockMutation;

    for(auto child: root->children){
        assignMutationsToNodes(child, currentIndex+1, nodes);
    }

}

PangenomeMAT::Tree::Tree(std::ifstream& fin){

    MAT::tree mainTree;

    if(!mainTree.ParseFromIstream(&fin)){
        throw std::invalid_argument("Could not read tree from input file.");
    }

    // Create tree
    root = createTreeFromNewickString(mainTree.newick());

    // printBfs();

    std::vector< MAT::node > storedNodes;
    for(int i = 0; i < mainTree.nodes_size(); i++){
        storedNodes.push_back(mainTree.nodes(i));
    }

    assignMutationsToNodes(root, 0, storedNodes);

    // Block sequence
    for(int i = 0; i < mainTree.blocks_size(); i++){
        blocks.emplace_back(mainTree.blocks(i));
        // std::cout << blocks[blocks.size() - 1].blockId << " ";
    }

    // Gap List
    for(int i = 0; i < mainTree.gaps().position_size(); i++){
        gaps.position.push_back(mainTree.gaps().position(i));
    }
    for(int i = 0; i < mainTree.gaps().condensed_size(); i++){
        gaps.condensed.push_back(mainTree.gaps().condensed(i));
    }

}

int getTotalParsimonyParallelHelper(PangenomeMAT::Node* root, PangenomeMAT::NucMutationType nucMutType, PangenomeMAT::BlockMutationType blockMutType){
    int totalMutations = 0;

    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->nucMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++){
            if((root->nucMutation[i].condensed & 0x3) == nucMutType){
                if(nucMutType == PangenomeMAT::NucMutationType::NS){
                    init += ((root->nucMutation[i].condensed) & (((1<<6)-1)<<2)); // Length of contiguous mutation in case of substitution
                } else {
                    init++;
                }
                // init++;
            }
        }
        return init;
    }, [&](int x, int y){
        return x + y;
    });

    if(blockMutType != PangenomeMAT::BlockMutationType::NONE){
        totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.condensedBlockMut.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
            for(int i = r.begin(); i != r.end(); i++){
                if((root->blockMutation.condensedBlockMut[i] & 0x1) == blockMutType){
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

int PangenomeMAT::Tree::getTotalParsimonyParallel(NucMutationType nucMutType, BlockMutationType blockMutType){

    return getTotalParsimonyParallelHelper(root, nucMutType, blockMutType);

}

int PangenomeMAT::Tree::getTotalParsimony(PangenomeMAT::NucMutationType nucMutType, PangenomeMAT::BlockMutationType blockMutType){
    int totalMutations = 0;

    std::queue<Node *> bfsQueue;

    bfsQueue.push(root);
    while(!bfsQueue.empty()){
        Node* current = bfsQueue.front();
        bfsQueue.pop();

        // Process children of current node
        for(auto nucMutation: current->nucMutation){
            if((nucMutation.condensed & 0x3) == nucMutType){
                if(nucMutType == PangenomeMAT::NucMutationType::NS){
                    totalMutations += (nucMutation.condensed & (((1<<6)-1)<<2)); // Length of contiguous mutation in case of substitution
                } else {
                    totalMutations++;
                }
            }
        }

        if(blockMutType != NONE){
            for(auto blockMutation: current->blockMutation.condensedBlockMut){
                if((blockMutation & 0x1) == blockMutType){
                    totalMutations++;
                }
            }
        }

        for(auto child: current->children){
            bfsQueue.push(child);
        }
    }

    return totalMutations;
}

void PangenomeMAT::Tree::printSummary(){
    // Traversal test

#ifdef DEBUG
    std::cout << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    std::cout << "Total Samples in Tree: " << m_numLeaves << std::endl;
    std::cout << "Total Substitutions: " << getTotalParsimony(PangenomeMAT::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimony(PangenomeMAT::NucMutationType::NI, PangenomeMAT::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimony(PangenomeMAT::NucMutationType::ND, PangenomeMAT::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP mutations: " << getTotalParsimony(PangenomeMAT::NucMutationType::NSNP) << std::endl;
    std::cout << "Max Tree Depth: " << m_maxDepth << std::endl;
    std::cout << "Mean Tree Depth: " << m_meanDepth << std::endl;

    std::cout << "\nParallel Results:\n";
    std::cout << "Total Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NI, PangenomeMAT::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::ND, PangenomeMAT::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP mutations: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNP) << std::endl;

#else
    std::cout << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    std::cout << "Total Samples in Tree: " << m_numLeaves << std::endl;
    std::cout << "Total Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NI, PangenomeMAT::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::ND, PangenomeMAT::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP mutations: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNP) << std::endl;
    std::cout << "Max Tree Depth: " << m_maxDepth << std::endl;
    std::cout << "Mean Tree Depth: " << m_meanDepth << std::endl;
#endif

}

void PangenomeMAT::Tree::printBfs(){
    // Traversal test
    std::queue<Node *> bfsQueue;
    size_t prevLev = 0;
    bfsQueue.push(root);
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
}

std::vector< std::string > getSequenceLines(std::vector< std::vector< std::pair< char, std::vector< char > > > >& sequence,\
    std::vector<bool>& blockExists, size_t lineSize, bool aligned){

    std::vector< std::string > lines;
    std::string line;

    for(size_t i = 1; i < blockExists.size(); i++){
        if(!blockExists[i]){
            if(aligned){
                for(size_t j = 0; j < sequence[i].size(); j++){
                    for(size_t k = 0; k < sequence[i][j].second.size(); k++){
                        line += '-';
                        if(line.length() == lineSize){
                            lines.push_back(line);
                            line = "";
                        }
                    }
                    if(sequence[i][j].first != 'x'){
                        line+='-';
                        if(line.length() == lineSize){
                            lines.push_back(line);
                            line = "";
                        }
                    }
                }
            }

            continue;
        }

        for(size_t j = 0; j < sequence[i].size(); j++){
            for(size_t k = 0; k < sequence[i][j].second.size(); k++){
                if(sequence[i][j].second[k] == '-'){
                    if(aligned){
                        line += '-';
                    }
                    
                } else {
                    line += sequence[i][j].second[k];
                }
                if(line.length() == lineSize){
                    lines.push_back(line);
                    line = "";
                }
            }

            if(sequence[i][j].first != 'x'){
                if(sequence[i][j].first == '-'){
                    if(aligned){
                        line += '-';
                    }
                } else {
                    line += sequence[i][j].first;
                }
                if(line.length() == lineSize){
                    lines.push_back(line);
                    line = "";
                }
            }
        }

    }

    if(line.length()){
        lines.push_back(line);
        line = "";
    }

    return lines;

}

void printFASTAHelper(PangenomeMAT::Node* root,\
    std::vector< std::vector< std::pair< char, std::vector< char > > > >& sequence, std::vector<bool>& blockExists,\
    std::ofstream& fout){
    
    // Apply mutations
    // Block mutations
    for(uint32_t mutation: root->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & ((1 << 24) - 1));
        int type = (mutation & 0x1);

        // std::cout << mutation << " " << bid  << " " << sequence.size() << " " << type << std::endl;

        if(type == 0){
            blockExists[bid] = true;
        } else {
            blockExists[bid] = false;
        }

    }

    // For backtracking. blockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++){

        int bid = ((root->nucMutation[i].condensed >> 8) & (((1 << 24) - 1)));

        // std::cout << bid << " " << blockExists[bid] << std::endl;

        // std::cout << '>' << root->identifier << " " << i << " " << root->nucMutation.size() << " " << bid << " " <<sequence.size() << std::endl;

        if(!blockExists[bid]){
            continue;
        }

        int pos = root->nucMutation[i].position;
        int gapPos = root->nucMutation[i].gapPosition;
        int type = ((root->nucMutation[i].condensed) & 3);

        if(type < 3){
            // Either S, I or D
            int len = (((root->nucMutation[i].condensed) >> 2) & 0x3F);
            char newVal;

            if(type == PangenomeMAT::NucMutationType::NS){
                for(int j = 0; j < len; j++){
                    char oldVal = sequence[bid][pos + j].first;
                    switch(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15){
                        case 1:
                            newVal = 'A';
                            break;
                        case 2:
                            newVal = 'C';
                            break;
                        case 4:
                            newVal = 'G';
                            break;
                        case 8:
                            newVal = 'T';
                            break;
                        case 15:
                            newVal = 'N';
                            break;
                    }


                    sequence[bid][pos + j].first = newVal;
                    mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, oldVal, newVal));
                }
            } else if(type == PangenomeMAT::NucMutationType::NI){
                
                if(root->nucMutation[i].gapPosition == -1){
                    for(int j = 0; j < len; j++){
                        switch(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15){
                            case 1:
                                newVal = 'A';
                                break;
                            case 2:
                                newVal = 'C';
                                break;
                            case 4:
                                newVal = 'G';
                                break;
                            case 8:
                                newVal = 'T';
                                break;
                            case 15:
                                newVal = 'N';
                                break;
                        }

                        sequence[bid][pos + j].first = newVal;
                        mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, '-', newVal));
                    }
                } else {
                    for(int j = 0; j < len; j++){
                        switch(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15){
                            case 1:
                                newVal = 'A';
                                break;
                            case 2:
                                newVal = 'C';
                                break;
                            case 4:
                                newVal = 'G';
                                break;
                            case 8:
                                newVal = 'T';
                                break;
                            case 15:
                                newVal = 'N';
                                break;
                        }
                        sequence[bid][pos].second[gapPos + j] = newVal;
                        mutationInfo.push_back(std::make_tuple(bid, pos, gapPos + j, '-', newVal));
                    }
                }
            } else if(type == PangenomeMAT::NucMutationType::ND){
                if(root->nucMutation[i].gapPosition == -1){
                    for(int j = 0; j < len; j++){
                        char oldVal = sequence[bid][pos + j].first;
                        sequence[bid][pos + j].first = '-';
                        mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, oldVal, '-'));
                    }
                } else {
                    for(int j = 0; j < len; j++){
                        char oldVal = sequence[bid][pos].second[gapPos + j];
                        sequence[bid][pos].second[gapPos + j] = '-';
                        mutationInfo.push_back(std::make_tuple(bid, pos, gapPos + j, oldVal, '-'));
                    }
                }
            }
        } else {
            // Todo: SNP

        }
    }

    if(root->children.size() == 0){
        // Print sequence
        fout << '>' << root->identifier << std::endl;
        auto lines = getSequenceLines(sequence, blockExists, 50, false);
        for(auto line: lines){
            fout << line << '\n';
        }
        fout << '\n';

    } else {
        // DFS on children
        for(PangenomeMAT::Node* child: root->children){
            printFASTAHelper(child, sequence, blockExists, fout);
        }
    }

    // Undo mutations
    for(auto mutation: mutationInfo){
        if(std::get<2>(mutation) == -1){
            sequence[std::get<0>(mutation)][std::get<1>(mutation)].first = std::get<3>(mutation);
        } else {
            sequence[std::get<0>(mutation)][std::get<1>(mutation)].second[std::get<2>(mutation)] = std::get<3>(mutation);
        }
    }

}

void PangenomeMAT::Tree::printFASTA(std::ofstream& fout){
    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    // First block is empty since the block IDs start at 1
    std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence(blocks.size() + 1);
    
    std::vector< bool > blockExists(blocks.size() + 1, false);

    // std::cout << blocks.size() << " " << blocks[0].consensusSeq.size() << std::endl;

    for(size_t i = 0; i < blocks.size(); i++){
        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            // std::cout << i << " " << j << std::endl;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                switch(nucCode){
                    case 1:
                        sequence[blocks[i].blockId].push_back({'A', {}});
                        break;
                    case 2:
                        sequence[blocks[i].blockId].push_back({'C', {}});
                        break;
                    case 4:
                        sequence[blocks[i].blockId].push_back({'G', {}});
                        break;
                    case 8:
                        sequence[blocks[i].blockId].push_back({'T', {}});
                        break;
                    case 15:
                        sequence[blocks[i].blockId].push_back({'N', {}});
                        break;
                }
            }
        }
        // std::cout << blocks[i].consensusSeq.size() << " " << sequence[i].size() << std::endl;

        // End character to incorporate for gaps at the end
        sequence[blocks[i].blockId].push_back({'x',{}});
    }

    // Assigning gaps
    for(size_t i = 0; i < gaps.position.size(); i++){
        int bId = ((gaps.condensed[i] >> 8) & (((1 << 24) - 1)));

        int len = ((gaps.condensed[i]) & 255);
        int pos = gaps.position[i];

        // if(bId == 53)
        //     std::cout << bId << " " << len << " " << pos << " " << sequence[bId].size() << std::endl;

        sequence[bId][pos].second.resize(len, '-');
    }

    printFASTAHelper(root, sequence, blockExists, fout);

}