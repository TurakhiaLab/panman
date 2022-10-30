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

void PangenomeMAT::Tree::assignMutationsToNodes(Node* root, size_t& currentIndex, std::vector< MAT::node >& nodes){
    std::vector< PangenomeMAT::NucMut > storedNucMutation;
    for(int i = 0; i < nodes[currentIndex].nuc_mutation_size(); i++){
        storedNucMutation.push_back( PangenomeMAT::NucMut(nodes[currentIndex].nuc_mutation(i)) );
    }

    PangenomeMAT::BlockMut storedBlockMutation;
    storedBlockMutation.loadFromProtobuf(nodes[currentIndex].block_mutation());

    root->nucMutation = storedNucMutation;
    root->blockMutation = storedBlockMutation;

    for(auto child: root->children){
        currentIndex++;
        assignMutationsToNodes(child, currentIndex, nodes);
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

    size_t initialIndex = 0;
    assignMutationsToNodes(root, initialIndex, storedNodes);

    // Block sequence
    for(int i = 0; i < mainTree.blocks_size(); i++){
        blocks.emplace_back(mainTree.blocks(i));
        // std::cout << blocks[blocks.size() - 1].blockId << " ";
    }

    // Gap List
    for(int i = 0; i < mainTree.gaps().position_size(); i++){
        gaps.position.push_back(mainTree.gaps().position(i));
    }
    // for(int i = 0; i < mainTree.gaps().condensed_size(); i++){
    //     gaps.condensed.push_back(mainTree.gaps().condensed(i));
    // }
    for(int i = 0; i < mainTree.gaps().block_id_size(); i++){
        gaps.blockId.push_back(mainTree.gaps().block_id(i));
    }
    for(int i = 0; i < mainTree.gaps().gap_length_size(); i++){
        gaps.gapLength.push_back(mainTree.gaps().gap_length(i));
    }

}

int getTotalParsimonyParallelHelper(PangenomeMAT::Node* root, PangenomeMAT::NucMutationType nucMutType, PangenomeMAT::BlockMutationType blockMutType){
    int totalMutations = 0;

    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->nucMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++){
            
            if((root->nucMutation[i].condensed & 0x7) == nucMutType){
                // init += (((root->nucMutation[i].condensed) >> 3) & 0x1F);
                if(nucMutType == PangenomeMAT::NucMutationType::NS){
                    init += (((root->nucMutation[i].condensed) >> 3) & 0x1F); // Length of contiguous mutation in case of substitution
                } else {
                    init++;
                }
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
            if((nucMutation.condensed & 0x7) == nucMutType){
                if(nucMutType == PangenomeMAT::NucMutationType::NS){
                    totalMutations += ((nucMutation.condensed >> 3) & 0x1F); // Length of contiguous mutation in case of substitution
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
    std::cout << "Total SNP Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPS) << std::endl;
    std::cout << "Total SNP Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPI) << std::endl;
    std::cout << "Total SNP Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPD) << std::endl;
    std::cout << "Max Tree Depth: " << m_maxDepth << std::endl;
    std::cout << "Mean Tree Depth: " << m_meanDepth << std::endl;

    std::cout << "\nParallel Results:\n";
    std::cout << "Total Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NI, PangenomeMAT::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::ND, PangenomeMAT::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPS) << std::endl;
    std::cout << "Total SNP Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPI) << std::endl;
    std::cout << "Total SNP Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPD) << std::endl;

#else
    std::cout << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    std::cout << "Total Samples in Tree: " << m_numLeaves << std::endl;
    std::cout << "Total Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NS) << std::endl;
    std::cout << "Total Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NI, PangenomeMAT::BlockMutationType::BI) << std::endl;
    std::cout << "Total Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::ND, PangenomeMAT::BlockMutationType::BD) << std::endl;
    std::cout << "Total SNP Substitutions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPS) << std::endl;
    std::cout << "Total SNP Insertions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPI) << std::endl;
    std::cout << "Total SNP Deletions: " << getTotalParsimonyParallel(PangenomeMAT::NucMutationType::NSNPD) << std::endl;
    std::cout << "Max Tree Depth: " << m_maxDepth << std::endl;
    std::cout << "Mean Tree Depth: " << m_meanDepth << std::endl;
#endif

}

void PangenomeMAT::Tree::printBfs(){
    // Traversal test
    std::queue<Node *> bfsQueue;
    size_t prevLev = 0;
    
    bfsQueue.push(root);
    // bfsQueue.push(allNodes["node_10"]);

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

void printSequenceLines(const std::vector< std::vector< std::pair< char, std::vector< char > > > >& sequence,\
    const std::vector<bool>& blockExists, size_t lineSize, bool aligned, std::ofstream& fout){

    // std::vector< std::string > lines;
    std::string line;

    for(size_t i = 1; i < blockExists.size(); i++){
        if(!blockExists[i]){
            if(aligned){
                for(size_t j = 0; j < sequence[i].size(); j++){
                    for(size_t k = 0; k < sequence[i][j].second.size(); k++){
                        line += '-';
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                    if(sequence[i][j].first != 'x'){
                        line+='-';
                        if(line.length() == lineSize){
                            fout << line << '\n';
                            line = "";
                        }
                    }
                }
            }

            continue;
        }

        for(size_t j = 0; j < sequence[i].size(); j++){
            // std::cout << "A" << std::endl;

            for(size_t k = 0; k < sequence[i][j].second.size(); k++){

                if(sequence[i][j].second[k] == '-'){
                    if(aligned){
                        line += '-';
                    }
                    
                } else {
                    // if(sequence[i][j].second[k] != 'A' && sequence[i][j].second[k] != 'C' && sequence[i][j].second[k] != 'G' && sequence[i][j].second[k] != 'T' && sequence[i][j].second[k] != 'N'){
                    //     std::cout << sequence[i][j].second[k] << std::endl;
                    // }
                    line += sequence[i][j].second[k];
                }
                if(line.length() == lineSize){
                    fout << line << '\n';
                    line = "";
                }

            }

            // std::cout << "D" << std::endl;

            if(sequence[i][j].first != 'x'){
                // std::cout << "E" << std::endl;

                if(sequence[i][j].first == '-'){
                    // std::cout << "F" << std::endl;
                    if(aligned){
                        line += '-';
                    }
                } else {
                    // if(sequence[i][j].first != 'A' && sequence[i][j].first != 'C' && sequence[i][j].first != 'G' && sequence[i][j].first != 'T' && sequence[i][j].first != 'N'){
                    //     std::cout << sequence[i][j].first << std::endl;
                    // }

                    // std::cout << "G" << std::endl;
                    line += sequence[i][j].first;
                    // std::cout << "G" << std::endl;
                }
                if(line.length() == lineSize){
                    // std::cout << "H" << std::endl;
                    fout << line << '\n';
                    // std::cout << "H" << std::endl;
                    line = "";
                    // std::cout << "H\n" << std::endl;
                }
            }

            // std::cout << "B" << std::endl;
        }

    }

    // std::cout << "C" << std::endl;

    if(line.length()){
        fout << line << '\n';
        line = "";
    }

    fout << '\n';

}

char getNucleotideFromCode(int code){
    switch(code){
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        default:
            return 'N';
    }
}

void printFASTAHelper(PangenomeMAT::Node* root,\
    std::vector< std::vector< std::pair< char, std::vector< char > > > >& sequence, std::vector<bool>& blockExists,\
    std::ofstream& fout){
    
    // Apply mutations
    // Block mutations - ignored for now since the block IDs don't seem right in the files

    std::vector< std::tuple< int, bool, bool > > blockMutationInfo;

    for(uint32_t mutation: root->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & ((1 << 24) - 1));
        int type = (mutation & 0x1);

        // std::cout << mutation << " " << bid  << " " << sequence.size() << " " << type << std::endl;

        if(type == 0){
            bool oldVal = blockExists[bid];
            blockExists[bid] = true;
            blockMutationInfo.push_back( std::make_tuple(bid, oldVal, true) );
        } else {
            bool oldVal = blockExists[bid];
            blockExists[bid] = false;
            blockMutationInfo.push_back( std::make_tuple(bid, oldVal, false) );
        }

    }

    // For backtracking. blockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++){

        int bid = ((root->nucMutation[i].condensed >> 8) & (((1 << 24) - 1)));

        int pos = root->nucMutation[i].position;
        // if(pos < sequence[bid].size()){
        //     std::cout << pos << " " << sequence[bid].size() << std::endl;
        // }

        int gapPos = root->nucMutation[i].gapPosition;
        int type = ((root->nucMutation[i].condensed) & 7);
        char newVal = '-';

        if(type < 3){
            // Either S, I or D

            int len = (((root->nucMutation[i].condensed) >> 3) & 0x1F);

            if(type == PangenomeMAT::NucMutationType::NS){
                for(int j = 0; j < len; j++){
                    char oldVal = sequence[bid][pos + j].first;
                    newVal = getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15);

                    sequence[bid][pos + j].first = newVal;
                    mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, oldVal, newVal));
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NI){
                
                // if(gapPos == -1){
                //     for(int j = 0; j < len; j++){
                //         newVal = getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15);
                //         if(bid > sequence.size() || pos + j > sequence[bid].size()){
                //             std::cout << bid << " " << sequence.size() << std::endl;
                //         }
                //         sequence[bid][pos + j].first = newVal;
                //         mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, '-', newVal));
                //     }
                // }
                // else {
                //     for(int j = 0; j < len; j++){
                //         newVal = getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(7-j))) & 15);

                //         sequence[bid][pos].second[gapPos + j] = newVal;
                //         mutationInfo.push_back(std::make_tuple(bid, pos, gapPos + j, '-', newVal));
                //     }
                // }
            }
            // else if(type == PangenomeMAT::NucMutationType::ND){
            //     if(gapPos == -1){
            //         for(int j = 0; j < len; j++){
            //             char oldVal = sequence[bid][pos + j].first;
            //             sequence[bid][pos + j].first = '-';
            //             mutationInfo.push_back(std::make_tuple(bid, pos + j, -1, oldVal, '-'));
            //         }
            //     } else {
            //         for(int j = 0; j < len; j++){
            //             char oldVal = sequence[bid][pos].second[gapPos + j];
            //             sequence[bid][pos].second[gapPos + j] = '-';
            //             mutationInfo.push_back(std::make_tuple(bid, pos, gapPos + j, oldVal, '-'));
            //         }
            //     }
            // }
        } 
    //     else {
    //         // Todo: SNP
    //         if(type == PangenomeMAT::NucMutationType::NSNPS){

    //             newVal = getNucleotideFromCode(((root->nucMutation[i].condensed) >> 3) & 0xF);
    //             char oldVal = sequence[bid][pos].first;

    //             sequence[bid][pos].first = newVal;

    //             mutationInfo.push_back(std::make_tuple(bid, pos, -1, oldVal, newVal));
    //         }
            
    //         // else if(type == PangenomeMAT::NucMutationType::NSNPI){
    //         //     newVal = getNucleotideFromCode(((root->nucMutation[i].condensed) >> 3) & 0xF);
    //         //     if(gapPos == -1){
    //         //         char oldVal = sequence[bid][pos].first;
    //         //         sequence[bid][pos].first = newVal;
    //         //         mutationInfo.push_back(std::make_tuple(bid, pos, -1, oldVal, newVal));
    //         //     } else {
    //         //         char oldVal = sequence[bid][pos].second[gapPos];
    //         //         sequence[bid][pos].second[gapPos] = newVal;
    //         //         mutationInfo.push_back(std::make_tuple(bid, pos, gapPos, oldVal, newVal));
    //         //     }
    //         // }
    //         else if(type == PangenomeMAT::NucMutationType::NSNPD){
    //             if(gapPos == -1){

    //                 char oldVal = sequence[bid][pos].first;

    //                 sequence[bid][pos].first = '-';
    //                 mutationInfo.push_back(std::make_tuple(bid, pos, -1, oldVal, '-'));
    //             } else {
    //                 char oldVal = sequence[bid][pos].second[gapPos];
    //                 sequence[bid][pos].second[gapPos] = '-';
    //                 mutationInfo.push_back(std::make_tuple(bid, pos, gapPos, oldVal, '-'));
    //             }
    //         }
    //     }
    }

    if(root->children.size() == 0){
        // Print sequence
        fout << '>' << root->identifier << std::endl;

        printSequenceLines(sequence, blockExists, 50, false, fout);

    } else {
        // DFS on children
        for(PangenomeMAT::Node* child: root->children){
            printFASTAHelper(child, sequence, blockExists, fout);
        }
    }

    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++){
        auto mutation = *it;
        blockExists[std::get<0>(mutation)] = std::get<1>(mutation);
    }

    // Undo mutations
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++){
        auto mutation = *it;
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

    for(size_t i = 0; i < blocks.size(); i++){
        
        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
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
                    default:
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
        // int bId = ((gaps.condensed[i] >> 8) & (0xFFFFFF));
        int bId = gaps.blockId[i];

        // int len = ((gaps.condensed[i]) & 255);
        int len = gaps.gapLength[i];
        int pos = gaps.position[i];

        // if(bId >= sequence.size() || pos >= sequence[bId].size()){
        //     std::cout << "Index: " << i << " Bid: " << bId << " Total Blocks: " << sequence.size() << std::endl;
        // }

        sequence[bId][pos].second.resize(len, '-');
    }

    printFASTAHelper(root, sequence, blockExists, fout);

}

// Merge parent node and child node into parent node
void mergeNodes(PangenomeMAT::Node* par, PangenomeMAT::Node* chi){
    
    par->identifier = chi->identifier;
    par->branchLength += chi->branchLength;
    par->children = chi->children;

    // For block mutations, we cancel out irrelevant mutations
    std::unordered_set< int > bidInserts;

    for(auto mutation: par->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & 0xFFFFFF);
        int type = (mutation & 0x1);
        if(type == PangenomeMAT::BlockMutationType::BI){
            bidInserts.insert(bid);
        } else {
            bidInserts.erase(bid);
        }
    }

    for(auto mutation: chi->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & 0xFFFFFF);
        int type = (mutation & 0x1);
        if(type == PangenomeMAT::BlockMutationType::BI){
            bidInserts.insert(bid);
        } else {
            bidInserts.erase(bid);
        }
    }

    PangenomeMAT::BlockMut newBlockMutation;
    for(auto bid: bidInserts){
        newBlockMutation.condensedBlockMut.push_back((( bid << 8 ) ^ 0x1));
    }

    par->blockMutation = newBlockMutation;

    for(auto mutation: chi->nucMutation){
        par->nucMutation.push_back(mutation);
    }

    delete chi;
}

void compressTree(PangenomeMAT::Node* node){
    if(node->children.size() == 0){
        return;
    }

    for(size_t i = 0; i < node->children.size(); i++){
        while(node->children[i]->children.size() == 1){
            mergeNodes(node->children[i], node->children[i]->children[0]);
        }
        compressTree(node->children[i]);
    }
}

PangenomeMAT::Node* subtreeExtractHelper(PangenomeMAT::Node* node, const std::unordered_map< PangenomeMAT::Node*, size_t >& ticks){
    if(ticks.find(node) == ticks.end()){
        return nullptr;
    }

    PangenomeMAT::Node* newNode = new PangenomeMAT::Node(node->identifier, node->branchLength);
    
    // We don't care about level anymore since that is just preprocessing for the summary
    newNode->level = -1;

    for(auto mutation: node->nucMutation){
        newNode->nucMutation.push_back(mutation);
    }

    for(auto mutation: node->blockMutation.condensedBlockMut){
        newNode->blockMutation.condensedBlockMut.push_back(mutation);
    }

    for(auto child: node->children){
        if(ticks.find(child) != ticks.end()){
            PangenomeMAT::Node* newChild = subtreeExtractHelper(child, ticks);
            newChild->parent = newNode;
            newNode->children.push_back(newChild);
        }
    }

    return newNode;

}

PangenomeMAT::Node* PangenomeMAT::Tree::subtreeExtract(std::vector< PangenomeMAT::Node* > requiredNodes){

    std::unordered_map< PangenomeMAT::Node*, size_t > ticks;
    for(auto node: requiredNodes){
        Node* current = node;

        while(current != nullptr){
            ticks[current]++;
            current = current->parent;
        }
    }

    return nullptr;

}

void PangenomeMAT::Tree::printFASTA_updated(std::ofstream& fout){

    // Categorize gaps by blockId
    std::map< int, std::vector< std::pair< int, int > > > gapSplit;
    
    for(size_t i = 0; i < gaps.position.size(); i++){

        // int bId = ((gaps.condensed[i] >> 8) & (0xFFFFFF));
        int bId = gaps.blockId[i];


        // int len = ((gaps.condensed[i]) & 255);
        int len = gaps.gapLength[i];

        int pos = gaps.position[i];
        gapSplit[bId].push_back( std::make_pair(pos, len) );

    }
    int lineCount = 1;
    int charCount = 1;

    // Get path from leaf to root
    for(auto leaf: allLeaves){
        std::vector< PangenomeMAT::Node* > path;
        Node* it = leaf;
        while(it != root){
            path.push_back(it);
            it = it->parent;
        }
        path.push_back(root);

        // Get all blocks on the path
        std::unordered_set< uint32_t > blockIds;
        for(auto node = path.rbegin(); node != path.rend(); node++){
            for(uint32_t mutation: (*node)->blockMutation.condensedBlockMut){
                int bid = ((mutation >> 8) & 0xFFFFFF);
                int type = (mutation & 0x1);

                if(type == PangenomeMAT::BlockMutationType::BI){
                    blockIds.insert(bid);
                } else {
                    blockIds.erase(bid);
                }
            }
        }

        // Create the required blocks
        std::map< int, std::vector< std::pair< char, std::vector< char > > > > sequence;
        for(auto bid: blockIds){
            if(blocks[bid - 1].blockId != bid){
                std::cout << "Block not in correct position in blocks array" << std::endl;
            }

            for(size_t i = 0; i < blocks[bid - 1].consensusSeq.size(); i++){
                for(int j = 0; j < 8; j++){
                    const int nucCode = (((blocks[bid - 1].consensusSeq[i]) >> (4*(7 - j))) & 0xF);
                    sequence[bid].push_back({ getNucleotideFromCode(nucCode), {} });
                }
            }
            sequence[bid].push_back({'x',{}});

            for(auto g: gapSplit[bid]){
                sequence[bid][g.first].second.resize(g.second, '-');
            }

        }

        // Apply nucleotide mutations
        for(auto node = path.rbegin(); node != path.rend(); node++){

            for(size_t i = 0; i < (*node)->nucMutation.size(); i++){

                int bid = (((*node)->nucMutation[i].condensed >> 8) & 0xFFFFFF);

                if(sequence.find(bid) == sequence.end()){
                    continue;
                }

                int pos = (*node)->nucMutation[i].position;
                int gapPos = (*node)->nucMutation[i].gapPosition;
                int type = (((*node)->nucMutation[i].condensed) & 7);
                char newVal = '-';

                

                if(type < 3){
                    // Either S, I or D

                    int len = ((((*node)->nucMutation[i].condensed) >> 3) & 0x1F);

                    if(type == PangenomeMAT::NucMutationType::NS){
                        for(int j = 0; j < len; j++){
                            newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(7-j))) & 15);
                            sequence[bid][pos + j].first = newVal;
                        }
                    } else if(type == PangenomeMAT::NucMutationType::NI){
                        
                        if(gapPos == -1){
                            for(int j = 0; j < len; j++){
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(7-j))) & 15);
                                
                                sequence[bid][pos + j].first = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(7-j))) & 15);

                                sequence[bid][pos].second[gapPos + j] = newVal;
                            }
                        }
                    } else if(type == PangenomeMAT::NucMutationType::ND){
                        if(gapPos == -1){
                            for(int j = 0; j < len; j++){
                                sequence[bid][pos + j].first = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++){
                                sequence[bid][pos].second[gapPos + j] = '-';
                            }
                        }
                    }
                } else {
                    // Todo: SNP
                    if(type == PangenomeMAT::NucMutationType::NSNPS){
                        newVal = getNucleotideFromCode((((*node)->nucMutation[i].condensed) >> 3) & 0xF);
                        sequence[bid][pos].first = newVal;
                    }
                    // else if(type == PangenomeMAT::NucMutationType::NSNPI){
                    //     newVal = getNucleotideFromCode(((root->nucMutation[i].condensed) >> 3) & 0xF);
                    //     if(gapPos == -1){
                    //         sequence[bid][pos].first = newVal;
                    //         // mutationInfo.push_back(std::make_tuple(bid, pos, -1, oldVal, newVal));
                    //     } else {
                    //         if((size_t)gapPos >= sequence[bid][pos].second.size()){
                    //             std::cout << (*node)->identifier << " " << i << " " << bid << " " << pos << " " << gapPos << " " << sequence[bid][pos].second.size() << std::endl;
                    //         }
                    //         sequence[bid][pos].second[gapPos] = newVal;
                    //         // mutationInfo.push_back(std::make_tuple(bid, pos, gapPos, oldVal, newVal));
                    //         if((size_t)gapPos >= sequence[bid][pos].second.size()){
                    //             std::cout << bid << " " << pos << " " << gapPos << " " << sequence[bid][pos].second.size() << std::endl;
                    //         }
                    //     }
                    // }
                    else if(type == PangenomeMAT::NucMutationType::NSNPD){
                        if(gapPos == -1){
                            sequence[bid][pos].first = '-';
                        } else {
                            sequence[bid][pos].second[gapPos] = '-';
                        }
                    }
                }
            }


        }

        // Print the sequence
        fout << '>' << leaf->identifier << '\n';
        lineCount++;
        std::string line;
        for(const auto &[bid, block]: sequence){

            for(size_t i = 0; i < block.size(); i++){
                for(size_t j = 0; j < block[i].second.size(); j++){
                    if(block[i].second[j] != '-'){

                        line += block[i].second[j];
                        charCount++;
                        if(line.size() == 50){
                            fout << line << '\n';
                            line = "";
                            lineCount++;
                            charCount=1;
                        }
                    }
                }
                if(block[i].first != 'x' && block[i].first != '-'){

                    line += block[i].first;
                    charCount++;
                    if(line.size() == 50){
                        fout << line << '\n';
                        line = "";
                        lineCount++;
                        charCount=1;
                    }
                }
            }
        }
        if(line.length()){
            fout << line << '\n';
            line = "";
            lineCount++;
            charCount=1;
        }
        fout << '\n';
        lineCount++;
        charCount=1;
    }

}

void getNodesPreorder(PangenomeMAT::Node* root, MAT::tree& treeToWrite){
    
    MAT::node n;
    
    MAT::block_mut bm;
    for(auto mutation: root->blockMutation.condensedBlockMut){
        bm.add_condensed_block_mut(mutation);
    }

    *n.mutable_block_mutation() = bm;

    for(size_t i = 0; i < root->nucMutation.size(); i++){
        const PangenomeMAT::NucMut& mutation = root->nucMutation[i];

        MAT::nuc_mut nm;
        nm.set_position(mutation.position);
        if(mutation.gapPosition != -1){
            nm.set_gap_position(mutation.gapPosition);
        }
        nm.set_condensed(mutation.condensed);
        nm.set_nucs(mutation.nucs);

        n.add_nuc_mutation();
        *n.mutable_nuc_mutation(i) = nm;
    }

    treeToWrite.add_nodes();
    *treeToWrite.mutable_nodes( treeToWrite.nodes_size() - 1 ) = n;

    for(auto child: root->children){
        getNodesPreorder(child, treeToWrite);
    }
}

void PangenomeMAT::Tree::writeToFile(std::ofstream& fout){
    MAT::tree treeToWrite;
    getNodesPreorder(root, treeToWrite);

    for(auto block: blocks){
        MAT::block b;
        b.set_block_id(block.blockId);
        b.set_chromosome_name(block.chromosomeName);
        for(auto n: block.consensusSeq){
            b.add_consensus_seq(n);
        }
        treeToWrite.add_blocks();
        *treeToWrite.mutable_blocks( treeToWrite.blocks_size() - 1 ) = b;
    }

    MAT::gap_list gl;
    for(size_t i = 0; i < gaps.position.size(); i++){
        gl.add_position(gaps.position[i]);
        gl.add_block_id(gaps.blockId[i]);
        gl.add_gap_length(gaps.gapLength[i]);
    }

    *treeToWrite.mutable_gaps() = gl;

    if (!treeToWrite.SerializeToOstream(&fout)) {
		std::cerr << "Failed to write output file." << std::endl;
    }

}

void sampleWriteToFileHelper(PangenomeMAT::Node* root, int& currentNode, std::ofstream& fout){
    fout << "node " << currentNode << ":\n";
    fout << "nuc_mutation:\n";
    for(auto mutation: root->nucMutation){
        fout << "(" << mutation.position << ", " << (mutation.gapPosition == -1? 0: mutation.gapPosition) << ", " << mutation.condensed << ", " << mutation.nucs << ") ";
    }
    fout << "\n";
	fout << "block_mutation:\n";

    for(auto mutation: root->blockMutation.condensedBlockMut){
        fout << mutation << " ";
    }
    fout << "\n";
    for(auto child: root->children){
        currentNode++;
        sampleWriteToFileHelper(child, currentNode, fout);
    }

}

void PangenomeMAT::Tree::sampleWriteToFile(std::ofstream& fout){
    fout << "nodes:\n";
    int initialNode = 0;
    sampleWriteToFileHelper(root, initialNode, fout);
}