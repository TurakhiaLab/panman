#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <ctime>
#include <iomanip>
#include <mutex>
#include <chrono>

#include "chaining.cpp"
#include "rotation.cpp"

#include "PangenomeMAT.hpp"

KSEQ_INIT(int, read)

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

char PangenomeMAT::getNucleotideFromCode(int code) {
    switch(code) {
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

char PangenomeMAT::getCodeFromNucleotide(char nuc) {
    switch(nuc) {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 4;
        case 'T':
            return 8;
        case 'R':
            return 5;
        case 'Y':
            return 10;
        case 'S':
            return 6;
        case 'W':
            return 9;
        case 'K':
            return 12;
        case 'M':
            return 3;
        case 'B':
            return 14;
        case 'D':
            return 13;
        case 'H':
            return 11;
        case 'V':
            return 7;
        case 'N':
            return 15;
        default:
            return 0;
    }
}

// For reverse complement
char PangenomeMAT::getComplementCharacter(char nuc) {
    switch(nuc) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 'R':
            return 'Y';
        case 'Y':
            return 'R';
        case 'S':
            return 'S';
        case 'W':
            return 'W';
        case 'K':
            return 'M';
        case 'M':
            return 'K';
        case 'B':
            return 'V';
        case 'D':
            return 'H';
        case 'H':
            return 'D';
        case 'V':
            return 'B';
        default:
            return 'N';
    }
}

std::string PangenomeMAT::getDate() {
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);
    std::string date;
    date += std::to_string(now->tm_year + 1900)
         + std::to_string(now->tm_mon + 1)
         +  std::to_string(now->tm_mday);
    return date;
}

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

int PangenomeMAT::Tree::nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states) {
    if(node->children.size() == 0) {
        if(states.find(node->identifier) == states.end()) {
            return states[node->identifier] = 0;
        }
        return states[node->identifier];
    }
    std::vector< int > childStates;
    for(auto child: node->children) {
        childStates.push_back(nucFitchForwardPass(child, states));
    }
    int orStates = 0, andStates = childStates[0];
    for(auto u: childStates) {
        orStates |= u;
        andStates &= u;
    }
    if(andStates) {
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;
}

void PangenomeMAT::Tree::nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultState) {
    if(node == root && defaultState != (1 << 28)) {
        states[node->identifier] = defaultState;
    } else {
        if(states[node->identifier] == 0) {
            return;
        }
        if(parentState & states[node->identifier]) {
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children) {
        nucFitchBackwardPass(child, states, states[node->identifier]);
    }
}

void PangenomeMAT::Tree::nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT::NucMutationType, char > >& mutations, int parentState) {
    if(states[node->identifier] == 0) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 1) {
            // insertion
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 1) {
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }
    for(auto child: node->children) {
        nucFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}

int PangenomeMAT::Tree::blockFitchForwardPassNew(Node* node, std::unordered_map< std::string, int >& states) {
    if(node->children.size() == 0) {
        if(states.find(node->identifier) == states.end()) {
            return states[node->identifier] = 0;
        }
        return states[node->identifier];
    }
    std::vector< int > childStates;
    for(auto child: node->children) {
        childStates.push_back(blockFitchForwardPassNew(child, states));
    }
    int orStates = 0, andStates = childStates[0];
    for(auto u: childStates) {
        orStates |= u;
        andStates &= u;
    }
    if(andStates) {
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;

    // if(node->children.size() == 0) {
    //     if(states.find(node->identifier) == states.end()) {
    //         std::cerr << "FATAL: State for leaf Node ID not found!" << std::endl;
    //         exit(-1);
    //         return 0;
    //     }
    //     return states[node->identifier];
    // }
    // std::vector< bool > stateExists(3,false);
    // for(auto child: node->children) {
    //     stateExists[blockFitchForwardPass(child, states)] = true;
    // }
    // if(stateExists[0] && stateExists[1]) {
    //     return states[node->identifier] = 2;
    // }
    // if(stateExists[0]) {
    //     return states[node->identifier] = 0;
    // }
    // if(stateExists[1]) {
    //     return states[node->identifier] = 1;
    // }
    // return states[node->identifier] = 2;
}

void PangenomeMAT::Tree::blockFitchBackwardPassNew(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultValue) {
    if(node == root && defaultValue != (1 << 28)) {
        states[node->identifier] = defaultValue;
    } else {
        if(states[node->identifier] == 0) {
            return;
        }
        if(parentState & states[node->identifier]) {
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children) {
        blockFitchBackwardPassNew(child, states, states[node->identifier]);
    }

}

void PangenomeMAT::Tree::blockFitchAssignMutationsNew(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT::BlockMutationType, bool > >& mutations, int parentState) {
    if(states[node->identifier] == 0) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 1) {
            // insertion
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;
            if(code == 2) {
                // insertion of inverted block
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, true);
            } else {
                // insertion of forward strand
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, false);
            }

        } else if(states[node->identifier] == 1) {
            // deletion
            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, false);
        } else {
            // inversion

            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, true);
        }
    }
    for(auto child: node->children) {
        blockFitchAssignMutationsNew(child, states, mutations, states[node->identifier]);
    }
}

int PangenomeMAT::Tree::blockFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states) {
    if(node->children.size() == 0) {
        if(states.find(node->identifier) == states.end()) {
            std::cerr << "FATAL: State for leaf Node ID not found!" << std::endl;
            exit(-1);
            return 0;
        }
        return states[node->identifier];
    }
    std::vector< bool > stateExists(3,false);
    for(auto child: node->children) {
        stateExists[blockFitchForwardPass(child, states)] = true;
    }
    if(stateExists[0] && stateExists[1]) {
        return states[node->identifier] = 2;
    }
    if(stateExists[0]) {
        return states[node->identifier] = 0;
    }
    if(stateExists[1]) {
        return states[node->identifier] = 1;
    }
    return states[node->identifier] = 2;
}

void PangenomeMAT::Tree::blockFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultValue) {
    if(defaultValue != 2) {
        states[node->identifier] = defaultValue;
    } else if(states[node->identifier] == 2) {
        states[node->identifier] = parentState;
    }
    for(auto child: node->children) {
        blockFitchBackwardPass(child, states, states[node->identifier]);
    }
}

void PangenomeMAT::Tree::blockFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, bool >& mutations, int parentState) {
    if(parentState == 0 && states[node->identifier] == 1) {
        mutations[node->identifier] = true;
    } else if(parentState == 1 && states[node->identifier] == 0) {
        mutations[node->identifier] = false;
    }
    for(auto child: node->children) {
        blockFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}

PangenomeMAT::Tree::Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype) {

    if(ftype == PangenomeMAT::FILE_TYPE::GFA) {
        std::map< std::string, std::string > nodes;
        std::map< std::string, std::vector< std::string > > paths;

        std::string line;
        while(getline(fin, line, '\n')) {
            std::vector< std::string > separatedLine;
            stringSplit(line, '\t', separatedLine);
            if(separatedLine[0] == "S") {
                nodes[separatedLine[1]] = separatedLine[2];
            } else if(separatedLine[0] == "P") {
                std::vector< std::string > v;
                stringSplit(separatedLine[2], ',', v);
                for(size_t i = 0; i < v.size(); i++) {
                    v[i].pop_back();
                }
                paths[separatedLine[1]] = v;
            }
        }
        std::vector< std::vector< std::string > > stringSequences;
        std::vector< std::string > sequenceIds;
        for(auto p: paths) {
            sequenceIds.push_back(p.first);
            stringSequences.push_back(p.second);
        }
        
        GFAGraph g(sequenceIds, stringSequences, nodes);
        std::cout << "Graph without cycles created" << std::endl;

        std::vector< size_t > topoArray = g.getTopologicalSort();
        std::vector< std::vector< int64_t > > alignedSequences = g.getAlignedSequences(topoArray);
        std::string newickString;
        secondFin >> newickString;
        root = createTreeFromNewickString(newickString);

        std::unordered_map< std::string, std::vector< int64_t > > pathIdToSequence;
        for(size_t i = 0; i < g.pathIds.size(); i++) {
            pathIdToSequence[g.pathIds[i]] = alignedSequences[i];
        }

        for(size_t i = 0; i < topoArray.size(); i++) {
            blocks.emplace_back(i, g.intNodeToSequence[topoArray[i]]);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, bool > > globalMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, bool > mutations;
            for(const auto& u: pathIdToSequence) {
                states[u.first] = (u.second[i] != -1);
            }
            blockFitchForwardPass(root, states);
            blockFitchBackwardPass(root, states, 1);
            blockFitchAssignMutations(root, states, mutations, 0);
            globalMutations[i] = mutations;
        });

        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes) {
            nodeMutexes[u.first];
        }

        tbb::parallel_for_each(globalMutations, [&](auto& pos) {
            auto& mutations = pos.second;
            for(const auto& node: allNodes) {
                if(mutations.find(node.first) != mutations.end()) {
                    nodeMutexes[node.first].lock();
                    node.second->blockMutation.emplace_back(pos.first, mutations[node.first]);
                    nodeMutexes[node.first].unlock();
                }
            }
        });

    } else if(ftype == PangenomeMAT::FILE_TYPE::PANGRAPH) {
        std::string newickString;
        secondFin >> newickString;
        Json::Value pangraphData;
        fin >> pangraphData;

        auto start = std::chrono::high_resolution_clock::now();

        PangenomeMAT::Pangraph pg(pangraphData);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds timing = end -start;
        std::cout << "Psuedo Root Calculated in: " << timing.count() << " nanoseconds \n";
        
        circularSequences = pg.circularSequences;
        sequenceInverted = pg.sequenceInverted;
        rotationIndexes = pg.rotationIndexes;

        std::vector< size_t > topoArray = pg.getTopologicalSort();
        std::cout << "Length of Pseudo Root: " << topoArray.size() << std::endl;
        std::unordered_map< std::string, std::vector< int > > alignedSequences = pg.getAlignedSequences(topoArray);
        std::unordered_map< std::string, std::vector< int > > alignedStrandSequences = pg.getAlignedStrandSequences(topoArray);
        
        root = createTreeFromNewickString(newickString);

        for(size_t i = 0; i < topoArray.size(); i++) {
            blocks.emplace_back(i, pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]]);
        }

        for(size_t i = 0; i < topoArray.size(); i++) {
            GapList g;
            g.primaryBlockId = i;
            g.secondaryBlockId = -1;
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++) {
                g.nucPosition.push_back(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first);
                g.nucGapLength.push_back(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second);
            }
            gaps.push_back(g);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, std::pair< BlockMutationType, bool > > > globalBlockMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;
            for(const auto& u: alignedSequences) {
                if(u.second[i] == -1) {
                    states[u.first] = 1;
                } else if(alignedStrandSequences[u.first][i]) {
                    // forward strand
                    states[u.first] = 2;
                } else {
                    // reverse strand
                    states[u.first] = 4;
                }
            }
            
            blockFitchForwardPassNew(root, states);
            blockFitchBackwardPassNew(root, states, 1);
            blockFitchAssignMutationsNew(root, states, mutations, 1);
            globalBlockMutations[i] = mutations;

        });

        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes) {
            nodeMutexes[u.first];
        }

        tbb::parallel_for_each(globalBlockMutations, [&](auto& pos) {
            auto& mutations = pos.second;
            for(const auto& node: allNodes) {
                if(mutations.find(node.first) != mutations.end()) {
                    nodeMutexes[node.first].lock();
                    node.second->blockMutation.emplace_back(pos.first, mutations[node.first]);
                    nodeMutexes[node.first].unlock();
                }
            }
        });

        std::unordered_map< std::string, std::vector< size_t > > blockCounts;
        for(const auto& u: alignedSequences) {
            blockCounts[u.first].resize(u.second.size(), 0);
        }

        tbb::parallel_for_each(alignedSequences, [&](const auto& u) {
            // std::unordered_map< std::string, size_t > currentCount;
            int currentPtr = 0;
            for(size_t i = 0; i < u.second.size(); i++) {
                if(u.second[i] != -1) {
                    blockCounts[u.first][i] = pg.blockNumbers[u.first][currentPtr];
                    currentPtr++;
                }
            }
        });

        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > gapMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            std::string consensusSeq = pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]];
            std::vector< std::pair< char, std::vector< char > > > sequence(consensusSeq.size()+1, {'-', {}});
            for(size_t j = 0; j < consensusSeq.length(); j++) {
                sequence[j].first = consensusSeq[j];
            }
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++) {
                sequence[pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first].second.resize(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second, '-');
            }
            tbb::concurrent_unordered_map< std::string, std::vector< std::pair< char, std::vector< char > > > > individualSequences;

            tbb::parallel_for_each(alignedSequences, [&](const auto& u) {
                if(u.second[i] == -1) {
                    return;
                }
                std::vector< std::pair< char, std::vector< char > > > currentSequence = sequence;

                for(const auto& v: pg.substitutions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]) {
                    currentSequence[v.first-1].first = v.second[0];
                }
                for(const auto& v: pg.insertions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]) {
                    for(size_t j = 0; j < std::get<2>(v).length(); j++) {
                        currentSequence[std::get<0>(v)].second[std::get<1>(v)+j] = std::get<2>(v)[j];
                    }
                }
                for(const auto& v: pg.deletions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]) {
                    for(size_t j = v.first; j < v.first + v.second; j++) {
                        currentSequence[j-1].first = '-';
                    }
                }
                individualSequences[u.first] = currentSequence;
            });
            tbb::parallel_for((size_t) 0, sequence.size(), [&](size_t j) {
                tbb::parallel_for((size_t)0, sequence[j].second.size(), [&](size_t k) {
                    std::unordered_map< std::string, int > states;
                    std::unordered_map< std::string, std::pair< PangenomeMAT::NucMutationType, char > > mutations;
                    for(const auto& u: individualSequences) {
                        if(u.second[j].second[k] != '-') {
                            states[u.first] = (1 << getCodeFromNucleotide(u.second[j].second[k]));
                        } else {
                            states[u.first] = 1;
                        }
                    }
                    nucFitchForwardPass(root, states);
                    nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                    nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                    for(auto mutation: mutations) {
                        nodeMutexes[mutation.first].lock();
                        gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                        nodeMutexes[mutation.first].unlock();
                    }
                });
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< PangenomeMAT::NucMutationType, char > > mutations;
                for(const auto& u: individualSequences) {
                    if(u.second[j].first != '-') {
                        states[u.first] = (1 << getCodeFromNucleotide(u.second[j].first));
                    } else {
                        states[u.first] = 1;
                    }
                }
                nucFitchForwardPass(root, states);
                nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].first)));
                nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].first)));
                for(auto mutation: mutations) {
                    nodeMutexes[mutation.first].lock();
                    nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                    nodeMutexes[mutation.first].unlock();
                }
            });
        });

        tbb::parallel_for_each(nonGapMutations, [&](auto& u) {
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++) {
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])) {
                    nodeMutexes[u.first].lock();
                    allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, i);
                    nodeMutexes[u.first].unlock();
                    currentStart = i;
                    continue;
                }
            }
            nodeMutexes[u.first].lock();
            allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, u.second.size());
            nodeMutexes[u.first].unlock();
        });
        tbb::parallel_for_each(gapMutations, [&](auto& u) {
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++) {
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1]) || std::get<3>(u.second[i]) != std::get<3>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])) {
                    nodeMutexes[u.first].lock();
                    allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, i);
                    nodeMutexes[u.first].unlock();
                    currentStart = i;
                    continue;
                }
            }
            nodeMutexes[u.first].lock();
            allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, u.second.size());
            nodeMutexes[u.first].unlock();
        });

    } else if(ftype == PangenomeMAT::FILE_TYPE::MSA) {
        std::string newickString;
        secondFin >> newickString;
        root = createTreeFromNewickString(newickString);

        std::map< std::string, std::string > sequenceIdsToSequences;
        std::string line;
        std::string currentSequence, currentSequenceId;
        size_t lineLength = 0;
        std::string consensusSeq;
        while(getline(fin,line,'\n')) {
            if(line.length() == 0) {
                continue;
            }
            if(line[0] == '>') {
                if(currentSequence.length()) {
                    if(lineLength == 0) {
                        lineLength = currentSequence.length();
                    } else if(lineLength != currentSequence.length()) {
                        std::cerr << "Error: sequence lengths don't match! " << currentSequenceId << std::endl;
                        exit(-1);
                    }
                    sequenceIdsToSequences[currentSequenceId] = currentSequence;
                }
                std::vector< std::string > splitLine;
                stringSplit(line,' ',splitLine);
                currentSequenceId = splitLine[0].substr(1);
                currentSequence = "";
            } else {
                currentSequence += line;
            }
        }
        if(currentSequence.length()) {
            if(lineLength != 0 && lineLength != currentSequence.length()) {
                std::cerr << "Error: sequence lengths don't match!" << std::endl;
                exit(-1);
            } else {
                lineLength = currentSequence.length();
            }
            sequenceIdsToSequences[currentSequenceId] = currentSequence;
        }
        std::set< size_t > emptyPositions;

        for(size_t i = 0; i < lineLength; i++) {
            bool nonGapFound = false;
            for(auto u: sequenceIdsToSequences) {
                if(u.second[i] != '-') {
                    consensusSeq += u.second[i];
                    nonGapFound = true;
                    break;
                }
            }
            if(!nonGapFound) {
                emptyPositions.insert(i);
            }
        }
        for(auto& u: sequenceIdsToSequences) {
            std::string sequenceString;
            for(size_t i = 0; i < u.second.length(); i++) {
                if(emptyPositions.find(i) == emptyPositions.end()) {
                    sequenceString += u.second[i];
                }
            }
            u.second = sequenceString;
        }
        blocks.emplace_back(0, consensusSeq);
        root->blockMutation.emplace_back(0, true);

        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes) {
            nodeMutexes[u.first];
        }

        tbb::parallel_for((size_t)0, consensusSeq.length(), [&](size_t i) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< PangenomeMAT::NucMutationType, char > > mutations;
            for(const auto& u: sequenceIdsToSequences) {
                if(u.second[i] != '-') {
                    states[u.first] = (1 << getCodeFromNucleotide(u.second[i]));
                } else {
                    states[u.first] = 1;
                }
            }
            nucFitchForwardPass(root, states);
            nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(consensusSeq[i])));
            nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSeq[i])));
            for(auto mutation: mutations) {
                nodeMutexes[mutation.first].lock();
                nonGapMutations[mutation.first].push_back(std::make_tuple(0, -1, i, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                nodeMutexes[mutation.first].unlock();
            }
        });

        tbb::parallel_for_each(nonGapMutations, [&](auto& u) {
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++) {
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])) {
                    nodeMutexes[u.first].lock();
                    allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, i);
                    nodeMutexes[u.first].unlock();
                    currentStart = i;
                    continue;
                }
            }
            nodeMutexes[u.first].lock();
            allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, u.second.size());
            nodeMutexes[u.first].unlock();
        });
        

    }

    // Setting up global coordinates for fast retrieval
    setupGlobalCoordinates();
}

void PangenomeMAT::Tree::protoMATToTree(const PanMAT::tree& mainTree) {
    // Create tree
    root = createTreeFromNewickString(mainTree.newick());

    std::vector< MAT::node > storedNodes;
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
    for(int i = 0; i < mainTree.gaps().position_size(); i++){
        gaps.position.push_back(mainTree.gaps().position(i));
    }

    for(int i = 0; i < mainTree.gaps().block_id_size(); i++){
        gaps.blockId.push_back(mainTree.gaps().block_id(i));
    }
    for(int i = 0; i < mainTree.gaps().gap_length_size(); i++){
        gaps.gapLength.push_back(mainTree.gaps().gap_length(i));
    }

}

PangenomeMAT::Tree::Tree(const PanMAT::tree& mainTree) {
    protoMATToTree(mainTree);
    setupGlobalCoordinates();

}

PangenomeMAT::Tree::Tree(std::istream& fin, FILE_TYPE ftype) {

    if(ftype == PangenomeMAT::FILE_TYPE::PANMAT) {
        PanMAT::tree mainTree;

        if(!mainTree.ParseFromIstream(&fin)) {
            throw std::invalid_argument("Could not read tree from input file.");
        }

        protoMATToTree(mainTree);
        setupGlobalCoordinates();

    }
    
}

int getTotalParsimonyParallelHelper(PangenomeMAT::Node* root, PangenomeMAT::NucMutationType nucMutType, PangenomeMAT::BlockMutationType blockMutType) {
    int totalMutations = 0;

    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->nucMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++){
            
            if((root->nucMutation[i].condensed & 0x7) == nucMutType){
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

void PangenomeMAT::Tree::printSummary(){

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

}

void PangenomeMAT::Tree::printBfs(Node* node){
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

std::string PangenomeMAT::getCurrentFastaSequence(const sequence_t& sequence, const blockExists_t& blockExists, blockStrand_t& blockStrand, bool aligned) {

    // String that stores the sequence to be printed
    std::string line;

    for(size_t i = 0; i < blockExists.size(); i++) {
        // Iterate through gap blocks - NOT BEING USED CURRENTLY
        for(size_t j = 0; j < blockExists[i].second.size(); j++) {
            // If block exists. Otherwise add gaps if MSA is to be printed
            if(blockExists[i].second[j]) {
                // If forward strand, iterare in forward direction
                if(blockStrand[i].second[j]) {
                    // Main nucs
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        // Gap nucs
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                            if(sequence[i].second[j][k].second[w] != '-') {
                                line += sequence[i].second[j][k].second[w];
                            } else if(aligned) {
                                line += '-';
                            }
                        }
                        // Main Nuc
                        if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x') {
                            line += sequence[i].second[j][k].first;
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                } else {
                    // If reverse strand, iterate backwards
                    for(size_t k = sequence[i].second[j].size()-1; k+1 > 0; k--) {
                        // Main nuc
                        if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x') {
                            line += getComplementCharacter(sequence[i].second[j][k].first);
                        } else if(aligned) {
                            line += '-';
                        }
                        // Gap nucs
                        for(size_t w = sequence[i].second[j][k].second.size()-1; w+1 > 0; w--) {
                            if(sequence[i].second[j][k].second[w] != '-') {
                                line += getComplementCharacter(sequence[i].second[j][k].second[w]);
                            } else if(aligned) {
                                line += '-';
                            }
                        }

                    }
                }
            } else if(aligned) {
                for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                        line += '-';
                    }
                    line += '-';
                }
            }
        }

        // Non-gap block - the only type being used currently
        if(blockExists[i].first) {
            // If forward strand
            if(blockStrand[i].first) {
                // Iterate through main nucs
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    // Gap nucs
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                        if(sequence[i].first[j].second[k] != '-') {
                            line += sequence[i].first[j].second[k];
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                    // Main nuc
                    if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                        line += sequence[i].first[j].first;
                    } else if(aligned) {
                        line += '-';
                    }
                }
            } else {
                // If reverse strand, iterate backwards
                for(size_t j = sequence[i].first.size()-1; j+1 > 0; j--) {
                    // Main nuc first since we are iterating in reverse direction
                    if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                        line += getComplementCharacter(sequence[i].first[j].first);
                    } else if(aligned) {
                        line += '-';
                    }

                    // Gap nucs
                    for(size_t k = sequence[i].first[j].second.size()-1; k+1 > 0; k--) {
                        if(sequence[i].first[j].second[k] != '-') {
                            line += getComplementCharacter(sequence[i].first[j].second[k]);
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                }   
            }
        } else if(aligned) {
            // If aligned sequence is required, print gaps instead if block does not exist
            for(size_t j = 0; j < sequence[i].first.size(); j++) {
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                    line+='-';
                }
                line+='-';
            }
        }
    }

    return line;

}

void PangenomeMAT::printSequenceLines(const sequence_t& sequence,\
    const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize, bool aligned, std::ofstream& fout, int offset, bool debug) {

    // String that stores the sequence to be printed
    std::string line;

    for(size_t i = 0; i < blockExists.size(); i++) {
        // Iterate through gap blocks - NOT BEING USED CURRENTLY
        for(size_t j = 0; j < blockExists[i].second.size(); j++) {
            // If block exists. Otherwise add gaps if MSA is to be printed
            if(blockExists[i].second[j]) {
                // If forward strand, iterare in forward direction
                if(blockStrand[i].second[j]) {
                    // Main nucs
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        // Gap nucs
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                            if(sequence[i].second[j][k].second[w] != '-') {
                                line += sequence[i].second[j][k].second[w];
                            } else if(aligned) {
                                line += '-';
                            }
                        }
                        // Main Nuc
                        if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x') {
                            line += sequence[i].second[j][k].first;
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                } else {
                    // If reverse strand, iterate backwards
                    for(size_t k = sequence[i].second[j].size()-1; k+1 > 0; k--) {
                        // Main nuc
                        if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x') {
                            line += getComplementCharacter(sequence[i].second[j][k].first);
                        } else if(aligned) {
                            line += '-';
                        }
                        // Gap nucs
                        for(size_t w = sequence[i].second[j][k].second.size()-1; w+1 > 0; w--) {
                            if(sequence[i].second[j][k].second[w] != '-') {
                                line += getComplementCharacter(sequence[i].second[j][k].second[w]);
                            } else if(aligned) {
                                line += '-';
                            }
                        }

                    }
                }
            } else if(aligned) {
                for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                        line += '-';
                    }
                    line += '-';
                }
            }
        }

        // Non-gap block - the only type being used currently
        if(blockExists[i].first) {
            // If forward strand
            if(blockStrand[i].first) {
                // Iterate through main nucs
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    // Gap nucs
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                        if(sequence[i].first[j].second[k] != '-') {
                            line += sequence[i].first[j].second[k];
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                    // Main nuc
                    if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                        line += sequence[i].first[j].first;
                    } else if(aligned) {
                        line += '-';
                    }
                }
            } else {
                // If reverse strand, iterate backwards
                for(size_t j = sequence[i].first.size()-1; j+1 > 0; j--) {
                    // Main nuc first since we are iterating in reverse direction
                    if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                        line += getComplementCharacter(sequence[i].first[j].first);
                    } else if(aligned) {
                        line += '-';
                    }

                    // Gap nucs
                    for(size_t k = sequence[i].first[j].second.size()-1; k+1 > 0; k--) {
                        if(sequence[i].first[j].second[k] != '-') {
                            line += getComplementCharacter(sequence[i].first[j].second[k]);
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                }   
            }
        } else if(aligned) {
            // If aligned sequence is required, print gaps instead if block does not exist
            for(size_t j = 0; j < sequence[i].first.size(); j++) {
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                    line+='-';
                }
                line+='-';
            }
        }

    }

    size_t ctr = 0;

    if(offset != 0) {
        for(size_t i = 0; i < line.length(); i++) {
            if(line[i] != '-') {
                if(ctr == (size_t)offset) {
                    // mark starting point
                    ctr = i;
                    break;
                }
                ctr++;
            }
        }
    }

    std::string currentLine = "";
    // From offset to end
    for(size_t i = ctr; i < line.length(); i++) {
        currentLine += line[i];
        if(currentLine.length() == lineSize) {
            fout << currentLine << '\n';
            currentLine = "";
        }
    }
    // From beginning to offset
    for(size_t i = 0; i < ctr; i++) {
        currentLine += line[i];
        if(currentLine.length() == lineSize) {
            fout << currentLine << '\n';
            currentLine = "";
        }
    }
    if(currentLine.length()) {
        fout << currentLine << '\n';
        currentLine = "";
    }

}

// Depth first traversal FASTA writer
void PangenomeMAT::Tree::printFASTAHelper(PangenomeMAT::Node* root, sequence_t& sequence,
    blockExists_t& blockExists, blockStrand_t& blockStrand, std::ofstream& fout, bool aligned) {

    // Apply mutations
    
    // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
    std::vector< std::tuple< int32_t, int32_t, bool, bool, bool, bool > > blockMutationInfo;

    // Block Mutations
    for(auto mutation: root->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        if(type == 1) {
            // insertion

            bool oldStrand;
            bool oldMut;
            if(secondaryBlockId != -1) {
                oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = true;

                // if insertion of inverted block takes place, the strand is backwards
                blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
            } else {
                oldStrand = blockStrand[primaryBlockId].first;
                oldMut = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = true;

                // if insertion of inverted block takes place, the strand is backwards
                blockStrand[primaryBlockId].first = !inversion;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion) );
        } else {
            bool oldMut;
            bool oldStrand;
            if(inversion) {
                // This means that this is not a deletion, but instead an inversion
                if(secondaryBlockId != -1) {
                    oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                    oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                    blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
                } else {
                    oldStrand = blockStrand[primaryBlockId].first;
                    oldMut = blockExists[primaryBlockId].first;
                    blockStrand[primaryBlockId].first = !oldStrand;
                }
                if(oldMut != true) {
                    std::cout << "Problem in PanMAT generation" << std::endl;
                }
                blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
            } else {
                // Actually a deletion

                if(secondaryBlockId != -1) {
                    oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                    oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                    blockExists[primaryBlockId].second[secondaryBlockId] = false;

                    // resetting strand to true during deletion
                    blockStrand[primaryBlockId].second[secondaryBlockId] = true;
                } else {
                    oldStrand = blockStrand[primaryBlockId].first;
                    oldMut = blockExists[primaryBlockId].first;
                    blockExists[primaryBlockId].first = false;

                    // resetting strand to true during deletion
                    blockStrand[primaryBlockId].first = true;
                }
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true) );
        }

    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++) {
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            // Either S, I or D

            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(type == PangenomeMAT::NucMutationType::NS) {
                // Substitution
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NI) {
                // Insertion
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::ND) {
                // Deletion
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }
                    }
                }
            }
        } 
        else {
            if(type == PangenomeMAT::NucMutationType::NSNPS) {
                // SNP Substitution
                newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1) {
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
            else if(type == PangenomeMAT::NucMutationType::NSNPI) {
                // SNP Insertion
                newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1) {
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
            else if(type == PangenomeMAT::NucMutationType::NSNPD) {
                // SNP Deletion
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    }
                } else {
                    if(nucGapPosition != -1) {
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
    
    sequence_t sequencePrint = sequence;
    blockExists_t blockExistsPrint = blockExists;
    blockStrand_t blockStrandPrint = blockStrand;

    if(root->children.size() == 0) {
        // Print sequence

        fout << '>' << root->identifier << std::endl;
        int offset = 0;
        if(!aligned && circularSequences.find(root->identifier) != circularSequences.end()) {
            // If MSA is to be printed, offset doesn't matter
            offset = circularSequences[root->identifier];
        }
        sequence_t sequencePrint = sequence;
        blockExists_t blockExistsPrint = blockExists;
        blockStrand_t blockStrandPrint = blockStrand;

        if(rotationIndexes.find(root->identifier) != rotationIndexes.end() && rotationIndexes[root->identifier] != 0) {
            int ctr = -1, rotInd = 0;
            for(size_t i = 0; i < blockExistsPrint.size(); i++) {
                if(blockExistsPrint[i].first) {
                    ctr++;
                }
                if(ctr == rotationIndexes[root->identifier]) {
                    rotInd = i;
                    break;
                }
            }
            rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
            rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
            rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
        }

        if(sequenceInverted.find(root->identifier) != sequenceInverted.end() && sequenceInverted[root->identifier]) {
            reverse(sequencePrint.begin(), sequencePrint.end());
            reverse(blockExistsPrint.begin(), blockExistsPrint.end());
            reverse(blockStrandPrint.begin(), blockStrandPrint.end());
        }

        PangenomeMAT::printSequenceLines(sequencePrint, blockExistsPrint, blockStrandPrint, 70, aligned, fout, offset);
    } else {
        // DFS on children
        for(PangenomeMAT::Node* child: root->children) {
            printFASTAHelper(child, sequence, blockExists, blockStrand, fout, aligned);
        }
    }


    // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<3>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
        }
    }

    // Undo nuc mutations when current node and its subtree have been processed
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }
}

void PangenomeMAT::Tree::printFASTA(std::ofstream& fout, bool aligned) {
    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(blocks.size() + 1, {false, {}});
    blockStrand_t blockStrand(blocks.size() + 1, {true, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        blockStrand[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], true);
    }

    int32_t maxBlockId = 0;

    // Create consensus sequence of blocks
    for(size_t i = 0; i < blocks.size(); i++) {

        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT::getNucleotideFromCode(nucCode);
                
                if(secondaryBlockId != -1) {
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }

            if(endFlag) {
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1) {
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    // Assigning nucleotide gaps in blocks
    for(size_t i = 0; i < gaps.size(); i++) {
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            if(secondaryBId != -1) {
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }

    // Run depth first traversal to extract sequences
    printFASTAHelper(root, sequence, blockExists, blockStrand, fout, aligned);

}

void PangenomeMAT::Tree::generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout) {
    for(auto u: allNodes) {
        if(u.second->children.size() == 0) {
            // sequence
            std::string sequenceId = u.first;
            fout << ">" << sequenceId << "\n";

            std::string line;
            std::map< int, std::string > positionToSequence;

            while(getline(fin, line, '\n')) {
                if(line.length() > 2 && line.substr(0,2) == "s\t") {
                    std::vector< std::string > words;
                    stringSplit(line, '\t', words);
                    if(words.size() != 7) {
                        std::cout << "Line not in correct format. Line size: " << words.size() << std::endl;
                        return;
                    }

                    int startPosition = std::stoll(words[2]);
                    if(words[1] != sequenceId) {
                        continue;
                    }

                    bool strand = (words[4] == "+"?true: false);
                    std::string sequence = words[words.size()-1];
                    std::string strippedSequence;
                    for(auto u: sequence) {
                        if(u != '-') {
                            if(strand) {
                                strippedSequence+=u;
                            } else {
                                strippedSequence+=getComplementCharacter(u);
                            }
                        }
                    }
                    if(!strand) {
                        std::reverse(strippedSequence.begin(), strippedSequence.end());
                    }
                    positionToSequence[startPosition] = strippedSequence;
                }
            }
            int nextExpected = 0;
            int endLength = 0;
            std::string fullSequence;
            for(auto u: positionToSequence) {
                if(nextExpected == 0 && u.first != nextExpected) {
                    nextExpected = u.first;
                    endLength = u.first;
                }
                if(u.first != nextExpected) {
                    std::cout << "Error in positions" << std::endl;
                }
                fullSequence+=u.second;
                nextExpected+=u.second.length();
            }
            fullSequence = fullSequence.substr(fullSequence.length() - endLength) + fullSequence.substr(0,fullSequence.length() - endLength);
            fout << fullSequence << '\n';
            fin.clear();
            fin.seekg(0);
        }
    }
}

void PangenomeMAT::Tree::printMAF(std::ofstream& fout) {
    std::vector< std::string > sequenceNames;

    std::map< std::vector< uint32_t >, std::vector< std::pair< int,int > > > blocksWithSameSequences;

    for(auto b: blocks) {
        blocksWithSameSequences[b.consensusSeq].push_back(std::make_pair(b.primaryBlockId, b.secondaryBlockId));
    }

    // sequence name, primary bid, secondary bid -> actual sequence offset
    std::map< std::tuple< std::string, int32_t, int32_t >, int32_t > sequenceBlockToStartPoint;
    std::map< std::string, int32_t > sequenceLengths;

    for(auto u: allNodes) {
        if(u.second->children.size() == 0) {
            sequenceNames.push_back(u.first);
            sequence_t sequence;
            blockExists_t blockExists;
            blockStrand_t blockStrand;
            int rotIndex;
            getSequenceFromReference(sequence, blockExists, blockStrand, u.first, &rotIndex);

            // List of block IDs to account for rotation and inversion
            std::vector< int > blockIds(sequence.size());
            for(size_t i = 0; i < sequence.size(); i++) {
                blockIds[i] = i;
            }
            if(rotationIndexes.find(u.first) != rotationIndexes.end()) {
                rotate(blockIds.begin(), blockIds.begin() + rotIndex, blockIds.end());
            }
            if(sequenceInverted.find(u.first) != sequenceInverted.end()
                && sequenceInverted[u.first]) {
                reverse(blockIds.begin(), blockIds.end());
            }

            int ctr = 0;

            for(size_t i = 0; i < sequence.size(); i++) {
                // for(size_t j = 0; j < sequence[i].second.size(); j++) {
                //     if(blockExists[i].second[j]) {
                //         sequenceBlockToStartPoint[std::make_tuple(u.first, i, j)] = ctr;
                //         for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                //             for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                //                 if(sequence[i].second[j][k].second[w] != '-' && sequence[i].second[j][k].second[w] != 'x') {
                //                     ctr++;
                //                 }
                //             }
                //             if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x') {
                //                 ctr++;
                //             }
                //         }
                //     }
                // }
                if(blockExists[i].first) {
                    sequenceBlockToStartPoint[std::make_tuple(u.first, blockIds[i], -1)] = ctr;
                    for(size_t j = 0; j < sequence[i].first.size(); j++) {
                        for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                            if(sequence[i].first[j].second[k] != '-' && sequence[i].first[j].second[k] != 'x') {
                                ctr++;
                            }
                        }
                        if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                            ctr++;
                        }
                    }
                }
            }
            sequenceLengths[u.first] = ctr;
            if(circularSequences.find(u.first) != circularSequences.end()) {
                int offset = circularSequences[u.first];
                for(auto& v: sequenceBlockToStartPoint) {
                    if(std::get<0>(v.first) != u.first) {
                        continue;
                    }
                    v.second -= offset;
                    if(v.second < 0) {
                        v.second += ctr;
                    }
                }
            }
        }
    }


    fout << "##maf version=1\n";

    for(auto& common: blocksWithSameSequences) {
        tbb::concurrent_unordered_map< std::string, std::pair< std::pair< int, int >, std::pair< std::string, bool > > > sequenceIdToSequence;
        fout << "a\n";
        for(auto& b: common.second) {
            int primaryBlockId = b.first;
            int secondaryBlockId = b.second;

            tbb::parallel_for_each(sequenceNames, [&](auto u) {
                block_t sequence;
                bool blockExists = false, blockStrand = true;
                getBlockSequenceFromReference(sequence, blockExists, blockStrand, u, primaryBlockId,
                    secondaryBlockId);
                if(!blockExists) {
                    return;
                }

                std::string stringSequence;
                for(size_t i = 0; i < sequence.size(); i++) {
                    for(size_t j = 0; j < sequence[i].second.size(); j++) {
                        if(sequence[i].second[j] != '-' && sequence[i].second[j] != 'x') {
                            stringSequence += sequence[i].second[j];
                        } else {
                            stringSequence += '-';
                        }
                    }
                    if(sequence[i].first != '-' && sequence[i].first != 'x') {
                        stringSequence += sequence[i].first;
                    } else {
                        stringSequence += '-';
                    }
                }
                sequenceIdToSequence[u] = std::make_pair(std::make_pair(primaryBlockId, secondaryBlockId), std::make_pair(stringSequence, blockStrand));
            });

            for(auto u: sequenceIdToSequence) {
                fout << "s\t" << u.first << "\t"<< sequenceBlockToStartPoint[std::make_tuple(u.first, u.second.first.first, u.second.first.second)] <<"\t" << u.second.second.first.length() << "\t" << (u.second.second.second? "+\t":"-\t") << sequenceLengths[u.first] << "\t" << u.second.second.first << "\n";
            }
        }
        fout << "\n";
    }
}

void PangenomeMAT::Tree::printVCFParallel(std::string reference, std::ofstream& fout) {

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "Error: Reference sequence with matching name not found!") {
        std::cerr << referenceSequence << std::endl;
        return;
    }

    size_t recordID = 0;

    std::mutex vcfMapMutex;
    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    tbb::parallel_for_each(allNodes, [&](auto& n) {
        if(n.second->children.size() == 0 && n.first != reference) {
            std::string altSequence = getStringFromReference(n.first);
            if(altSequence.length() != referenceSequence.length()) {
                std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
                return;
            }

            std::string currentRefString, currentAltString;
            int currentCoordinate = 0;

            int diffStart = 0;

            for(size_t i = 0; i < referenceSequence.length(); i++) {

                if(referenceSequence[i] == '-' && altSequence[i] == '-') {
                    continue;
                } else if(referenceSequence[i] != '-' && altSequence[i] == '-') {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }

                    currentRefString += referenceSequence[i];
                } else if(referenceSequence[i] == '-' && altSequence[i] != '-') {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }

                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] != altSequence[i]) {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }
                    if(currentRefString == currentAltString) {
                        currentRefString = "";
                        currentAltString = "";
                        diffStart = currentCoordinate;
                    }
                    currentRefString += referenceSequence[i];
                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] == altSequence[i]) {
                    if(currentRefString == currentAltString) {
                        // Reset
                        diffStart = currentCoordinate;
                        currentRefString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString = currentRefString;
                    } else {
                        // Create VCF record at position i
                        if(currentRefString == "") {
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

                if(referenceSequence[i] != '-') {
                    currentCoordinate++;
                }
            }

            if(currentRefString != currentAltString) {
                vcfMapMutex.lock();
                vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                vcfMapMutex.unlock();

                // Reset
                diffStart = referenceSequence.size();
                currentRefString = "";
                currentAltString = currentRefString;
            }
        }
    }

    if(root->children.size() == 0){
        // Print sequence
        fout << '>' << root->identifier << std::endl;

        printSequenceLines(sequence, blockExists, 80, aligned, fout);

    } else {
        // DFS on children
        for(PangenomeMAT::Node* child: root->children){
            printFASTAHelper(child, sequence, blockExists, fout, aligned);
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

void PangenomeMAT::Tree::dfsExpansion(PangenomeMAT::Node* node, std::vector< PangenomeMAT::Node* >& vec) {
    vec.push_back(node);
    for(auto child: node->children) {
        dfsExpansion(child, vec);
    }
}

std::string PangenomeMAT::Tree::getNewickString(Node* node) {

    std::vector< PangenomeMAT::Node* > traversal;
    dfsExpansion(node, traversal);

    std::string newick;

    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<int32_t> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        int32_t branch_length = (n->nucMutation).size();

        if(curr_level < level) {
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
                    newick += std::to_string(branch_length);
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
                    newick += std::to_string(branch_length_stack.top());
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->children.size() == 0) {
                
                newick += ',';
                newick += n->identifier;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
            }
            if(endFlag){
                break;
            }
        }

        // End character to incorporate for gaps at the end
        sequence[blocks[i].blockId].push_back({'x',{}});
    }

    // Assigning gaps
    for(size_t i = 0; i < gaps.position.size(); i++){
        int bId = gaps.blockId[i];

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
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
            newick += std::to_string(branch_length_stack.top());
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

}

// Merge parent node and child node into parent node
void mergeNodes(PangenomeMAT::Node* par, PangenomeMAT::Node* chi){
    
    par->identifier = chi->identifier;
    par->branchLength += chi->branchLength;
    par->children = chi->children;

    // For block mutations, we cancel out irrelevant mutations
    std::unordered_map< int, PangenomeMAT::BlockMutationType > bidMutations;

    for(auto mutation: par->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & 0xFFFFFF);
        int type = (mutation & 0x1);
        if(type == PangenomeMAT::BlockMutationType::BI){
            bidMutations[bid] = PangenomeMAT::BlockMutationType::BI;
        } else {
            if(bidMutations.find(bid) != bidMutations.end()){
                if(bidMutations[bid] == PangenomeMAT::BlockMutationType::BI){
                    // If it was insertion earlier, cancel out
                    bidMutations.erase(bid);
                }
                // Otherwise, it remains deletion
            } else {
                bidMutations[bid] = PangenomeMAT::BlockMutationType::BD;
            }
        }
    }

    for(auto mutation: chi->blockMutation.condensedBlockMut){
        int bid = ((mutation >> 8) & 0xFFFFFF);
        int type = (mutation & 0x1);
        if(type == PangenomeMAT::BlockMutationType::BI){
            bidMutations[bid] = PangenomeMAT::BlockMutationType::BI;
        } else {
            if(bidMutations.find(bid) != bidMutations.end()){
                if(bidMutations[bid] == PangenomeMAT::BlockMutationType::BI){
                    // If it was insertion earlier, cancel out
                    bidMutations.erase(bid);
                }
                // Otherwise, it remains deletion
            } else {
                bidMutations[bid] = PangenomeMAT::BlockMutationType::BD;
            }
        }
    }

    PangenomeMAT::BlockMut newBlockMutation;
    for(auto mutation: bidMutations){
        if(mutation.second == PangenomeMAT::BlockMutationType::BI){
            newBlockMutation.condensedBlockMut.push_back((( mutation.first << 8 ) ^ 0x2));
        } else {
            newBlockMutation.condensedBlockMut.push_back((( mutation.first << 8 ) ^ 0x1));
        }
    }

    par->blockMutation = newBlockMutation;

    for(auto mutation: chi->nucMutation){
        par->nucMutation.push_back(mutation);
    }

    delete chi;
}

// Replace old type, char pair with new type char pair
std::pair< int, int > replaceMutation(std::pair<int,int> oldMutation, std::pair<int, int> newMutation){
    std::pair<int, int> ans = newMutation;
    if(oldMutation.first == newMutation.first){
        ans = newMutation;
    } else if(oldMutation.first == PangenomeMAT::NucMutationType::NSNPS){
        // Insertion after substitution (doesn't make sense but just in case)
        if(newMutation.first == PangenomeMAT::NucMutationType::NSNPI){
            ans.first = PangenomeMAT::NucMutationType::NSNPS;
        } else if(newMutation.first == PangenomeMAT::NucMutationType::NSNPD){
            ans = newMutation;
        }
    } else if(oldMutation.first == PangenomeMAT::NucMutationType::NSNPI){
        if(newMutation.first == PangenomeMAT::NucMutationType::NSNPS){
            ans.first = PangenomeMAT::NucMutationType::NSNPI;
        } else if(newMutation.first == PangenomeMAT::NucMutationType::NSNPD){
            // Cancel out the two mutations if deletion after insertion
            ans = std::make_pair(404, 404);
        }
    } else if(oldMutation.first == PangenomeMAT::NucMutationType::NSNPD){
        if(newMutation.first == PangenomeMAT::NucMutationType::NSNPI){
            ans.first = PangenomeMAT::NucMutationType::NSNPS;
        } else if(newMutation.first == PangenomeMAT::NucMutationType::NSNPS){
            // Substitution after deletion. Doesn't make sense but still
            ans.first = PangenomeMAT::NucMutationType::NSNPI;
        }
    }
    return ans;
}

bool debugSimilarity(const std::vector< PangenomeMAT::NucMut > array1, const std::vector< PangenomeMAT::NucMut > array2){
    std::map< std::tuple< int, int, int >, std::pair< int, int > > mutationRecords1, mutationRecords2;

    for(auto mutation: array1){
        int bid = ((mutation.condensed) >> 8);
        int pos = mutation.position;
        int gapPos = mutation.gapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.condensed) & 0x7);
        int len = (((mutation.condensed) >> 3) & 0x1F);

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT::NucMutationType::NS:
                newType = PangenomeMAT::NucMutationType::NSNPS;
                break;
            case PangenomeMAT::NucMutationType::ND:
                newType = PangenomeMAT::NucMutationType::NSNPD;
                break;
            case PangenomeMAT::NucMutationType::NI:
                newType = PangenomeMAT::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar;
            if(type < 3){
                newChar = (((mutation.nucs) >> (4*(15 - i))) & 0xF);
            } else {
                // SNP
                newChar = ((mutation.condensed >> 3) & 0xF);
            }

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords1.find(std::make_tuple( bid, pos, gapPos + i )) == mutationRecords1.end()){
                    mutationRecords1[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( bid, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords1[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( bid, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords1.find(std::make_tuple( bid, pos + i, gapPos )) == mutationRecords1.end()){
                    mutationRecords1[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( bid, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords1[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( bid, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    for(auto mutation: array2){
        int bid = ((mutation.condensed) >> 8);
        int pos = mutation.position;
        int gapPos = mutation.gapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.condensed) & 0x7);
        int len = (((mutation.condensed) >> 3) & 0x1F);

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT::NucMutationType::NS:
                newType = PangenomeMAT::NucMutationType::NSNPS;
                break;
            case PangenomeMAT::NucMutationType::ND:
                newType = PangenomeMAT::NucMutationType::NSNPD;
                break;
            case PangenomeMAT::NucMutationType::NI:
                newType = PangenomeMAT::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar;
            if(type < 3){
                newChar = (((mutation.nucs) >> (4*(15 - i))) & 0xF);
            } else {
                // SNP
                newChar = ((mutation.condensed >> 3) & 0xF);
            }

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords2.find(std::make_tuple( bid, pos, gapPos + i )) == mutationRecords2.end()){
                    mutationRecords2[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( bid, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords2[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( bid, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords2.find(std::make_tuple( bid, pos + i, gapPos )) == mutationRecords2.end()){
                    mutationRecords2[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( bid, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords2[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( bid, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    std::vector< std::tuple< int, int, int, int, int > > mutationArray1, mutationArray2;
    for(auto u: mutationRecords1){
        mutationArray1.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), u.second.first, u.second.second ) );
    }
    for(auto u: mutationRecords2){
        mutationArray2.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), u.second.first, u.second.second ) );
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

std::vector< PangenomeMAT::NucMut > consolidateNucMutations(const std::vector< PangenomeMAT::NucMut >& nucMutation){
    // bid, pos, gap_pos -> type, nuc
    std::map< std::tuple< int, int, int >, std::pair< int, int > > mutationRecords;
    for(auto mutation: nucMutation){
        int bid = ((mutation.condensed) >> 8);
        int pos = mutation.position;
        int gapPos = mutation.gapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.condensed) & 0x7);
        int len = (((mutation.condensed) >> 3) & 0x1F);

        if(type >= 3){
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type){
            case PangenomeMAT::NucMutationType::NS:
                newType = PangenomeMAT::NucMutationType::NSNPS;
                break;
            case PangenomeMAT::NucMutationType::ND:
                newType = PangenomeMAT::NucMutationType::NSNPD;
                break;
            case PangenomeMAT::NucMutationType::NI:
                newType = PangenomeMAT::NucMutationType::NSNPI;
                break;
        }

        for(int i = 0; i < len; i++){
            int newChar;
            if(type < 3){
                newChar = (((mutation.nucs) >> (4*(15 - i))) & 0xF);
            } else {
                // SNP
                newChar = ((mutation.condensed >> 3) & 0xF);
            }

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1){
                if(mutationRecords.find(std::make_tuple( bid, pos, gapPos + i )) == mutationRecords.end()){
                    mutationRecords[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( bid, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords[std::make_tuple( bid, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords.erase(std::make_tuple( bid, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords.find(std::make_tuple( bid, pos + i, gapPos )) == mutationRecords.end()){
                    mutationRecords[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( bid, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404){
                        mutationRecords[std::make_tuple( bid, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords.erase(std::make_tuple( bid, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    // bid, pos, gapPos, type, char
    std::vector< std::tuple< int, int, int, int, int > > mutationArray;
    for(auto u: mutationRecords){
        mutationArray.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), u.second.first, u.second.second ) );
    }
    
    // mutation array is already sorted since mutationRecord was sorted
    std::vector< PangenomeMAT::NucMut > consolidatedMutationArray;

    for(size_t i = 0; i < mutationArray.size(); i++){
        size_t j = i + 1;
        for(; j < std::min(i + 16, mutationArray.size()); j++){
            if(std::get<2>(mutationArray[i]) != -1){
                // gapPos exists
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && std::get<1>(mutationArray[i]) == std::get<1>(mutationArray[j])
                    && std::get<3>(mutationArray[i]) == std::get<3>(mutationArray[j]) && (size_t)(std::get<2>(mutationArray[j]) - std::get<2>(mutationArray[i])) == j - i)){
                    break;
                }
            } else {
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && (size_t)(std::get<1>(mutationArray[j]) - std::get<1>(mutationArray[i])) == j - i
                    && std::get<3>(mutationArray[i]) == std::get<3>(mutationArray[j]) && std::get<2>(mutationArray[j]) == std::get<2>(mutationArray[i]))){
                    break;
                }
            }
        }

        if(j - i <= 1){
            consolidatedMutationArray.push_back(PangenomeMAT::NucMut(mutationArray[i]));

            continue;
        }
        // combine mutations from i to j
        auto newMutation = PangenomeMAT::NucMut(mutationArray, i, j);

        consolidatedMutationArray.push_back(newMutation);

        i = j - 1;
    }

    return consolidatedMutationArray;

}

void dfsExpansion(PangenomeMAT::Node* node, std::vector< PangenomeMAT::Node* >& vec){
    vec.push_back(node);
    for(auto child: node->children){
        dfsExpansion(child, vec);
    }
}

std::string PangenomeMAT::Tree::getNewickString(Node* node){
    std::vector< PangenomeMAT::Node* > traversal;
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

    return newick;

}

void compressTreeParallel(PangenomeMAT::Node* node, size_t level){
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

            debugSimilarity(oldVector, node->children[i]->nucMutation);

            compressTreeParallel(node->children[i], level + 1);
        }
    });
}

PangenomeMAT::Node* subtreeExtractParallelHelper(PangenomeMAT::Node* node, const tbb::concurrent_unordered_map< PangenomeMAT::Node*, size_t >& ticks){
    if(ticks.find(node) == ticks.end()){
        return nullptr;
    }

    PangenomeMAT::Node* newNode = new PangenomeMAT::Node(node->identifier, node->branchLength);

    for(auto mutation: node->nucMutation){
        newNode->nucMutation.push_back(mutation);
    }

    for(auto mutation: node->blockMutation.condensedBlockMut){
        newNode->blockMutation.condensedBlockMut.push_back(mutation);
    }

    newNode->children.resize(node->children.size(), nullptr);

    tbb::parallel_for(tbb::blocked_range(0, (int)node->children.size()), [&](tbb::blocked_range<int> r){
        for(int i = r.begin(); i < r.end(); i++){
            PangenomeMAT::Node* child = node->children[i];
            if(ticks.find(child) != ticks.end()){

                PangenomeMAT::Node* newChild = subtreeExtractParallelHelper(child, ticks);

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

PangenomeMAT::Node* PangenomeMAT::Tree::subtreeExtractParallel(std::vector< std::string > nodeIds){
    tbb::concurrent_vector< PangenomeMAT::Node* > requiredNodes;
    
    tbb::parallel_for_each(nodeIds.begin(), nodeIds.end(), [&]( std::string& id ) {
        requiredNodes.push_back(allNodes[id]);
    });

    tbb::concurrent_unordered_map< PangenomeMAT::Node*, size_t > ticks;

    tbb::parallel_for_each(requiredNodes.begin(), requiredNodes.end(), [&](PangenomeMAT::Node*& node){
        Node* current = node;

        while(current != nullptr){
            ticks[current]++;
            current = current->parent;
        }
    });

    PangenomeMAT::Node* newTreeRoot = subtreeExtractParallelHelper(root, ticks);

    compressTreeParallel(newTreeRoot, 1);

    return newTreeRoot;
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

void PangenomeMAT::Tree::writeToFile(std::ofstream& fout, Node* node){
    if(node == nullptr){
        node = root;
    }

    MAT::tree treeToWrite;
    getNodesPreorder(node, treeToWrite);

    std::string newick = getNewickString(node);

    treeToWrite.set_newick(newick);

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

std::string PangenomeMAT::Tree::getStringFromReference(std::string reference){
    Node* referenceNode = nullptr;
    
    for(auto u: allNodes){
        if(u.second->children.size() == 0 && u.first == reference){
            referenceNode = u.second;
            break;
        }
    }

    if(referenceNode == nullptr){
        return "No such leaf node found.";
    }

    // Categorize gaps by blockId
    std::map< int, std::vector< std::pair< int, int > > > gapSplit;
    
    for(size_t i = 0; i < gaps.position.size(); i++){

        int bId = gaps.blockId[i];

        int len = gaps.gapLength[i];

        int pos = gaps.position[i];
        gapSplit[bId].push_back( std::make_pair(pos, len) );

    }

    std::vector< PangenomeMAT::Node* > path;
    Node* it = referenceNode;

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
            int type = (((*node)->nucMutation[i].condensed) & 0x7);
            char newVal = '-';

            if(type < 3){
                // Either S, I or D

                int len = ((((*node)->nucMutation[i].condensed) >> 3) & 0x1F);

                if(type == PangenomeMAT::NucMutationType::NS){
                    for(int j = 0; j < len; j++){
                        newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(15-j))) & 15);

                        sequence[bid][pos + j].first = newVal;
                    }
                }
                else if(type == PangenomeMAT::NucMutationType::NI){

                    if(gapPos == -1){
                        for(int j = 0; j < len; j++){
                            newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(15-j))) & 15);
                            
                            sequence[bid][pos + j].first = newVal;
                        }
                    } else {
                        for(int j = 0; j < len; j++){
                            newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(15-j))) & 15);

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
                if(type == PangenomeMAT::NucMutationType::NSNPS){
                    newVal = getNucleotideFromCode((((*node)->nucMutation[i].condensed) >> 3) & 0xF);
                    sequence[bid][pos].first = newVal;
                }
                else if(type == PangenomeMAT::NucMutationType::NSNPI){
                    newVal = getNucleotideFromCode((((*node)->nucMutation[i].condensed) >> 3) & 0xF);

                    if(gapPos == -1){
                        sequence[bid][pos].first = newVal;
                    } else {
                        sequence[bid][pos].second[gapPos] = newVal;
                    }
                }
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

    std::string sequenceString;
    for(size_t i = 0; i < blocks.size(); i++){
        if(sequence.find(blocks[i].blockId) == sequence.end()){
            for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
                for(int k = 0; k < 8; k++){
                    sequenceString += '-';
                }
            }
            sequenceString += '-'; // For last character that might have been inserted (x)

            for(auto g: gapSplit[blocks[i].blockId]){
                for(int j = 0; j < g.second; j++){
                    sequenceString += '-';
                }
            }
        } else {
            for(size_t j = 0; j < sequence[blocks[i].blockId].size(); j++){
                for(size_t k = 0; k < sequence[blocks[i].blockId][j].second.size(); k++){
                    sequenceString += sequence[blocks[i].blockId][j].second[k];
                }

                if(sequence[blocks[i].blockId][j].first != 'x'){
                    sequenceString += sequence[blocks[i].blockId][j].first;
                } else {
                    sequenceString += '-';
                }
            }
        }
    }

    return sequenceString;

}

std::string stripString(std::string s){
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

std::string stripGaps(std::string sequenceString){
    std::string result;
    for(auto u: sequenceString){
        if(u != '-'){
            result+=u;
        }
    }
    return result;
}

std::string PangenomeMAT::Tree::getSequenceFromVCF(std::string sequenceId, std::ifstream& fin){
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

    // column headers
    std::getline(fin, line);

    std::vector< std::string > columnWords;
    std::string word;

    for(size_t i = 0; i < line.size(); i++){
        if(line[i] != ' '){
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
        std::cout << "sequence not found!" << std::endl;
        return "";
    }

    std::vector< std::pair< char, std::vector< char > > > alteredSequence;
    for(auto u: referenceSequence){
        alteredSequence.push_back({u, {}});
    }

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

        // if(position >= alteredSequence.size()){
        //     std::cout << "position too high " << position << " " << alteredSequence.size() << std::endl;
        // }

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

        // if(choice >= altChoices.size()){
        //     std::cout << altChoices.size() << " " << choice << " " << altStrings << " " << position << std::endl;
        // }

        std::string alt = altChoices[choice];

        if(ref != "."){
            int len = ref.length();
            for(int i = position; i < position + len; i++){
                // if(i >= alteredSequence.size()){
                //     std::cout << "index exceeding!!!" << std::endl;
                // }
                alteredSequence[i].first = '-';
            }
        }

        if(alt != "."){
            if(alt.length() && alteredSequence[position].second.size()){
                std::cout << "alternate sequence already exists!" << std::endl;
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

    std::cout << (alteredSequenceOriginal.length() == finalSequence.length()) << (alteredSequenceOriginal == finalSequence) << std::endl;

    return finalSequence;

}

void PangenomeMAT::Tree::printVCF(std::string reference, std::ofstream& fout){

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "No such leaf node found."){
        std::cerr << referenceSequence << std::endl;
        return;
    }

    size_t recordID = 0;

    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    for(auto n: allNodes){
        if(n.second->children.size() == 0 && n.first != reference){
            std::string altSequence = getStringFromReference(n.first);
            if(altSequence.length() != referenceSequence.length()){
                std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
                return;
            }

            auto blockExists = blockExistsGlobal;

            std::vector< PangenomeMAT::Node* > path;

            Node* it = u.second;
            while(it != root) {
                path.push_back(it);
                it = it->parent;
            }
            path.push_back(root);
            std::reverse(path.begin(), path.end());
            for(auto node: path) {
                for(auto mutation: node->blockMutation) {
                    int primaryBlockId = mutation.primaryBlockId;
                    int secondaryBlockId = mutation.secondaryBlockId;
                    int type = (mutation.blockMutInfo);
                    
                    if(type == PangenomeMAT::BlockMutationType::BI) {
                        if(secondaryBlockId != -1) {
                            blockExists[primaryBlockId].second[secondaryBlockId] = true;
                        } else {
                            blockExists[primaryBlockId].first = true;
                        }
                    } else {
                        if(secondaryBlockId != -1) {
                            blockExists[primaryBlockId].second[secondaryBlockId] = false;
                        } else {
                            blockExists[primaryBlockId].first = false;
                        }
                    }
                }
            }
            std::vector< std::pair< int32_t, int32_t > > currentPath;
            for(size_t i = 0; i < blockExists.size(); i++) {
                if(blockExists[i].first) {
                    currentPath.push_back(std::make_pair(i, -1));
                }
                for(size_t j = 0; j < blockExists[i].second.size(); j++) {
                    if(blockExists[i].second[j]) {
                        currentPath.push_back(std::make_pair(i, j));
                    }
                }
            }
            paths[u.second->identifier] = currentPath;
            for(size_t i = 1; i < currentPath.size(); i++) {
                edges.insert(std::make_pair(currentPath[i-1], currentPath[i]));
            }
        });
        std::map< std::pair< int32_t, int32_t >, uint64_t > nodeIds;
        uint64_t ctr = 0;
        for(auto u: nodes) {
            nodeIds[u.first] = ctr;
            ctr++;
        }
        for(auto u: nodes) {
            fout << "S\t" << nodeIds[u.first] << "\t" << u.second << "\n";
        }
        for(auto u: edges) {
            fout << "L\t" << nodeIds[u.first] << "\t+\t" << nodeIds[u.second] << "\t+\t0M\n";
        }
        for(auto u: paths) {
            fout << "P\t" << u.first << "\t";
            for(size_t i = 0; i < u.second.size(); i++) {
                fout << nodeIds[u.second[i]] << "+";
                if(i != u.second.size() - 1) {
                    fout << ",";
                }
            }
            fout << "\t*\n";
        }
    } else {
        size_t autoIncrId = 0;
        std::map< std::pair< std::tuple< size_t, size_t, size_t >, std::string >,
            std::pair< size_t, bool > > allSequenceNodes;
        std::mutex allSequenceNodeMutex;
        tbb::concurrent_unordered_map< std::string, std::vector< size_t > > paths;
        tbb::concurrent_unordered_map< std::string, std::vector< bool > > strandPaths;

        for(const auto& u: allNodes){
        // tbb::parallel_for_each(allNodes, [&](const auto& u) {
            if(u.second->children.size() != 0) {
                // return;
                continue;
            }

            sequence_t sequence;
            blockExists_t blockExists;
            blockStrand_t blockStrand;
            getSequenceFromReference(sequence, blockExists, blockStrand, u.first);

            // if(u.first == "NZ_CP006027.1") {
            //     for(int i = 0; i < sequence.size(); i++) {
            //         if(blockExists[i].first) {
            //             for(int j = sequence[i].first.size()-1; j>=0; j--) {
            //                 if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
            //                     std::cout << sequence[i].first[j].first << std::endl;
            //                 }
            //                 for(int k = sequence[i].first[j].second.size()-1; k>=0; k--) {
            //                     if(sequence[i].first[j].second[k] != '-' && sequence[i].first[j].second[k] != 'x') {
            //                         std::cout << sequence[i].first[j].second[k] << std::endl;
            //                     }
            //                 }
            //             }
            //             break;
            //         }
            //     }
            // }

            std::string currentSequence;
            std::vector< size_t > sequenceNodeIds;
            std::vector< bool > sequenceStrands;

            for(size_t i = 0; i < sequence.size(); i++) {
                if(blockExists[i].first) {
                    if(blockStrand[i].first) {
                        std::tuple< size_t, size_t, size_t > currentStart;
                        for(size_t j = 0; j < sequence[i].first.size(); j++) {
                            for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                                if(currentSequence.length() == 0) {
                                    currentStart = std::make_tuple(i,j,k);
                                }
                                currentSequence += sequence[i].first[j].second[k];
                                if(currentSequence.length() == 32) {
                                    currentSequence = stripGaps(currentSequence);
                                    if(currentSequence.length()) {
                                        allSequenceNodeMutex.lock();
                                        if(allSequenceNodes.find(std::make_pair(currentStart,
                                            currentSequence)) == allSequenceNodes.end()) {
                                            allSequenceNodes[std::make_pair(currentStart,
                                                currentSequence)] = std::make_pair(autoIncrId, true);
                                            sequenceNodeIds.push_back(autoIncrId);
                                            sequenceStrands.push_back(true);
                                            autoIncrId++;
                                        } else {
                                            sequenceNodeIds.push_back(
                                                allSequenceNodes[std::make_pair(currentStart,
                                                currentSequence)].first);
                                            sequenceStrands.push_back(true);
                                        }
                                        allSequenceNodeMutex.unlock();
                                    }
                                    currentSequence = "";
                                }
                            }
                            if(currentSequence.length() == 0) {
                                currentStart = std::make_tuple(i,j,-1);
                            }
                            currentSequence += sequence[i].first[j].first;
                            if(currentSequence.length() == 32) {
                                currentSequence = stripGaps(currentSequence);
                                if(currentSequence.length()) {
                                    allSequenceNodeMutex.lock();
                                    if(allSequenceNodes.find(std::make_pair(currentStart,
                                        currentSequence)) == allSequenceNodes.end())  {
                                        allSequenceNodes[std::make_pair(currentStart,
                                            currentSequence)] = std::make_pair(autoIncrId, true);
                                        sequenceNodeIds.push_back(autoIncrId);
                                        sequenceStrands.push_back(true);
                                        autoIncrId++;
                                    } else {
                                        sequenceNodeIds.push_back(
                                            allSequenceNodes[std::make_pair(currentStart,
                                            currentSequence)].first);
                                        sequenceStrands.push_back(true);
                                    }
                                    allSequenceNodeMutex.unlock();
                                }
                                currentSequence = "";
                            }
                        }
                        if(currentSequence.length()) {
                            currentSequence = stripGaps(currentSequence);
                            if(currentSequence.length()) {
                                allSequenceNodeMutex.lock();
                                if(allSequenceNodes.find(std::make_pair(currentStart,
                                    currentSequence)) == allSequenceNodes.end()) {
                                    allSequenceNodes[std::make_pair(currentStart, currentSequence)]
                                        = std::make_pair(autoIncrId, true);
                                    sequenceNodeIds.push_back(autoIncrId);
                                    sequenceStrands.push_back(true);
                                    autoIncrId++;
                                } else {
                                    sequenceNodeIds.push_back(
                                        allSequenceNodes[std::make_pair(currentStart,
                                        currentSequence)].first);
                                    sequenceStrands.push_back(true);
                                }
                                allSequenceNodeMutex.unlock();
                            }
                            currentSequence = "";
                        }
                    } else {
                        std::tuple< size_t, size_t, size_t > currentStart;
                        for(size_t j = sequence[i].first.size()-1; j + 1 > 0; j--) {
                            currentSequence += sequence[i].first[j].first;
                            currentStart = std::make_tuple(i,j,-1);
                            if(currentSequence.length() == 32) {
                                currentSequence = stripGaps(currentSequence);
                                if(currentSequence.length()) {
                                    // Since the GFA stores the strand parameter, the reverse
                                    // complement will be computed anyway
                                    std::reverse(currentSequence.begin(), currentSequence.end());
                                    // if(u.first == "NZ_CP006027.1") {
                                    //     std::cout << currentSequence << std::endl;
                                    //     return;
                                    // }
                                    allSequenceNodeMutex.lock();
                                    if(allSequenceNodes.find(std::make_pair(currentStart,
                                        currentSequence)) == allSequenceNodes.end()) {                                    
                                        allSequenceNodes[std::make_pair(currentStart, currentSequence)]
                                            = std::make_pair(autoIncrId, false);
                                        sequenceNodeIds.push_back(autoIncrId);
                                        sequenceStrands.push_back(false);
                                        autoIncrId++;
                                    } else {
                                        sequenceNodeIds.push_back(allSequenceNodes[
                                            std::make_pair(currentStart, currentSequence)].first);
                                        sequenceStrands.push_back(false);
                                    }
                                    allSequenceNodeMutex.unlock();
                                    currentSequence = "";
                                }
                            }
                            for(size_t k = sequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                                currentSequence += sequence[i].first[j].second[k];
                                currentStart = std::make_tuple(i,j,k);
                                if(currentSequence.length() == 32) {
                                    currentSequence = stripGaps(currentSequence);
                                    if(currentSequence.length()) {
                                        // Since the GFA stores the strand parameter, the reverse
                                        // complement will be computed anyway
                                        std::reverse(currentSequence.begin(), currentSequence.end());
                                        // if(u.first == "NZ_CP006027.1") {
                                        //     std::cout << currentSequence << std::endl;
                                        //     return;
                                        // }
                                        allSequenceNodeMutex.lock();
                                        if(allSequenceNodes.find(std::make_pair(currentStart,
                                            currentSequence)) == allSequenceNodes.end()) {
                                            allSequenceNodes[std::make_pair(currentStart,
                                                currentSequence)] = std::make_pair(autoIncrId, false);
                                            sequenceNodeIds.push_back(autoIncrId);
                                            sequenceStrands.push_back(false);
                                            autoIncrId++;
                                        } else {
                                            sequenceNodeIds.push_back(allSequenceNodes[
                                                std::make_pair(currentStart, currentSequence)].first);
                                            sequenceStrands.push_back(false);
                                        }
                                        allSequenceNodeMutex.unlock();
                                    }
                                    currentSequence = "";
                                }
                            }
                        }
                        if(currentSequence.length()) {
                            currentSequence = stripGaps(currentSequence);
                            if(currentSequence.length()) {
                                // Since the GFA stores the strand parameter, the reverse
                                // complement will be computed anyway
                                std::reverse(currentSequence.begin(), currentSequence.end());
                                allSequenceNodeMutex.lock();
                                if(allSequenceNodes.find(std::make_pair(currentStart, currentSequence))
                                    == allSequenceNodes.end()) {
                                    allSequenceNodes[std::make_pair(currentStart, currentSequence)]
                                        = std::make_pair(autoIncrId, false);
                                    sequenceNodeIds.push_back(autoIncrId);
                                    sequenceStrands.push_back(false);
                                    autoIncrId++;
                                } else {
                                    sequenceNodeIds.push_back(allSequenceNodes[
                                        std::make_pair(currentStart, currentSequence)].first);
                                    sequenceStrands.push_back(false);
                                }
                                allSequenceNodeMutex.unlock();
                            }
                            currentSequence = "";
                        }
                    }
                }
            }

            paths[u.first] = sequenceNodeIds;
            strandPaths[u.first] = sequenceStrands;
        // });
        }

        std::map< std::pair< size_t, bool >, std::string > finalNodes;
        for(const auto& u: allSequenceNodes) {
            finalNodes[u.second] = u.first.second;
        }

        // Graph and its transpose
        // Maps from < blockID, strand > pair to list of neighbours stored as
        // < sequenceId, < blockId, strand > >
        std::map< std::pair< size_t, bool >, std::vector< std::pair< std::string,
            std::pair< size_t, bool > > > > G;
        std::map< std::pair< size_t, bool >, std::vector< std::pair< std::string,
            std::pair< size_t, bool > > > > GT;
        // std::vector< std::pair< std::string, size_t > > G[autoIncrId];
        // std::vector< std::pair< std::string, size_t > > GT[autoIncrId];

        for(const auto& u: paths) {
            for(size_t i = 1; i < u.second.size(); i++) {
                G[std::make_pair(u.second[i-1], strandPaths[u.first][i-1])].push_back(
                    std::make_pair(u.first, std::make_pair(u.second[i], strandPaths[u.first][i])));
                GT[std::make_pair(u.second[i], strandPaths[u.first][i])].push_back(
                    std::make_pair(u.first, std::make_pair(u.second[i-1],
                    strandPaths[u.first][i-1])));
            }
        }

        for(auto& u: G) {
            // Sort so we can compare the sequence IDs of incoming and outgoing edges for equality
            tbb::parallel_sort(u.second.begin(), u.second.end());
        }
        for(auto& u: GT) {
            // Sort so we can compare the sequence IDs of incoming and outgoing edges for equality
            tbb::parallel_sort(u.second.begin(), u.second.end());
        }

        for(auto& u: G) {
            if(finalNodes.find(u.first) == finalNodes.end()) {
                continue;
            }

            while(true) {
                // check if a node's edges only go to one next next node and there is an
                // outgoing edge for every incoming edge
                if(u.second.size() != GT[u.first].size() || u.second.size() == 0) {
                    break;
                }
                bool check = true;

                for(size_t j = 0; j < u.second.size(); j++) {
                    if(u.second[j].first != GT[u.first][j].first) {
                        check = false;
                        break;
                    } else if(j>0 && u.second[j].second != u.second[j-1].second) {
                        check = false;
                        break;
                    }
                }

                if(!check) {
                    break;
                }

                std::pair< size_t, bool > dest = u.second[0].second;
                if(u.second.size() != GT[dest].size()) {
                    break;
                }
                // Combine only if strands match
                if(u.first.second != dest.second) {
                    break;
                }
                for(size_t j = 0; j < u.second.size(); j++) {
                    if(u.second[j].first != GT[dest][j].first && GT[dest][j].second != u.first) {
                        check = false;
                        break;
                    }
                }
                if(G[dest].size() != GT[dest].size()) {
                    break;
                }
                for(size_t j = 0; j < GT[dest].size(); j++) {
                    if(G[dest][j].first != GT[dest][j].first) {
                        check = false;
                        break;
                    }
                }
                if(!check) {
                    break;
                }

                // combine src and dest
                if(u.first.second){
                    // forward strand
                    finalNodes[u.first] += finalNodes[dest];
                } else {
                    // reverse strand
                    finalNodes[u.first] = finalNodes[dest] + finalNodes[u.first];
                }
                finalNodes.erase(dest);
                u.second.clear();
                u.second = G[dest];
            }
        }

        std::set< std::pair<std::pair< size_t, bool >, std::pair< size_t, bool > > > edges;

        for(const auto& u: G) {
            if(finalNodes.find(u.first) == finalNodes.end()) {
                continue;
            }
            for(auto edge: u.second) {
                edges.insert(std::make_pair(u.first, edge.second));
            }
        }

        for(auto& p: paths) {
            std::vector< size_t > newPath;
            std::vector< bool > newStrandPath;
            for(size_t j = 0; j < p.second.size(); j++) {
                if(finalNodes.find(std::make_pair(p.second[j], strandPaths[p.first][j]))
                    != finalNodes.end()) {
                    newPath.push_back(p.second[j]);
                    newStrandPath.push_back(strandPaths[p.first][j]);
                }
            }

            if(currentRefString != currentAltString){
                vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                // Reset
                diffStart = referenceSequence.size();
                currentRefString = "";
                currentAltString = currentRefString;
            }
        }
    }

    std::map< std::string, size_t > sequenceIds;
    for(auto u: allNodes){
        if(u.second->children.size() == 0 && u.first != reference){
            sequenceIds[u.first] = 0;
        }
    }


    fout << "##fileformat=VCFv" << VCF_VERSION << '\n';
    fout << "##fileDate=" << getDate() << '\n';
    fout << "##source=PanMATv" << PMAT_VERSION << '\n';
    fout << "##reference=" << reference << '\n';
    fout << std::left << std::setw(20) << "#CHROM " << std::setw(20) << "POS " << std::setw(20) << "ID " << std::setw(20) << "REF " << std::setw(20) << "ALT " << std::setw(20) << "QUAL " << std::setw(20) << "FILTER " << std::setw(20) << "INFO " << std::setw(20) << "FORMAT ";
    for(auto u: sequenceIds){
        fout << std::left << std::setw(20) << u.first + " ";
    }
    fout << '\n';

    for(auto u: vcfMap){
        for(auto v: u.second){
            if(v.first == ""){
                fout << std::left << std::setw(20) << ". " << std::setw(20) << u.first << " " << std::setw(20) << recordID++ << " " << std::setw(20) << ". ";
            } else {
                fout << std::left << std::setw(20) << ". " << std::setw(20) << u.first << " " << std::setw(20) << recordID++ << " " << std::setw(20) << v.first << " ";
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
            // altStrings += "\t.\t.\t.\t.\t";

            fout << std::left << std::setw(20) << altStrings << " " << std::setw(20) << ". " << std::setw(20) << ". " << std::setw(20) << ". " << std::setw(20) << ". ";

            for(auto w: tempSequenceIds){
                fout << std::left << std::setw(20) << w.second << " ";
            }

            fout << '\n';
        }
    }
}

void PangenomeMAT::Tree::printVCFParallel(std::string reference, std::ofstream& fout){

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "No such leaf node found."){
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
    fout << "##fileDate=" << getDate() << '\n';
    fout << "##source=PanMATv" << PMAT_VERSION << '\n';
    fout << "##reference=" << reference << '\n';
    fout << std::left << std::setw(20) << "#CHROM " << std::setw(20) << "POS " << std::setw(20) << "ID " << std::setw(20) << "REF " << std::setw(20) << "ALT " << std::setw(20) << "QUAL " << std::setw(20) << "FILTER " << std::setw(20) << "INFO " << std::setw(20) << "FORMAT ";
    for(auto u: sequenceIds){
        fout << std::left << std::setw(20) << u.first + " ";
    }
    fout << '\n';

    for(auto u: vcfMap){
        for(auto v: u.second){
            if(v.first == ""){
                fout << std::left << std::setw(20) << ". " << std::setw(20) << u.first << " " << std::setw(20) << recordID++ << " " << std::setw(20) << ". ";
            } else {
                fout << std::left << std::setw(20) << ". " << std::setw(20) << u.first << " " << std::setw(20) << recordID++ << " " << std::setw(20) << v.first << " ";
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
            // altStrings += "\t.\t.\t.\t.\t";

            fout << std::left << std::setw(20) << altStrings << " " << std::setw(20) << ". " << std::setw(20) << ". " << std::setw(20) << ". " << std::setw(20) << ". ";

            for(auto w: tempSequenceIds){
                fout << std::left << std::setw(20) << w.second << " ";
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["insert"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["insert"][(int)j][1].size(); k++) {
                std::string mutationString = pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                insertions[blockId][seqName][number].push_back( std::make_tuple( pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][0].asInt(), pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][1].asInt(), mutationString ) );
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["delete"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["delete"][(int)j][1].size(); k++) {
                deletions[blockId][seqName][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][0].asInt(), pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][1].asInt() ) );
            }
        }
        
    }

    // Rotation
    // Testing data structure 
    std::unordered_map<std::string, std::vector<std::string>> test;
    if (circular)
    {

        std::vector<std::string> sample_base = {};
        int seq_count = 0;
        std::string sample_base_string;
        
        std::vector<std::string> sample_new = {};
        for(const auto& p: paths) 
        {
            // std::cout << p.first << "\n";
            
            test[p.first] = p.second;
            if (seq_count == 0)
            {
                std::unordered_map< std::string, size_t > baseBlockNumber;
                sequenceInverted[p.first] = false;
                rotationIndexes[p.first] = 0;

                for(const auto& block: p.second)
                {
                    blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                    baseBlockNumber[block]++;
                    sample_base.push_back(block);
                }
            }
            else
            {
                // Assigning block numbers
                std::unordered_map< std::string, size_t > baseBlockNumber;
                for(const auto& block: p.second)
                {
                    blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                    baseBlockNumber[block]++;
                }

                std::vector<std::string> sample_dumy = {};
                sample_new.clear();
                for(const auto& block: p.second)
                {
                    // std::cout << block << ",";
                    sample_dumy.push_back(block);
                }
                int rotation_index;
                bool invert = false;
                sample_new= rotate_sample(sample_base, sample_dumy, strandPaths[p.first], blockNumbers[p.first], blockSizeMap, rotation_index, invert);

                std::cout << p.first << "\n";
                // std::vector<string> temp1({"a","b","c","d","e","f"});
                // std::vector<string> temp2({"a","b","c","d","g","h"});
                // std::vector<int> temp3({1,1,1,1,1,1});
                // int temp4;
                // bool temp5;
                // temp2 = rotate_sample(temp1, temp2, temp3, blockSizeMap, temp4, temp5);
                // std::cout << "ROTATED" << std::endl;
                // for(auto u: temp2) {
                //     std::cout << u << " ";
                // }
                // std::cout << std::endl;

                // // Testing
                // for (auto i = 0; i < sample_dumy.size(); i++)
                // {
                //     if (sample_dumy[(i+rotation_index)%sample_dumy.size()] != sample_new[i])
                //     {
                //         std::cout << "Error\n";
                //         // break;
                //     }
                // }

                sequenceInverted[p.first] = invert;
                rotationIndexes[p.first] = rotation_index;

                paths[p.first] = sample_new;
            }
            seq_count++;
        }

        std::cout << "All Seqeunces Rotated\n";
    } else {
        for(auto p: paths) {
            std::unordered_map< std::string, size_t > baseBlockNumber;
            sequenceInverted[p.first] = false;
            rotationIndexes[p.first] = 0;

            for(const auto& block: p.second)
            {
                blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                baseBlockNumber[block]++;
            }
        }
    }

    // Auto increment ID to assign to nodes
    numNodes = 0;

    // Old string Node ID to new integer Node ID
    std::unordered_map< std::string, size_t > stringToNodeId;
    std::unordered_map< std::string, std::vector< size_t > > stringToNodeIds;

    std::unordered_map<int,std::string> intToString; // Locally stored -> Later on mapped to intIdToStringId
    int seqCount = 0; // Current Sequence ID
    std::vector<std::string> consensus = {};
    std::vector<std::string> sample = {};
    std::vector<std::string> consensus_new = {};
    std::vector<int> intSequenceConsensus={};
    std::vector<int> intSequenceSample={};
    std::vector<int> intSequenceConsensus_new={};

    for(const auto& p: paths) 
    {
        if (seqCount == 0)// Load first sequence path
        {
            for(const auto& block: p.second)
            {
                consensus.push_back(block);
                // sample_base.push_back(block);
                intToString[numNodes] = block;
                intSequences[p.first].push_back(numNodes);
                intSequenceConsensus.push_back(numNodes);
                numNodes++;
            }
        }
        else
        {
            intSequenceSample.clear();
            intSequenceConsensus_new.clear();
            sample.clear();
            consensus_new.clear();

            for(const auto& block: p.second)
            {
                sample.push_back(block);
                // std::cout << block << " ";
            }
            // std::cout << "\n";

            chain_align (consensus, 
                sample, 
                intSequenceConsensus,
                intSequenceSample,
                numNodes, 
                consensus_new,
                intSequenceConsensus_new,
                intToString);
            for (auto &b: intSequenceSample)
            {
                intSequences[p.first].push_back(b);
            }
            consensus.clear();
            intSequenceConsensus.clear();
            for (auto &b: consensus_new)
            {
                consensus.push_back(b);
            }

            for (auto &b: intSequenceConsensus_new)
            {
                intSequenceConsensus.push_back(b);
            }
            
        }
        seqCount++;
        std::cout << seqCount << " " << intSequenceConsensus_new.size() << std::endl;
    }

    // // re-assigning IDs in fixed order
    int reorder = 0;
    std::unordered_map<int,int> order_map = {};
    for (auto &i: intSequenceConsensus)
    {
        order_map[i] = reorder;
        intIdToStringId[reorder] = intToString[i];
        topoSortedIntSequences.push_back(reorder);
        reorder++;
    }

    for (auto &m: intSequences)
    {
        for (auto &s: m.second)
        {
            s = order_map[s];
        }
    }

}

std::unordered_map< std::string,std::vector< int > > PangenomeMAT::Pangraph::getAlignedStrandSequences(const std::vector< size_t >& topoArray) {
    std::unordered_map< std::string, std::vector< int > > alignedStrandSequences;
    for(auto p: intSequences) {
        size_t p1 = 0, p2 = 0;
        while(p1 < topoArray.size() && p2 < p.second.size()) {
            if(topoArray[p1] == p.second[p2]) {
                alignedStrandSequences[p.first].push_back(strandPaths[p.first][p2]);
                p2++;
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

/* Syncmer indexing */

// Read alignment 
extern "C" {
    void align_reads(const char *reference, int n_reads, const char **reads, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions, size_t seed_k, size_t seed_s);
}

void PangenomeMAT::Tree::setupGlobalCoordinates(){
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


size_t PangenomeMAT::Tree::getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition, int nucGapPosition){
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


// Helpers
std::set<kmer_t> PangenomeMAT::syncmersFromFastq(std::string fastqPath,  std::vector<read_t> &reads, size_t k, size_t s) {
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
    float est_coverage = 1; //TODO change this to 1
    bool open = false;
    
    std::set<kmer_t> syncmers;
    std::unordered_map<std::string, int> counts;
    std::unordered_map<std::string, int> counts_rc;

    std::cerr << "length: " << input.size() << "\n";
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

 //               m.pos = m.pos + k - 1;
//                m.reversed = false;
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, false});
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
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, true});
            }
        }
        reads[i] = this_read;
    }
   
    return syncmers;
}

template <typename INT, typename T> void PangenomeMAT::removeIndices(std::vector<T>& v, std::stack<INT>& rm)
{
    if (rm.size() < 1) {
        return;
    }
  int32_t rmVal = rm.top();
  rm.pop();
  v.erase(
    std::remove_if(std::begin(v), std::end(v), [&](T& elem)
    {
        if (rmVal == -1) {
            return false;
        }
        if (&elem - &v[0] == rmVal) {
            if (!rm.empty()) {
                rmVal = rm.top();
                rm.pop();
            } else {
                rmVal = -1;
            }
            return true;
        }
        return false;

    }),
    std::end(v)
  );
}

std::pair<int32_t, int32_t> PangenomeMAT::getRecomputePositions(std::pair<int32_t, int32_t> p, std::string &gappedSequence, int32_t k) {
    int32_t mutPos = p.first;
    int32_t mutLen = p.second;

    int32_t numSeen = 0;
    int32_t start = mutPos;

    while(numSeen < k && start > 0) {
        if (gappedSequence[start] != '-') {
            numSeen++;
            if (numSeen == k) {
                break;
            }
        }
        start--;
    }

    numSeen = 0;

    int32_t stop = mutPos + mutLen - 1;
    while(numSeen < k && stop < gappedSequence.size()) {
        if (gappedSequence[stop] != '-') {
            numSeen++;
            if (numSeen == k) {
                break;
            }
        }
        stop++;
    }
    stop = std::min(stop, (int32_t) gappedSequence.size() - 1);

    return std::make_pair(start, stop);
}

int32_t PangenomeMAT::alignedEndPos(int32_t pos, int32_t k, std::string &gappedSequence) {
    int32_t numSeen = 0;
    int32_t stop = pos;
    while(numSeen < k && stop < gappedSequence.size()) {
        if (gappedSequence[stop] != '-') {
            numSeen++;
            if (numSeen == k) {
                break;
            }
        }
        stop++;
    }
    stop = std::min(stop, (int32_t) gappedSequence.size() - 1);
    return stop;

}

void updateJaccard(dynamicJaccard &dj, std::unordered_map<std::string, bool> &readSyncmers, std::vector<kmer_t> &deletedSyncmers, std::vector<kmer_t> &insertedSyncmers) {
    for (const kmer_t &syncmer : deletedSyncmers) {
        //std::cout << "jac del: " << syncmer.seq << "\n";
        if (readSyncmers.find(syncmer.seq) != readSyncmers.end()) {
            dj.intersectionSize -= 1;
        } else {
            dj.unionSize -= 1;
        }
    }
    for (const kmer_t &syncmer : insertedSyncmers) {
        //std::cout << "jac ins: " << syncmer.seq << "\n";
        if (readSyncmers.find(syncmer.seq) != readSyncmers.end()) {
            dj.intersectionSize += 1;
        } else {
            dj.unionSize += 1;
        }
    }
    dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
    //std::cout << "inter: " << dj.intersectionSize << " union: " << dj.unionSize << " ji: " << dj.jaccardIndex << "\n";

}

struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

// Index construction
void PangenomeMAT::discardSyncmers(std::vector<kmer_t> &inSyncmers, const std::vector<std::pair<int32_t, int32_t>>& B, std::string &gappedSequence, std::unordered_map<std::string, kmer_t> &to_insert, std::unordered_map<std::string, bool> &variable_syncmers, seedIndex &index, std::string nid, size_t k) {
    //std::cout << "discardSyncmers\n";
    for (int32_t i = inSyncmers.size() - 1; i >= 0; i--) {
        
        const kmer_t s = inSyncmers[i];
    //    std::cout << "checking " << s.seq << "\n";
        for (const auto& b : B) {
            // std::cout << "range " << b.first << " " << b.second << "\n";
            // std::cout << "start pos " << s.pos << "\n";
            // std::cout << "end pos " << PangenomeMAT::alignedEndPos(s.pos, k, gappedSequence) << "\n";

            //todo make s.end a thing
            if (s.pos >= b.first && PangenomeMAT::alignedEndPos(s.pos, k, gappedSequence) <= b.second) {

                auto it = to_insert.find(s.seq);
          
                if (it == to_insert.end()) {
                 //   std::cout << "can delete " << s.seq << "\n";
                    index.deletions[nid].push_back(kmer_t{s.seq, s.pos, i}); 
                    variable_syncmers[s.seq] = true; // track variant sites
                    inSyncmers.erase(inSyncmers.begin() + i);
                } else {
                //    std::cout << "dont insert " << s.seq << "\n";
                    to_insert.erase(it->first);
                }
                break;
            }
        }
    }
}

void PangenomeMAT::Tree::condenseTree(PangenomeMAT::Node *node) {
    for (auto child: node->children) {
        condenseTree(child);
    }
    if (node->nucMutation.size() == 0 && node->blockMutation.size() == 0 && node->parent != root) {
        node->parent->children.erase(std::remove(node->parent->children.begin(), node->parent->children.end(), node), node->parent->children.end());
        node->parent->identifier += "_" + node->identifier;
        for (auto child: node->children) {
            child->parent = node->parent;
            node->parent->children.push_back(child);
        }
        delete node;
    }
}


void PangenomeMAT::indexSyncmersHelper(PangenomeMAT::Tree *T, PangenomeMAT::Node *root,\
    sequence_t &sequence, blockExists_t &blockExists, blockStrand_t &blockStrand, seedIndex &index, std::vector<kmer_t> &syncmers,\
    std::unordered_map<std::string, int32_t> &counts, std::unordered_map<std::string, bool> &variable_syncmers, size_t k, size_t s){
    

    // Apply mutations    
    // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
    std::vector< std::tuple< int32_t, int32_t, bool, bool, bool, bool > > blockMutationInfo;
    std::vector<std::pair<int32_t, int32_t>> mutPositions;
    // Block Mutations
    for(auto mutation: root->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        if(type == 1) {
            // insertion

            bool oldStrand;
            bool oldMut;
            if(secondaryBlockId != -1) {
                oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                blockExists[primaryBlockId].second[secondaryBlockId] = true;

                // if insertion of inverted block takes place, the strand is backwards
                blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
            } else {
                oldStrand = blockStrand[primaryBlockId].first;
                oldMut = blockExists[primaryBlockId].first;
                blockExists[primaryBlockId].first = true;

                // if insertion of inverted block takes place, the strand is backwards
                blockStrand[primaryBlockId].first = !inversion;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion) );
        } else {
            bool oldMut;
            bool oldStrand;
            if(inversion) {
                // This means that this is not a deletion, but instead an inversion
                if(secondaryBlockId != -1) {
                    oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                    oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                    blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
                } else {
                    oldStrand = blockStrand[primaryBlockId].first;
                    oldMut = blockExists[primaryBlockId].first;
                    blockStrand[primaryBlockId].first = !oldStrand;
                }
                if(oldMut != true) {
                    std::cout << "Problem in PanMAT generation" << std::endl;
                }
                blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
            } else {
                // Actually a deletion

                if(secondaryBlockId != -1) {
                    oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
                    oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
                    blockExists[primaryBlockId].second[secondaryBlockId] = false;

                    // resetting strand to true during deletion
                    blockStrand[primaryBlockId].second[secondaryBlockId] = true;
                } else {
                    oldStrand = blockStrand[primaryBlockId].first;
                    oldMut = blockExists[primaryBlockId].first;
                    blockExists[primaryBlockId].first = false;

                    // resetting strand to true during deletion
                    blockStrand[primaryBlockId].first = true;
                }
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true) );
        }

    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++) {
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';
        size_t globalCoord = T->getGlobalCoordinate(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition);


        if(type < 3) {
            // Either S, I or D

            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(type == PangenomeMAT::NucMutationType::NS) {
                // Substitution
                mutPositions.push_back(std::make_pair(globalCoord, len));
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));   
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NI) {
                // Insertion
                mutPositions.push_back(std::make_pair(globalCoord, len));
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            const int nucCode = ((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF;
                            //std::cout << "nucCode: " << nucCode << std::endl;
                            newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));   
                        }
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::ND) {
                // Deletion
                mutPositions.push_back(std::make_pair(globalCoord, 0));
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        }
                    }
                }
            }
        } 
        else {
            if(type == PangenomeMAT::NucMutationType::NSNPS) {
                // SNP Substitution
                mutPositions.push_back(std::make_pair(globalCoord, 0));
                newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1) {
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
            else if(type == PangenomeMAT::NucMutationType::NSNPI) {
                // SNP Insertion
                mutPositions.push_back(std::make_pair(globalCoord, 1));
                newVal = PangenomeMAT::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    }
                } else {
                    if(nucGapPosition != -1) {
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
            else if(type == PangenomeMAT::NucMutationType::NSNPD) {
                // SNP Deletion
                mutPositions.push_back(std::make_pair(globalCoord, 0));
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    } else {
                        char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
                        sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    }
                } else {
                    if(nucGapPosition != -1) {
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

   
    // Aligned sequence of the current node
    std::string currNodeSequence = getCurrentFastaSequence(sequence, blockExists, blockStrand, true);
	//std::coutt << currNodeSequence << "\n";

    std::vector<std::pair<int32_t, int32_t>> recompute;
    //std::coutt << "mut positions: \n";
    // for (const std::pair<int32_t, int32_t> &p : mutPositions) {
    //     //std::coutt << p.first << ", " << p.second << "\n";
    // }
    for (const std::pair<int32_t, int32_t> &p : mutPositions) {
        auto r = getRecomputePositions(p, currNodeSequence, k);
        recompute.push_back(r);
    }

    std::unordered_map<std::string, kmer_t> to_insert;
    	////std::coutt << "\n***** " << root->identifier << "\n";

    for (auto &range : recompute) {
        if (range.first >= range.second) {
            continue;
        }
      //  std::cout << "redo range: " << range.first << ", " << range.second << "\n";

        int32_t seqLen = currNodeSequence.size();
        //  xx(xxxxxxx)
        //  01.23456789
        //    
        std::string redo = currNodeSequence.substr(std::max(0, range.first), std::min(seqLen - range.first, 1 + range.second - range.first)); 
        auto redone = syncmerize(redo, k, s, false, true, std::max(0,range.first));
        for (const kmer_t &syncmer : redone) {
     //       	std::cout << " => redone: " << syncmer.seq << " " << syncmer.pos << " idx: " << syncmer.idx << "\n";
            to_insert[syncmer.seq] = syncmer;
        }
    }
    PangenomeMAT::discardSyncmers(syncmers, recompute, currNodeSequence, to_insert, variable_syncmers, index, root->identifier, k); //modifies mutated_syncmers
    for (auto &s : to_insert) {
        syncmers.push_back(s.second);
        index.insertions[root->identifier].push_back(s.second);
    }

     //std::cout << root->identifier << "\t";
    //  for (kmer_t s : syncmers) {
    //      std::cout << s.seq << "\t";
    //  }
    //  std::cout << "\n";
   for(PangenomeMAT::Node* child: root->children){
        indexSyncmersHelper(T, child, sequence, blockExists, blockStrand, index, syncmers, counts, variable_syncmers, k, s);
    }
    
    	//std::cout << "undo by deleting last " << index.insertions[root->identifier].size() << "\n";
    	//std::cout << "size before " << syncmers.size() << "\n";
    syncmers.erase(syncmers.end() - index.insertions[root->identifier].size(), syncmers.end());
    	//std::cout << "size after " << syncmers.size() << "\n";

    for (int32_t i = index.deletions[root->identifier].size() - 1; i >= 0; i--) {
        	//std::cout << "undo by inserting " << index.deletions[root->identifier][i].seq << " at " << index.deletions[root->identifier][i].idx << "\n";
        syncmers.insert(syncmers.begin() + index.deletions[root->identifier][i].idx, index.deletions[root->identifier][i]);
    }


    // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<3>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
        }
    }

    // Undo nuc mutations when current node and its subtree have been processed
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }
}

void PangenomeMAT::indexSyncmers(PangenomeMAT::Tree *T, std::ofstream& fout, size_t k, size_t s) {

   // T->condenseTree(T->root);

    PangenomeMAT::Node *root = T->root;
    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    sequence_t sequence(T->blocks.size() + 1);
    BlockGapList blockGaps = T->blockGaps;
    std::vector< GapList > gaps = T->gaps;
    std::vector< Block > blocks = T->blocks;
    blockExists_t blockExists(blocks.size() + 1, {false, {}});
    blockStrand_t blockStrand(blocks.size() + 1, {true, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        blockStrand[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], true);
    }

    int32_t maxBlockId = 0;

    // Create consensus sequence of blocks
    for(size_t i = 0; i < blocks.size(); i++) {

        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = PangenomeMAT::getNucleotideFromCode(nucCode);
                
                if(secondaryBlockId != -1) {
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }

            if(endFlag) {
                break;
            }
        }

        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1) {
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    // Assigning nucleotide gaps in blocks
    for(size_t i = 0; i < gaps.size(); i++) {
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            if(secondaryBId != -1) {
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }
    
    blockExists_t blockExistsAllOn(blockExists.size(), {true, {}});   
    std::string consensusSequence = getCurrentFastaSequence(sequence, blockExistsAllOn, blockStrand, true);

    //std::cout << "consensus sequence: " << consensusSequence << "\n";
    std::vector<kmer_t> initialSyncmers;
    std::unordered_map<std::string, int> counts;


    //TODO this is bad
    initialSyncmers = syncmerize(consensusSequence, k, s, false, true, 0);
    std::set<kmer_t> seedSet = std::set<kmer_t>(initialSyncmers.begin(), initialSyncmers.end());
    std::vector<kmer_t> dynamicSeeds = std::vector<kmer_t>(seedSet.begin(), seedSet.end());

    seedIndex index;
    std::unordered_map<std::string, bool> variable_syncmers;

    std::cout << "indexing syncmers...\n";
    // main DFS
    indexSyncmersHelper(T, T->root, sequence, blockExists, blockStrand, index, dynamicSeeds, counts, variable_syncmers, k, s);

    std::unordered_map<std::string, bool> invariants;
    std::vector<kmer_t> initialShaved;
    for (const kmer_t &syncmer : dynamicSeeds) {
        //std::cout << syncmer.seq << "\n";
        if (variable_syncmers.find(syncmer.seq) == variable_syncmers.end()) {
            invariants[syncmer.seq] = true;
        } else {
            initialShaved.push_back(kmer_t{syncmer.seq, syncmer.pos, -1});
        }
    }
    std::cout << "total consensus seeds: " << dynamicSeeds.size() << "\n";
    std::cout << "# invariant seeds: " << invariants.size() << "\n";
    index.rootSeeds = dynamicSeeds;
    seedIndex shaved = PangenomeMAT::shaveIndex(T->root, index, invariants, initialShaved);
    shaved.rootSeeds = initialShaved;
    PangenomeMAT::writeIndex(T->root, fout, index);
}

void PangenomeMAT::shaveDFS(PangenomeMAT::Node *currNode, std::vector<kmer_t> &currNodeSyncmers, seedIndex &index_orig, std::unordered_map<std::string, std::vector<kmer_t>> &deletions, std::unordered_map<std::string, bool> &invariants) {
    
    std::vector<int32_t> di;
    // descending order of idx
    for (kmer_t d : index_orig.deletions[currNode->identifier]) {
        auto r = std::find(currNodeSyncmers.begin(), currNodeSyncmers.end(), d);
        if (r != currNodeSyncmers.end()) {
            deletions[currNode->identifier].push_back(kmer_t{r->seq, r->pos, static_cast<int32_t>(r - currNodeSyncmers.begin())});
            di.push_back(static_cast<int32_t>(r - currNodeSyncmers.begin()));
        }
    }

    // output index dels in descending order
    std::sort(deletions[currNode->identifier].begin(), deletions[currNode->identifier].end(),
         [](const kmer_t &x, const kmer_t &y){ return (x.idx > y.idx); });

    std::sort(di.begin(), di.end());

    std::stack<int32_t> delIndices;
    for (auto i = di.begin(); i != di.end(); i++) {
        delIndices.push(*i);
    }
    
    PangenomeMAT::removeIndices(currNodeSyncmers, delIndices);

    for (const kmer_t &s : index_orig.insertions[currNode->identifier]) {
        currNodeSyncmers.push_back(s);
    }
    
    for (Node *child : currNode->children) {
        PangenomeMAT::shaveDFS(child, currNodeSyncmers, index_orig, deletions, invariants);
    }

    currNodeSyncmers.erase(currNodeSyncmers.end() - index_orig.insertions[currNode->identifier].size(), currNodeSyncmers.end());
    
    for (int32_t i = deletions[currNode->identifier].size() - 1; i >= 0; i--) {
        currNodeSyncmers.insert(currNodeSyncmers.begin() + deletions[currNode->identifier][i].idx, deletions[currNode->identifier][i]);
    }
    
}

seedIndex PangenomeMAT::shaveIndex(PangenomeMAT::Node *root, seedIndex &index, std::unordered_map<std::string, bool> &invariants, std::vector<kmer_t> &initialSyncmers) {
    seedIndex shaved;
    shaved.insertions = index.insertions;
    std::unordered_map<std::string, std::vector<kmer_t>> deletions;
    PangenomeMAT::shaveDFS(root, initialSyncmers, index, deletions, invariants);
    shaved.deletions = deletions;
    return shaved;
}

// Sample placement
void PangenomeMAT::placeDFS(Node *currNode, std::vector<kmer_t> &currNodeSyncmers, std::unordered_map<std::string, bool> &querySyncmers, seedIndex &index, dynamicJaccard dj, std::unordered_map<std::string, float> &scores) {
    

    std::stack<int32_t> delIndices;
    // TODO make not silly
    for (auto i = index.deletions[currNode->identifier].begin(); i != index.deletions[currNode->identifier].end(); i++) {
        delIndices.push(i->idx);
        i->seq = currNodeSyncmers[i->idx].seq;
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

void PangenomeMAT::placeSample(PangenomeMAT::Tree *T, std::string fastqPath, seedIndex &index, size_t k, size_t s){
    //T->condenseTree(T->root);
    
    PangenomeMAT::Node *root = T->root;
    std::vector<read_t> reads;

    auto fastq_start = std::chrono::high_resolution_clock::now();
    std::set<kmer_t> readSyncmers = syncmersFromFastq(fastqPath, reads, k, s);
    auto fastq_end = std::chrono::high_resolution_clock::now();

    std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

    auto place_start = std::chrono::high_resolution_clock::now();

    std::set<kmer_t> rootSyncmers = std::set<kmer_t>(index.rootSeeds.begin(), index.rootSeeds.end());

    std::cerr << "\n";
    std::cerr << "Placing sample...\n";


    struct dynamicJaccard dj;
 
    dj.intersectionSize = intersection_size(rootSyncmers, readSyncmers);
    dj.unionSize = rootSyncmers.size() + readSyncmers.size() - dj.intersectionSize;
    dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
    
    
    std::cout << "root seeds: " << rootSyncmers.size() << "\n";
    std::cout << "read seeds: " << readSyncmers.size() << "\n";
    for (const auto &k : readSyncmers) {
        std::cout << k.seq << "\n";
    }
    std::cout << "initial jaccard: " << dj.jaccardIndex << "\n";

    std::unordered_map<std::string, float> scores;
    std::unordered_map<std::string, bool> readSyncmersMap;
    for (const auto &k : readSyncmers) {
        readSyncmersMap[k.seq] = true;
    }

    placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);

    auto place_end = std::chrono::high_resolution_clock::now();

    std::cout << "place time: " << std::chrono::duration_cast<std::chrono::milliseconds>(place_end - place_start).count() << "\n";


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
    auto aln_start = std::chrono::high_resolution_clock::now();

    std::string ref_seq = T->getStringFromReference(best_match, false);


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

//    std::cout << "\n" << ref_seq << std::endl;


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
    align_reads(reference, n_reads, read_strings, r_lens, seed_counts, reversed, ref_positions, qry_positions, k, s);


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
    auto aln_end = std::chrono::high_resolution_clock::now();
    std::cout << "aln time: " << std::chrono::duration_cast<std::chrono::milliseconds>(aln_end - aln_start).count() << "\n";

}

// Index saving and loading
void PangenomeMAT::writeIndexDFS(Node *currNode, seedIndex &index, std::stringstream &ss, std::vector<kmer_t> &seeds) {
    ss << currNode->identifier << "\t";
  

    std::stack<int32_t> delIndices;
    // TODO make not silly
    for (auto i = index.deletions[currNode->identifier].begin(); i != index.deletions[currNode->identifier].end(); i++) {
        delIndices.push(i->idx);
    }
    for (const kmer_t &s : index.deletions[currNode->identifier]) {
        ss << s.idx << " ";
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

void PangenomeMAT::writeIndex(PangenomeMAT::Node *root, std::ofstream &fout, seedIndex &index) {
    for (const kmer_t &s : index.rootSeeds) {
        fout << s.seq << " ";
    }
    fout << "\n";
    std::stringstream ss;
    writeIndexDFS(root, index, ss, index.rootSeeds);
    fout << ss.str();
}

void PangenomeMAT::loadIndex(PangenomeMAT::Node *root, std::ifstream &indexFile, seedIndex &index) {
    std::string rootSyncmersString;
    std::vector<kmer_t> rootSyncmers;
    std::getline(indexFile, rootSyncmersString);
    std::vector< std::string > rootSplt;
    PangenomeMAT::stringSplit(rootSyncmersString, ' ', rootSplt);
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