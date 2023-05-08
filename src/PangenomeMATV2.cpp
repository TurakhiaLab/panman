#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
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

PangenomeMAT2::Block::Block(MATNew::block b, const std::vector< uint32_t >& blockConsensusSeq){
    primaryBlockId = (b.blockid() >> 32);
    if(b.blockgapexist()){
        secondaryBlockId = (b.blockid() & 0xFFFFFFFF);
    } else {
        secondaryBlockId = -1;
    }
    
    chromosomeName = b.chromosomename();
    consensusSeq = blockConsensusSeq;
    // for(int i = 0; i < b.consensusseq_size(); i++){
    //     consensusSeq.push_back(b.consensusseq(i));
    // }
}

PangenomeMAT2::Block::Block(size_t blockId, std::string seq){
    primaryBlockId = blockId;
    secondaryBlockId = -1;
    for(size_t i = 0; i < seq.length(); i+=8){
        uint32_t currentConsensusSeq = 0;
        for(size_t j = i; j < std::min(i+8, seq.length()); j++){
            int code = PangenomeMAT2::getCodeFromNucleotide(seq[j]);
            currentConsensusSeq^=(code << (4*(7-(j-i))));
        }
        consensusSeq.push_back(currentConsensusSeq);
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
        for(auto nucMut: nodes[currentIndex].mutations(i).nucmutation()){
            storedNucMutation.push_back( PangenomeMAT2::NucMut(nucMut, nodes[currentIndex].mutations(i).blockid(), nodes[currentIndex].mutations(i).blockgapexist()) );
        }
    }

    std::vector< PangenomeMAT2::BlockMut > storedBlockMutation;
    for(int i = 0; i < nodes[currentIndex].mutations_size(); i++){
        PangenomeMAT2::BlockMut tempBlockMut;
        if(nodes[currentIndex].mutations(i).blockmutexist()){
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

std::string generateRandomNewickString(std::vector< std::string > sequenceIds){
    if(sequenceIds.size() == 0){
        return "()";
    }

    while(sequenceIds.size() > 1){
        size_t p1 = rand()%sequenceIds.size();
        size_t p2 = rand()%sequenceIds.size();
        while(p2 == p1){
            p2 = rand()%sequenceIds.size();
        }
        std::string str1 = sequenceIds[p1], str2 = sequenceIds[p2];
        sequenceIds[p1] = "("+sequenceIds[p1]+","+sequenceIds[p2]+")";
        sequenceIds.erase(sequenceIds.begin()+p2);
    }
    if(sequenceIds[0][0] != '('){
        sequenceIds[0] = "(" + sequenceIds[0] + ")";
    }
    return sequenceIds[0];
}

int PangenomeMAT2::Tree::nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states){
    if(node->children.size() == 0){
        if(states.find(node->identifier) == states.end()){
            return states[node->identifier] = 0;
        }
        return states[node->identifier];
    }
    std::vector< int > childStates;
    for(auto child: node->children){
        childStates.push_back(nucFitchForwardPass(child, states));
    }
    int orStates = 0, andStates = childStates[0];
    for(auto u: childStates){
        orStates |= u;
        andStates &= u;
    }
    if(andStates){
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;
}

void PangenomeMAT2::Tree::nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultState){
    if(node == root && defaultState != (1 << 28)){
        states[node->identifier] = defaultState;
    } else {
        if(states[node->identifier] == 0){
            return;
        }
        if(parentState & states[node->identifier]){
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)){
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children){
        nucFitchBackwardPass(child, states, states[node->identifier]);
    }
}

void PangenomeMAT2::Tree::nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > >& mutations, int parentState){
    if(states[node->identifier] == 0){
        return;
    }
    if(parentState != states[node->identifier]){
        if(parentState == 1){
            // insertion
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0){
                currentState >>= 1;
                code++;
            }
            code--;
            // if(code != 1 && code != 2 && code != 4 && code != 8 && code != 15){
            //     std::cout << code << std::endl;
            // }
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 1){
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0){
                currentState >>= 1;
                code++;
            }
            code--;
            // if(code != 1 && code != 2 && code != 4 && code != 8 && code != 15){
            //     std::cout << code << std::endl;
            // }
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }
    for(auto child: node->children){
        nucFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}

int PangenomeMAT2::Tree::blockFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states){
    if(node->children.size() == 0){
        if(states.find(node->identifier) == states.end()){
            std::cerr << "FATAL: State for leaf Node ID not found!" << std::endl;
            exit(-1);
            return 0;
        }
        return states[node->identifier];
    }
    std::vector< bool > stateExists(3,false);
    for(auto child: node->children){
        stateExists[blockFitchForwardPass(child, states)] = true;
    }
    if(stateExists[0] && stateExists[1]){
        return states[node->identifier] = 2;
    }
    if(stateExists[0]){
        return states[node->identifier] = 0;
    }
    if(stateExists[1]){
        return states[node->identifier] = 1;
    }
    return states[node->identifier] = 2;
}

void PangenomeMAT2::Tree::blockFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultValue){
    if(defaultValue != 2){
        states[node->identifier] = defaultValue;
    } else if(states[node->identifier] == 2){
        states[node->identifier] = parentState;
    }
    for(auto child: node->children){
        blockFitchBackwardPass(child, states, states[node->identifier]);
    }
}

void PangenomeMAT2::Tree::blockFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, bool >& mutations, int parentState){
    if(parentState == 0 && states[node->identifier] == 1){
        mutations[node->identifier] = true;
    } else if(parentState == 1 && states[node->identifier] == 0){
        mutations[node->identifier] = false;
    }
    for(auto child: node->children){
        blockFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}

PangenomeMAT2::Tree::Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype){
    if(ftype == PangenomeMAT2::FILE_TYPE::GFA){
        std::map< std::string, std::string > nodes;
        std::map< std::string, std::vector< std::string > > paths;
        // std::map< std::string, std::vector< size_t > > newPaths;
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
        std::vector< std::vector< std::string > > stringSequences;
        std::vector< std::string > sequenceIds;
        for(auto p: paths){
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
        for(size_t i = 0; i < g.pathIds.size(); i++){
            pathIdToSequence[g.pathIds[i]] = alignedSequences[i];
        }

        // std::cout << "Path ID to sequence generated " << topoArray.size() << std::endl;

        for(size_t i = 0; i < topoArray.size(); i++){
            blocks.emplace_back(i, g.intNodeToSequence[topoArray[i]]);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, bool > > globalMutations;


        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i){
            // std::cout << i << std::endl;
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, bool > mutations;
            for(const auto& u: pathIdToSequence){
                states[u.first] = (u.second[i] != -1);
            }
            blockFitchForwardPass(root, states);
            blockFitchBackwardPass(root, states, 1);
            blockFitchAssignMutations(root, states, mutations, 0);
            globalMutations[i] = mutations;
        });

        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes){
            nodeMutexes[u.first];
        }

        tbb::parallel_for_each(globalMutations, [&](auto& pos){
            auto& mutations = pos.second;
            for(const auto& node: allNodes){
                if(mutations.find(node.first) != mutations.end()){
                    nodeMutexes[node.first].lock();
                    node.second->blockMutation.emplace_back(pos.first, mutations[node.first]);
                    nodeMutexes[node.first].unlock();
                }
            }
        });
    } else if(ftype == PangenomeMAT2::FILE_TYPE::PANGRAPH){
        std::string newickString;
        secondFin >> newickString;
        Json::Value pangraphData;
        fin >> pangraphData;
        PangenomeMAT2::Pangraph pg(pangraphData);

        std::vector< size_t > topoArray = pg.getTopologicalSort();
        std::unordered_map< std::string, std::vector< int > > alignedSequences = pg.getAlignedSequences(topoArray);
        // for(int i = 0; i < alignedSequences["ON647447.1"].size(); i++){
        //     std::cout << alignedSequences["ON647447.1"][i] << " ";
        // }
        // std::cout << std::endl;
        
        root = createTreeFromNewickString(newickString);

        for(size_t i = 0; i < topoArray.size(); i++){
            // if(i == 4){
            //     std::cout << pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]] << std::endl;
            // }
            blocks.emplace_back(i, pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]]);
        }

        for(size_t i = 0; i < topoArray.size(); i++){
            GapList g;
            g.primaryBlockId = i;
            g.secondaryBlockId = -1;
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++){
                g.nucPosition.push_back(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first);
                g.nucGapLength.push_back(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second);
            }
            gaps.push_back(g);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, bool > > globalBlockMutations;


        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i){
            // std::cout << i << std::endl;
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, bool > mutations;
            for(const auto& u: alignedSequences){
                states[u.first] = (u.second[i] != -1);
            }
            
            blockFitchForwardPass(root, states);
            blockFitchBackwardPass(root, states, 1);
            blockFitchAssignMutations(root, states, mutations, 0);
            globalBlockMutations[i] = mutations;

        });


        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes){
            nodeMutexes[u.first];
        }

        tbb::parallel_for_each(globalBlockMutations, [&](auto& pos){
            auto& mutations = pos.second;
            for(const auto& node: allNodes){
                if(mutations.find(node.first) != mutations.end()){
                    nodeMutexes[node.first].lock();
                    node.second->blockMutation.emplace_back(pos.first, mutations[node.first]);
                    nodeMutexes[node.first].unlock();
                }
            }
        });

        std::unordered_map< std::string, std::vector< size_t > > blockCounts;
        for(auto u: alignedSequences){
            blockCounts[u.first].resize(u.second.size(), 0);
        }

        tbb::parallel_for_each(alignedSequences, [&](const auto& u){
            std::unordered_map< std::string, size_t > currentCount;
            for(size_t i = 0; i < u.second.size(); i++){
                if(u.second[i] != -1){
                    blockCounts[u.first][i] = currentCount[pg.intIdToStringId[u.second[i]]] + 1;
                    currentCount[pg.intIdToStringId[u.second[i]]]++;
                }
            }
        });

        // tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, bool > > globalBlockMutations;
        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > gapMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i){
            std::string consensusSeq = pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]];
            std::vector< std::pair< char, std::vector< char > > > sequence(consensusSeq.size()+1, {'-', {}});
            for(size_t j = 0; j < consensusSeq.length(); j++){
                sequence[j].first = consensusSeq[j];
            }
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++){
                sequence[pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first].second.resize(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second, '-');
            }
            tbb::concurrent_unordered_map< std::string, std::vector< std::pair< char, std::vector< char > > > > individualSequences;

            tbb::parallel_for_each(alignedSequences, [&](const auto& u){
                if(u.second[i] == -1){
                    return;
                }
                std::vector< std::pair< char, std::vector< char > > > currentSequence = sequence;

                for(const auto& v: pg.substitutions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]){
                    currentSequence[v.first-1].first = v.second[0];
                }
                for(const auto& v: pg.insertions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]){
                    for(size_t j = 0; j < std::get<2>(v).length(); j++){
                        currentSequence[std::get<0>(v)].second[std::get<1>(v)+j] = std::get<2>(v)[j];
                    }
                }
                for(const auto& v: pg.deletions[pg.intIdToStringId[topoArray[i]]][u.first][blockCounts[u.first][i]]){
                    for(size_t j = v.first; j < v.first + v.second; j++){
                        currentSequence[j-1].first = '-';
                    }
                }
                individualSequences[u.first] = currentSequence;
            });
            tbb::parallel_for((size_t) 0, sequence.size(), [&](size_t j){
                tbb::parallel_for((size_t)0, sequence[j].second.size(), [&](size_t k){
                    std::unordered_map< std::string, int > states;
                    std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
                    for(const auto& u: individualSequences){
                        if(u.second[j].second[k] != '-'){
                            states[u.first] = (1 << getCodeFromNucleotide(u.second[j].second[k]));
                            // if(states[u.first] != 2 && states[u.first] != 4 && states[u.first] != 16 && states[u.first] != 256){
                            //     std::cout << "Found something! " << states[u.first] << " " << (int)u.second[j].second[k] << std::endl;
                            // }
                        } else {
                            states[u.first] = 1;
                        }
                    }
                    nucFitchForwardPass(root, states);
                    nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                    nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                    for(auto mutation: mutations){
                        nodeMutexes[mutation.first].lock();
                        gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                        // allNodes[mutation.first]->nucMutation.emplace_back( std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)) );
                        nodeMutexes[mutation.first].unlock();
                    }
                });
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
                for(const auto& u: individualSequences){
                    if(u.second[j].first != '-'){
                        states[u.first] = (1 << getCodeFromNucleotide(u.second[j].first));
                        // if(states[u.first] != 2 && states[u.first] != 4 && states[u.first] != 16 && states[u.first] != 256){
                        //         std::cout << "Found something! " << states[u.first] << " " << (int)u.second[j].first << std::endl;
                        //     }
                    } else {
                        states[u.first] = 1;
                    }
                }
                nucFitchForwardPass(root, states);
                nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].first)));
                nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].first)));
                for(auto mutation: mutations){
                    nodeMutexes[mutation.first].lock();
                    nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                    // allNodes[mutation.first]->nucMutation.emplace_back( std::make_tuple((int)i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)) );
                    nodeMutexes[mutation.first].unlock();
                }
            });
        });

        tbb::parallel_for_each(nonGapMutations, [&](auto& u){
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++){
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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
        tbb::parallel_for_each(gapMutations, [&](auto& u){
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++){
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1]) || std::get<3>(u.second[i]) != std::get<3>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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
    } else if(ftype == PangenomeMAT2::FILE_TYPE::MSA){
        std::string newickString;
        secondFin >> newickString;
        root = createTreeFromNewickString(newickString);

        std::map< std::string, std::string > sequenceIdsToSequences;
        std::string line;
        std::string currentSequence, currentSequenceId;
        size_t lineLength = 0;
        std::string consensusSeq;
        while(getline(fin,line,'\n')){
            if(line.length() == 0){
                continue;
            }
            if(line[0] == '>'){
                if(currentSequence.length()){
                    if(lineLength == 0){
                        lineLength = currentSequence.length();
                    } else if(lineLength != currentSequence.length()){
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
        if(currentSequence.length()){
            if(lineLength != 0 && lineLength != currentSequence.length()){
                std::cerr << "Error: sequence lengths don't match!" << std::endl;
                exit(-1);
            } else {
                lineLength = currentSequence.length();
            }
            sequenceIdsToSequences[currentSequenceId] = currentSequence;
        }
        std::set< size_t > emptyPositions;

        for(size_t i = 0; i < lineLength; i++){
            bool nonGapFound = false;
            for(auto u: sequenceIdsToSequences){
                if(u.second[i] != '-'){
                    consensusSeq += u.second[i];
                    nonGapFound = true;
                    break;
                }
            }
            if(!nonGapFound){
                emptyPositions.insert(i);
            }
        }
        for(auto& u: sequenceIdsToSequences){
            std::string sequenceString;
            for(size_t i = 0; i < u.second.length(); i++){
                if(emptyPositions.find(i) == emptyPositions.end()){
                    sequenceString += u.second[i];
                }
            }
            u.second = sequenceString;
        }
        // std::cout << consensusSeq.length() << std::endl;
        blocks.emplace_back(0, consensusSeq);
        root->blockMutation.emplace_back(0, true);

        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes){
            nodeMutexes[u.first];
        }

        tbb::parallel_for((size_t)0, consensusSeq.length(), [&](size_t i){
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
            for(const auto& u: sequenceIdsToSequences){
                if(u.second[i] != '-'){
                    states[u.first] = (1 << getCodeFromNucleotide(u.second[i]));
                } else {
                    states[u.first] = 1;
                }
            }
            nucFitchForwardPass(root, states);
            nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(consensusSeq[i])));
            nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSeq[i])));
            for(auto mutation: mutations){
                nodeMutexes[mutation.first].lock();
                nonGapMutations[mutation.first].push_back(std::make_tuple(0, -1, i, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                // allNodes[mutation.first]->nucMutation.emplace_back( std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)) );
                nodeMutexes[mutation.first].unlock();
            }
        });

        tbb::parallel_for_each(nonGapMutations, [&](auto& u){
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++){
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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
        
    } else if(ftype == PangenomeMAT2::FILE_TYPE::MSA_OPTIMIZE){
        // std::string newickString;
        // secondFin >> newickString;
        // root = createTreeFromNewickString(newickString);

        // std::map< std::string, std::string > sequenceIdsToSequences;
        // std::string line;
        // std::string currentSequence, currentSequenceId;
        // size_t lineLength = 0;
        // std::string consensusSeq;
        // while(getline(fin,line,'\n')){
        //     if(line.length() == 0){
        //         continue;
        //     }
        //     if(line[0] == '>'){
        //         if(currentSequence.length()){
        //             if(lineLength == 0){
        //                 lineLength = currentSequence.length();
        //             } else if(lineLength != currentSequence.length()){
        //                 std::cerr << "Error: sequence lengths don't match! " << currentSequenceId << std::endl;
        //                 exit(-1);
        //             }
        //             sequenceIdsToSequences[currentSequenceId] = currentSequence;
        //         }
        //         std::vector< std::string > splitLine;
        //         stringSplit(line,' ',splitLine);
        //         currentSequenceId = splitLine[0].substr(1);
        //         currentSequence = "";
        //     } else {
        //         currentSequence += line;
        //     }
        // }
        // if(currentSequence.length()){
        //     if(lineLength != 0 && lineLength != currentSequence.length()){
        //         std::cerr << "Error: sequence lengths don't match!" << std::endl;
        //         exit(-1);
        //     } else {
        //         lineLength = currentSequence.length();
        //     }
        //     sequenceIdsToSequences[currentSequenceId] = currentSequence;
        // }
        // std::set< size_t > emptyPositions;

        // for(size_t i = 0; i < lineLength; i++){
        //     bool nonGapFound = false;
        //     for(auto u: sequenceIdsToSequences){
        //         if(u.second[i] != '-'){
        //             consensusSeq += u.second[i];
        //             nonGapFound = true;
        //             break;
        //         }
        //     }
        //     if(!nonGapFound){
        //         emptyPositions.insert(i);
        //     }
        // }
        // for(auto& u: sequenceIdsToSequences){
        //     std::string sequenceString;
        //     for(size_t i = 0; i < u.second.length(); i++){
        //         if(emptyPositions.find(i) == emptyPositions.end()){
        //             sequenceString += u.second[i];
        //         }
        //     }
        //     u.second = sequenceString;
        // }
    } else if(ftype == PangenomeMAT2::FILE_TYPE::FASTA){
        std::string newickString;
        secondFin >> newickString;
        root = createTreeFromNewickString(newickString);

        std::map< std::string, std::string > sequenceIdsToSequences;
        std::string line;
        std::string currentSequence, currentSequenceId;

        std::string consensusSeq;
        while(getline(fin,line,'\n')){
            if(line.length() == 0){
                continue;
            }
            if(line[0] == '>'){
                if(currentSequence.length()){
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
        if(currentSequence.length()){
            sequenceIdsToSequences[currentSequenceId] = currentSequence;
        }
        std::map< std::string, std::vector< std::string > > sequenceIdsToBlockSequences;
        size_t maxBlockId = 0;

        int BLOCKSIZE = 1000;

        for(auto u: sequenceIdsToSequences){
            for(size_t i = 0; i < u.second.size(); i+=BLOCKSIZE){
                sequenceIdsToBlockSequences[u.first].push_back(u.second.substr(i,BLOCKSIZE));
                maxBlockId = std::max(maxBlockId, i/BLOCKSIZE);
            }
        }

        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< size_t, std::string > > sequenceIdsToBlockMSA;
        tbb::concurrent_unordered_map< size_t, size_t > blockSizes;

        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< size_t, bool > > blockExists;

        // std::cout << "Reached here" << std::endl;
        // std::cout << sequenceIdsToBlockSequences.size() << std::endl;

        tbb::parallel_for((size_t)0,maxBlockId+1,[&](size_t i){
            std::vector< std::string > stringSequences;
            // int emptyCtr = 0;
            for(auto u: sequenceIdsToBlockSequences){
                if(i < u.second.size()){
                    stringSequences.push_back(u.second[i]);
                    blockExists[u.first][i] = true;
                } else {
                    blockExists[u.first][i] = false;
                }
            }
            auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -3, -10, -5);

            // if(stringSequences.size() != sequenceIdsToBlockSequences.size()){
            //     std::cout << "Grave error!!" << std::endl;
            // }

            spoa::Graph graph{};

            for (const auto& it : stringSequences) {
                auto alignment = alignment_engine->Align(it, graph);
                graph.AddAlignment(alignment, it);
            }

            auto consensus = graph.GenerateConsensus();
            auto msa = graph.GenerateMultipleSequenceAlignment();
            // if(msa.size() != stringSequences.size()){
            //     std::cout << "Grave error!!!!!! " << msa.size() <<" " << stringSequences.size() << " " << emptyCtr << " "<< msa.size()+emptyCtr << std::endl;
            // }

            int ctr = 0;
            blockSizes[i] = msa[0].size();
            // std::cout << i << " " << blockSizes[i] << std::endl;

            for(auto u: sequenceIdsToBlockSequences){
                if(i < u.second.size()){
                    sequenceIdsToBlockMSA[u.first][i] = msa[ctr];
                    ctr++;
                }
            }
        });

        // std::cout << "Reached here!" << std::endl;

        tbb::concurrent_unordered_map< size_t, std::string > consensusSequences;

        for(size_t i = 0; i <= maxBlockId; i++){
            std::string consensusSeq;
            std::set< size_t > emptyPositions;
            for(size_t j = 0; j < blockSizes[i]; j++){
                bool nonGapFound = false;
                for(auto& u: sequenceIdsToBlockMSA){
                    if(blockExists[u.first][i] && u.second[i][j] != '-'){
                        consensusSeq += u.second[i][j];
                        nonGapFound = true;
                        break;
                    }
                }
                if(!nonGapFound){
                    emptyPositions.insert(j);
                }
            }
            // if(emptyPositions.size()){
            //     std::cout << "Empty positions found!!!" << std::endl;
            // }
            consensusSequences[i] = consensusSeq;
            for(auto& u: sequenceIdsToBlockMSA){
                if(!blockExists[u.first][i]){
                    continue;
                }
                std::string currentSequence;
                for(size_t j = 0; j < blockSizes[i]; j++){
                    if(emptyPositions.find(j) == emptyPositions.end()){
                        currentSequence+=u.second[i][j];
                    }
                }
                u.second[i] = currentSequence;
            }
        }
        
        for(auto u: consensusSequences){
            blocks.emplace_back(u.first, u.second);
            // root->blockMutation.emplace_back(u.first, true);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string, bool > > globalBlockMutations;

        tbb::parallel_for((size_t)0, maxBlockId+1, [&](size_t i){
            // std::cout << i << std::endl;
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, bool > mutations;
            for(auto& u: blockExists){
                states[u.first] = u.second[i];
            }
            
            blockFitchForwardPass(root, states);
            blockFitchBackwardPass(root, states, 1);
            blockFitchAssignMutations(root, states, mutations, 0);
            globalBlockMutations[i] = mutations;

        });

        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes){
            nodeMutexes[u.first];
        }

        tbb::parallel_for_each(globalBlockMutations, [&](auto& pos){
            auto& mutations = pos.second;
            for(const auto& node: allNodes){
                if(mutations.find(node.first) != mutations.end()){
                    nodeMutexes[node.first].lock();
                    node.second->blockMutation.emplace_back(pos.first, mutations[node.first]);
                    nodeMutexes[node.first].unlock();
                }
            }
        });

        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        // std::unordered_map< std::string, std::mutex > nodeMutexes;

        // for(auto u: allNodes){
        //     nodeMutexes[u.first];
        // }


        tbb::parallel_for((size_t)0, maxBlockId+1, [&](size_t i){
            size_t blockSize = consensusSequences[i].size();
            tbb::parallel_for((size_t)0, blockSize, [&](size_t j){
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
                for(auto& u: sequenceIdsToBlockMSA){
                    if(blockExists[u.first][i]){
                        if(u.second[i][j] != '-'){
                            states[u.first] = (1 << getCodeFromNucleotide(u.second[i][j]));
                        } else {
                            states[u.first] = 1;
                        }
                    }
                }
                nucFitchForwardPass(root, states);
                nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(consensusSequences[i][j])));
                nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSequences[i][j])));
                for(auto mutation: mutations){
                    nodeMutexes[mutation.first].lock();
                    nonGapMutations[mutation.first].push_back(std::make_tuple(i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                    // allNodes[mutation.first]->nucMutation.emplace_back( std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)) );
                    nodeMutexes[mutation.first].unlock();
                }
            });
        });

        tbb::parallel_for_each(nonGapMutations, [&](auto& u){
            nodeMutexes[u.first].lock();
            std::sort(u.second.begin(), u.second.end());
            nodeMutexes[u.first].unlock();
            size_t currentStart = 0;
            for(size_t i = 1; i < u.second.size(); i++){
                if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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
}

void PangenomeMAT2::Tree::protoMATToTree(const MATNew::tree& mainTree){
    // Create tree
    root = createTreeFromNewickString(mainTree.newick());
    invertTree(root);

    std::map< int64_t, std::vector< uint32_t > > blockIdToConsensusSeq;

    for(int i = 0; i < mainTree.consensusseqmap_size(); i++){
        std::vector< uint32_t > seq;
        for(int j = 0; j < mainTree.consensusseqmap(i).consensusseq_size(); j++){
            seq.push_back(mainTree.consensusseqmap(i).consensusseq(j));
        }
        for(int j = 0; j < mainTree.consensusseqmap(i).blockid_size(); j++){
            blockIdToConsensusSeq[mainTree.consensusseqmap(i).blockid(j)] = seq;
        }
    }

    std::vector< MATNew::node > storedNodes;
    for(int i = 0; i < mainTree.nodes_size(); i++){
        storedNodes.push_back(mainTree.nodes(i));
    }

    size_t initialIndex = 0;

    assignMutationsToNodes(root, initialIndex, storedNodes);

    // Block sequence
    for(int i = 0; i < mainTree.blocks_size(); i++){
        blocks.emplace_back(mainTree.blocks(i), blockIdToConsensusSeq[mainTree.blocks(i).blockid()]);
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

PangenomeMAT2::Tree::Tree(const MATNew::tree& mainTree){
    protoMATToTree(mainTree);
}

PangenomeMAT2::Tree::Tree(std::ifstream& fin, FILE_TYPE ftype){

    if(ftype == PangenomeMAT2::FILE_TYPE::PANMAT){
        MATNew::tree mainTree;

        if(!mainTree.ParseFromIstream(&fin)){
            throw std::invalid_argument("Could not read tree from input file.");
        }

        protoMATToTree(mainTree);

        // // Rerooting debug
        // for(auto u: allNodes){
        //     if(u.second->children.size() == 0){
        //         std::cout << "Rerooting to " << u.first << std::endl;
        //         reroot(u.first);
        //         return;
        //     }
        // }
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

                    if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x'){
                        line += sequence[i].second[j][k].first;
                    } else if(aligned){
                        line += '-';
                    }
                    if(line.length() == lineSize){
                        fout << line << '\n';
                        line = "";
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

                    line += '-';
                    if(line.length() == lineSize){
                        fout << line << '\n';
                        line = "";
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

                if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x'){
                    line += sequence[i].first[j].first;
                } else if(aligned){
                    line += '-';
                }
                if(line.length() == lineSize){
                    fout << line << '\n';
                    line = "";
                }
            }
        } else if(aligned) {
            for(size_t j = 0; j < sequence[i].first.size(); j++){
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++){
                    line+='-';
                    if(line.length() == lineSize){
                        fout << line << '\n';
                        line = "";
                    }
                }
                line+='-';
                if(line.length() == lineSize){
                    fout << line << '\n';
                    line = "";
                }
            }
        }

    }


    if(line.length()){
        fout << line << '\n';
        line = "";
    }

}

char PangenomeMAT2::getCodeFromNucleotide(char nuc){
    switch(nuc){
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

        // if(primaryBlockId == 4){
        //     std::cout << i << " " << blocks[i].consensusSeq.size() << std::endl;
        // }

        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                // if(primaryBlockId == 4){
                //     std::cout << nucCode << std::endl;
                // }

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
            // if(primaryBId > sequence.size()){
            //     std::cout << "Exceeded 1 " << primaryBId << " " << sequence.size() << std::endl;
            // }
            // if(pos >= sequence[primaryBId].first.size()){
            //     std::cout << "Exceeded 2 " << pos << " " << sequence[primaryBId].first.size() << " " << primaryBId << std::endl;
            // }
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
    std::map< std::pair< int32_t, int32_t >, std::pair< std::vector< MATNew::nucMut >, int > > blockToMutations;

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
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].first.push_back(nm);
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].second = 2;
    }

    for(size_t i = 0; i < root->blockMutation.size(); i++){
        const PangenomeMAT2::BlockMut& mutation = root->blockMutation[i];
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].second = mutation.blockMutInfo;
    }

    for(auto u: blockToMutations){
        MATNew::mutation mutation;
        mutation.set_blockmutexist((u.second.second != 2));
        mutation.set_blockmutinfo(u.second.second);

        int32_t primaryBlockId = u.first.first;
        int32_t secondaryBlockId = u.first.second;
        if(secondaryBlockId != -1){
            mutation.set_blockid(((int64_t)primaryBlockId << 32) + secondaryBlockId);
            mutation.set_blockgapexist(true);
        } else {
            mutation.set_blockid(((int64_t)primaryBlockId << 32));
            mutation.set_blockgapexist(false);
        }
        for(auto v: u.second.first){
            mutation.add_nucmutation();
            *mutation.mutable_nucmutation(mutation.nucmutation_size() - 1) = v;
        }
        n.add_mutations();
        *n.mutable_mutations(n.mutations_size() - 1) = mutation;
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

    std::map< std::vector< uint32_t >, std::vector< int64_t > > consensusSeqToBlockIds;

    for(auto block: blocks){
        MATNew::block b;
        int64_t blockId;
        if(block.secondaryBlockId != -1){
            blockId = ((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId;
            b.set_blockid(blockId);
            b.set_blockgapexist(true);
        } else {
            blockId = ((int64_t)block.primaryBlockId << 32);
            b.set_blockid(blockId);
            b.set_blockgapexist(false);
        }
        b.set_chromosomename(block.chromosomeName);
        consensusSeqToBlockIds[block.consensusSeq].push_back(blockId);
        // for(auto n: block.consensusSeq){
        //     b.add_consensusseq(n);
        // }
        treeToWrite.add_blocks();
        *treeToWrite.mutable_blocks( treeToWrite.blocks_size() - 1 ) = b;
    }

    for(auto u: consensusSeqToBlockIds){
        MATNew::consensusSeqToBlockIds c;
        for(auto v: u.first){
            c.add_consensusseq(v);
        }
        for(auto v: u.second){
            c.add_blockid(v);
        }
        treeToWrite.add_consensusseqmap();
        *treeToWrite.mutable_consensusseqmap( treeToWrite.consensusseqmap_size() - 1 ) = c;
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
		std::cerr << "Failed to write to output file." << std::endl;
    }

}

void PangenomeMAT2::Tree::getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists, std::string reference){
    Node* referenceNode = nullptr;

    for(auto u: allNodes){
        if(u.second->children.size() == 0 && u.first == reference){
            referenceNode = u.second;
            break;
        }
    }

    if(referenceNode == nullptr){
        std::cerr << "Error: Reference sequence with matching name not found!" << std::endl;
        return;
    }

    std::vector< PangenomeMAT2::Node* > path;
    Node* it = referenceNode;

    while(it != root){
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(root);

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    sequence.resize(blocks.size() + 1);
    blockExists.resize(blocks.size() + 1, {false, {}});

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

    // First we check if there are any nucleotide mutations. If there are no nuc mutations, we can simply construct a GFA of blocks.
    bool nucMutationFlag = false;
    for(auto u: allNodes){
        if(u.second->nucMutation.size() != 0){
            nucMutationFlag = true;
        }
    }

    if(!nucMutationFlag){
        // get nodes
        std::map<std::pair<int32_t, int32_t>, std::string> nodes;
        for(auto block: blocks){
            int64_t primaryBlockId = block.primaryBlockId;
            int64_t secondaryBlockId = block.secondaryBlockId;
            std::string sequenceString;
            for(auto u: block.consensusSeq){
                for(size_t k = 0; k < 8; k++){
                    const int nucCode = (((u) >> (4*(7 - k))) & 15);
                    if(nucCode == 0){
                        break;
                    }
                    const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                    sequenceString += nucleotide;
                }
            }
            nodes[std::make_pair(primaryBlockId, secondaryBlockId)] = sequenceString;
        }


        // block presense map
        std::vector< std::pair< bool, std::vector< bool > > > blockExistsGlobal(blocks.size() + 1, {false, {}});
        // Assigning block gaps
        for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
            blockExistsGlobal[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        }
        
        tbb::concurrent_unordered_set< std::pair< std::pair<int32_t, int32_t>, std::pair<int32_t, int32_t> > > edges;
        tbb::concurrent_unordered_map< std::string, std::vector< std::pair<int32_t, int32_t> > > paths;

        // get all paths
        tbb::parallel_for_each(allNodes, [&](auto u){
            if(u.second->children.size()){
                return;
            }

            auto blockExists = blockExistsGlobal;

            std::vector< PangenomeMAT2::Node* > path;

            Node* it = u.second;
            while(it != root){
                path.push_back(it);
                it = it->parent;
            }
            path.push_back(root);
            std::reverse(path.begin(), path.end());
            for(auto node: path){
                for(auto mutation: node->blockMutation){
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
            std::vector< std::pair< int32_t, int32_t > > currentPath;
            for(size_t i = 0; i < blockExists.size(); i++){
                if(blockExists[i].first){
                    currentPath.push_back(std::make_pair(i, -1));
                }
                for(size_t j = 0; j < blockExists[i].second.size(); j++){
                    if(blockExists[i].second[j]){
                        currentPath.push_back(std::make_pair(i, j));
                    }
                }
            }
            paths[u.second->identifier] = currentPath;
            for(size_t i = 1; i < currentPath.size(); i++){
                edges.insert(std::make_pair(currentPath[i-1], currentPath[i]));
            }
        });
        std::map< std::pair< int32_t, int32_t >, uint64_t > nodeIds;
        uint64_t ctr = 0;
        for(auto u: nodes){
            nodeIds[u.first] = ctr;
            ctr++;
        }
        for(auto u: nodes){
            fout << "S\t" << nodeIds[u.first] << "\t" << u.second << "\n";
        }
        for(auto u: edges){
            fout << "L\t" << nodeIds[u.first] << "\t+\t" << nodeIds[u.second] << "\t+\t0M\n";
        }
        for(auto u: paths){
            fout << "P\t" << u.first << "\t";
            for(size_t i = 0; i < u.second.size(); i++){
                fout << nodeIds[u.second[i]] << "+";
                if(i != u.second.size() - 1){
                    fout << ",";
                }
            }
            fout << "\t*\n";
        }
    } else {
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
                if(G[dest].size() != GT[dest].size()){
                    break;
                }
                for(size_t j = 0; j < GT[dest].size(); j++){
                    if(G[dest][j].first != GT[dest][j].first){
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

        // convert node IDs to consecutive node IDs after path compression
        int currentID = 0;
        std::map< size_t, size_t > oldToNew;
        for(auto u: finalNodes){
            oldToNew[u.first] = currentID;
            currentID++;
        }

        fout << "H\tVN:Z:1.1\n";
        for(auto u: finalNodes){
            fout << "S\t" << oldToNew[u.first] << "\t" << u.second << "\n";
        }

        for(auto u: edges){
            fout << "L\t" << oldToNew[u.first] << "\t+\t" << oldToNew[u.second] << "\t+\t0M\n";
        }

        for(auto u: paths){
            fout << "P\t" << u.first << "\t";
            for(size_t i = 0; i < u.second.size(); i++){
                fout << oldToNew[u.second[i]] << "+";
                if(i != u.second.size() - 1){
                    fout << ",";
                }
            }
            fout << "\t*\n";
        }
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

std::tuple< int, int, int, int > PangenomeMAT2::Tree::globalCoordinateToBlockCoordinate(int64_t globalCoordinate, const sequence_t& sequence, const blockExists_t& blockExists){
    int ctr = 0;
    for(size_t i = 0; i < blockExists.size(); i++){
        for(size_t j = 0; j < blockExists[i].second.size(); j++){
            if(!blockExists[i].second[j]){
                continue;
            } else {
                for(size_t k = 0; k < sequence[i].second[j].size(); k++){
                    for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++){
                        if(sequence[i].second[j][k].second[w] != '-' && sequence[i].second[j][k].second[w] != 'x'){
                            // useful character found
                            if(ctr == globalCoordinate){
                                return std::make_tuple(i,j,k,w);
                            }
                            ctr++;
                        }
                    }
                    if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x'){
                        // useful character found
                        if(ctr == globalCoordinate){
                            return std::make_tuple(i,j,k,-1);
                        }
                        ctr++;
                    }
                }
            }
        }
        if(!blockExists[i].first){
            continue;
        }
        for(size_t k = 0; k < sequence[i].first.size(); k++){
            for(size_t w = 0; w < sequence[i].first[k].second.size(); w++){
                if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x'){
                    if(ctr == globalCoordinate){
                        return std::make_tuple(i, -1, k, w);
                    }
                    ctr++;
                }
            }
            if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x'){
                if(ctr == globalCoordinate){
                    return std::make_tuple(i, -1, k, -1);
                }
                ctr++;
            }
        }
    }
    return std::make_tuple(-1,-1,-1,-1);
}

void PangenomeMAT2::Tree::adjustLevels(Node* node){
    if(node->parent == nullptr){
        node->level = 1;
    } else {
        node->level = node->parent->level + 1;
    }
    for(auto u: node->children){
        adjustLevels(u);
    }
}

PangenomeMAT2::Node* PangenomeMAT2::Tree::transformHelper(Node* node){
    if(node == root){
        if(node->children.size() > 1){
            node->branchLength = 0;
            return node;
        } else {
            node = node->children[0];
            node->branchLength = 0;
            delete root;
            allNodes.erase(root->identifier);
            return node;
        }
    }
    Node* par = node->parent;

    // erase node from parent's children
    for(size_t i = 0; i < par->children.size(); i++){
        if(par->children[i]->identifier == node->identifier){
            par->children.erase((par->children).begin()+i);
            break;
        }
    }

    node->parent = nullptr;
    size_t oldBranchLen = node->branchLength;
    node->branchLength = 0;

    // make transformed parent a new child of the node
    Node* newChild = transformHelper(par);
    node->children.push_back(newChild);
    newChild->parent = node;
    newChild->branchLength = oldBranchLen;

    return node;
}

void PangenomeMAT2::Tree::transform(Node* node){
    Node* par = node->parent;
    if(par == nullptr){
        // already root
        return;
    }
    if(par == root){
        // Parent already root. The root will contain the same sequence as the node
        node->branchLength = 0;
        return;
    }

    // remove node from parent's children
    for(size_t i = 0; i < par->children.size(); i++){
        if(par->children[i]->identifier == node->identifier){
            par->children.erase((par->children).begin()+i);
            break;
        }
    }
    node->parent = nullptr;

    size_t oldBranchLen = node->branchLength;

    Node* newRoot = new Node(newInternalNodeId(), 0);
    newRoot->children.push_back(node);
    node->parent = newRoot;
    node->level = 2;
    node->branchLength = 0;
    
    // transform node's parent into node's sibling
    Node* newSibling = transformHelper(par);
    newRoot->children.push_back(newSibling);
    newSibling->parent = newRoot;
    newSibling->branchLength = oldBranchLen;

    root = newRoot;
    allNodes[root->identifier] = root;
    adjustLevels(root);

}

void PangenomeMAT2::Tree::reroot(std::string sequenceName){
    if(allNodes.find(sequenceName) == allNodes.end()){
        std::cerr << "Sequence with name " << sequenceName << " not found!" << std::endl;
        return;
    }
    PangenomeMAT2::Node* newRoot = allNodes[sequenceName];
    if(newRoot->children.size()){
        std::cerr << "Node with id " << sequenceName << " is not a tip!" << std::endl;
        return;
    }
    sequence_t sequence;
    blockExists_t blockExists;
    getSequenceFromReference(sequence, blockExists, sequenceName);

    std::unordered_map<std::string, sequence_t> nodeIdToSequence;
    std::unordered_map<std::string, blockExists_t> nodeIdToBlockExists;

    for(const auto& u: allNodes){
        if(u.second->children.size() == 0){
            sequence_t tempSequence;
            blockExists_t tempBlockExists;
            getSequenceFromReference(tempSequence, tempBlockExists, u.first);
            nodeIdToSequence[u.first] = tempSequence;
            nodeIdToBlockExists[u.first] = tempBlockExists;
        }
    }


    // Transform tree topology
    transform(newRoot);

    std::cout << "Transformation complete!" << std::endl;

    // clear previous block mutations
    for(auto u: allNodes){
        u.second->blockMutation.clear();
    }

    std::unordered_map< std::string, std::mutex > nodeMutexes;

    for(auto u: allNodes){
        nodeMutexes[u.first];
    }

    // make new block mutations
    tbb::parallel_for((size_t)0, blockExists.size(), [&](size_t i){
        tbb::parallel_for((size_t)0, blockExists[i].second.size(), [&](size_t j){
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, bool > mutations;
            // block gaps
            for(const auto& u: nodeIdToBlockExists){
                states[u.first] = u.second[i].second[j];
            }

            blockFitchForwardPass(root, states);
            blockFitchBackwardPass(root, states, 1, blockExists[i].second[j]);
            blockFitchAssignMutations(root, states, mutations, 0);
            for(auto u: mutations){
                nodeMutexes[u.first].lock();
                allNodes[u.first]->blockMutation.emplace_back(i, u.second, j);
                nodeMutexes[u.first].unlock();
            }
        });
        std::unordered_map< std::string, int > states;
        std::unordered_map< std::string, bool > mutations;
        // main block
        for(const auto& u: nodeIdToBlockExists){
            states[u.first] = u.second[i].first;
        }

        blockFitchForwardPass(root, states);
        blockFitchBackwardPass(root, states, 1, blockExists[i].first);
        blockFitchAssignMutations(root, states, mutations, 0);
        for(auto u: mutations){
            nodeMutexes[u.first].lock();
            allNodes[u.first]->blockMutation.emplace_back(i, u.second, -1);
            nodeMutexes[u.first].unlock();
        }
    });

    std::cout << "Block Mutations Ready!" << std::endl;

    // clear previous nuc mutations
    for(auto u: allNodes){
        u.second->nucMutation.clear();
    }

    tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
    tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > gapMutations;

    tbb::parallel_for((size_t)0, sequence.size(), [&](size_t i){
    //     tbb::parallel_for((size_t)0, sequence[i].second.size(), [&](size_t j){
    //         // gap blocks
    //         std::vector< char > consensusSeq;
    //         bool endFlag = false;
    //         for(size_t t1 = 0; t1 < blocks.size(); t1++){
    //             if(blocks[t1].primaryBlockId == (int)i && blocks[t1].secondaryBlockId == (int)j){
    //                 for(size_t t2 = 0; t2 < blocks[t1].consensusSeq.size(); t2++){
    //                     for(size_t t3 = 0; t3 < 8; t3++){
    //                         const int nucCode = (((blocks[t1].consensusSeq[t2]) >> (4*(7 - t3))) & 15);
    //                         if(nucCode == 0){
    //                             endFlag = true;
    //                             break;
    //                         }
    //                         const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                            
    //                         consensusSeq.push_back(nucleotide);
    //                     }
    //                     if(endFlag){
    //                         break;
    //                     }
    //                 }
    //                 endFlag = true;
    //                 break;
    //             }
    //         }

    //         if(!endFlag){
    //             std::cerr << "FATAL: Block with id " << i << " " << j << " not found!" << std::endl;
    //             exit(-1);
    //         }

    //         consensusSeq.push_back('-');
    //         if(consensusSeq.size() != sequence[i].second[j].size()){
    //             std::cerr << "FATAL: consenusSeq length doesn't match with sequence length" << std::endl;
    //             exit(-1);
    //         }

    //         std::cout << "BID " << i << " " << j << " started" << std::endl;

    //         tbb::parallel_for((size_t)0, sequence[i].second[j].size(), [&](size_t k){
    //             tbb::parallel_for((size_t)0, sequence[i].second[j][k].second.size(), [&](size_t w){
    //                 // gap nuc
    //                 std::unordered_map< std::string, int > states;
    //                 std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
    //                 for(auto u: nodeIdToSequence){
    //                     if(u.second[i].second[j][k].second[w] != '-' && u.second[i].second[j][k].second[w] != 'x'){
    //                         states[u.first] = (1 << getCodeFromNucleotide(u.second[i].second[j][k].second[w]));
    //                     }  else {
    //                         states[u.first] = 1;
    //                     }
    //                 }
    //                 char nucleotideCode = 1;
    //                 if(sequence[i].second[j][k].second[w] != '-' && sequence[i].second[j][k].second[w] != 'x'){
    //                     nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].second[j][k].second[w]));
    //                 }
    //                 nucFitchForwardPass(root, states);
    //                 nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
    //                 nucFitchAssignMutations(root, states, mutations, 1);
    //                 for(auto mutation: mutations){
    //                     nodeMutexes[mutation.first].lock();
    //                     gapMutations[mutation.first].push_back(std::make_tuple((int)i, j, k, w, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
    //                     nodeMutexes[mutation.first].unlock();
    //                 }
    //             });
    //             // main nuc
    //             std::unordered_map< std::string, int > states;
    //             std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;

    //             for(auto u: nodeIdToSequence){
    //                 if(u.second[i].second[j][k].first != '-' && u.second[i].second[j][k].first != 'x'){
    //                     states[u.first] = (1 << getCodeFromNucleotide(u.second[i].second[j][k].first));
    //                 }  else {
    //                     states[u.first] = 1;
    //                 }
    //             }
    //             char nucleotideCode = 1;
    //             if(sequence[i].second[j][k].first != '-' && sequence[i].second[j][k].first != 'x'){
    //                 nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].second[j][k].first));
    //             }
    //             nucFitchForwardPass(root, states);
    //             nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
    //             nucFitchAssignMutations(root, states, mutations, (consensusSeq[k] == '-'? 1: (1 << getCodeFromNucleotide(consensusSeq[k]))));
    //             for(auto mutation: mutations){
    //                 nodeMutexes[mutation.first].lock();
    //                 nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, j, k, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
    //                 nodeMutexes[mutation.first].unlock();
    //             }
    //         });

    //         std::cout << "BID " << i << " " << j << " ended" << std::endl;
    //     });

        // main block
        std::vector< char > consensusSeq;
        bool endFlag = false;
        for(size_t t1 = 0; t1 < blocks.size(); t1++){
            if(blocks[t1].primaryBlockId == (int)i && blocks[t1].secondaryBlockId == -1){
                for(size_t t2 = 0; t2 < blocks[t1].consensusSeq.size(); t2++){
                    for(size_t t3 = 0; t3 < 8; t3++){
                        const int nucCode = (((blocks[t1].consensusSeq[t2]) >> (4*(7 - t3))) & 15);
                        if(nucCode == 0){
                            endFlag = true;
                            break;
                        }
                        const char nucleotide = PangenomeMAT2::getNucleotideFromCode(nucCode);
                        
                        consensusSeq.push_back(nucleotide);
                    }
                    if(endFlag){
                        break;
                    }
                }
                endFlag = true;
                break;
            }
        }

        if(!endFlag){
            std::cerr << "FATAL: Block with id " << i << " " << -1 << " not found!" << std::endl;
            exit(-1);
        }

        consensusSeq.push_back('-');
        if(consensusSeq.size() != sequence[i].first.size()){
            std::cerr << "FATAL: consenusSeq length doesn't match with sequence length" << std::endl;
            exit(-1);
        }

        std::cout << "BID " << i << " " << -1 << " started" << std::endl;

        tbb::parallel_for((size_t)0, sequence[i].first.size(), [&](size_t k){

            tbb::parallel_for((size_t)0, sequence[i].first[k].second.size(), [&](size_t w){
                // gap nuc
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;
                for(const auto& u: nodeIdToSequence){
                    if(u.second[i].first[k].second[w] != '-' && u.second[i].first[k].second[w] != 'x'){
                        states[u.first] = (1 << getCodeFromNucleotide(u.second[i].first[k].second[w]));
                    }  else {
                        states[u.first] = 1;
                    }
                }
                int nucleotideCode = 1;
                if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x'){
                    nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].first[k].second[w]));
                }

                // nucFitchForwardPass(root, states);
                // nucFitchBackwardPass(root, states, 1);
                // nucFitchAssignMutations(root, states, mutations, 1);

                nucFitchForwardPass(root, states);
                nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
                nucFitchAssignMutations(root, states, mutations, 1);
                
                for(auto mutation: mutations){
                    nodeMutexes[mutation.first].lock();
                    gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, k, w, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                    nodeMutexes[mutation.first].unlock();
                }
            });
            // main nuc
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > > mutations;

            for(const auto& u: nodeIdToSequence){
                if(u.second[i].first[k].first != '-' && u.second[i].first[k].first != 'x'){
                    states[u.first] = (1 << getCodeFromNucleotide(u.second[i].first[k].first));
                }  else {
                    states[u.first] = 1;
                }
            }
            int nucleotideCode = 1;
            if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x'){
                nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].first[k].first));
            }
            
            // nucFitchForwardPass(root, states);
            // nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(consensusSeq[k])));
            // nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSeq[k])));

            nucFitchForwardPass(root, states);
            nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
            nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSeq[k])));

            for(auto mutation: mutations){
                nodeMutexes[mutation.first].lock();
                nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, k, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                nodeMutexes[mutation.first].unlock();
            }
        });

        std::cout << "BID " << i << " " << -1 << " ended" << std::endl;
    });

    tbb::parallel_for_each(nonGapMutations, [&](auto& u){
        nodeMutexes[u.first].lock();
        std::sort(u.second.begin(), u.second.end());
        nodeMutexes[u.first].unlock();
        size_t currentStart = 0;
        for(size_t i = 1; i < u.second.size(); i++){
            if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<1>(u.second[i]) != std::get<1>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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
    tbb::parallel_for_each(gapMutations, [&](auto& u){
        nodeMutexes[u.first].lock();
        std::sort(u.second.begin(), u.second.end());
        nodeMutexes[u.first].unlock();
        size_t currentStart = 0;
        for(size_t i = 1; i < u.second.size(); i++){
            if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<1>(u.second[i]) != std::get<1>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1]) || std::get<3>(u.second[i]) != std::get<3>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])){
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

bool PangenomeMAT2::GFAGraph::pathExists(size_t nId1, size_t nId2, std::vector< bool >& visited){
    if(nId1 == nId2){
        return true;
    }
    visited[nId1] = true;

    for(auto u: adj[nId1]){
        if(!visited[u] && pathExists(u, nId2, visited)){
            return true;
        }
    }
    return false;
}

PangenomeMAT2::GFAGraph::GFAGraph(const std::vector< std::string >& pathNames, const std::vector< std::vector< std::string > >& sequences, std::map< std::string, std::string >& nodes){
    pathIds = pathNames;

    // New integer Node ID to old string Node ID
    std::vector< std::string > nodeIdToString;

    // Old string Node ID to new integer Node ID
    std::unordered_map< std::string, size_t > stringToNodeId;

    // Auto increment ID to assign to nodes
    numNodes = 0;

    for(auto sequence: sequences){
        for(auto blockId: sequence){
            if(stringToNodeId.find(blockId) == stringToNodeId.end()){
                nodeIdToString.push_back(blockId);
                stringToNodeId[blockId] = numNodes;
                numNodes++;
            }
        }
    }

    adj.resize(numNodes);
    std::set< std::pair< size_t, size_t > > edges;

    for(auto sequence: sequences){
        for(size_t i = 1; i < sequence.size(); i++){
            size_t nId1 = stringToNodeId[sequence[i-1]];
            size_t nId2 = stringToNodeId[sequence[i]];
            if(edges.find( std::make_pair(nId1, nId2) ) == edges.end()){
                edges.insert(std::make_pair(nId1, nId2));
                adj[nId1].push_back(nId2);
            }
        }
    }

    intSequences.resize(sequences.size());
    for(size_t i = 0; i < sequences.size(); i++){
        const std::vector< std::string >& sequence = sequences[i];
        for(size_t j = 0; j < sequence.size(); j++){
            intSequences[i].push_back(stringToNodeId[sequence[j]]);
        }
    }

    for(size_t i = 0; i < nodeIdToString.size(); i++){
        intNodeToSequence.push_back(nodes[nodeIdToString[i]]);
    }

    if(!checkForCycles()){
        return;
    }

    std::cout << "Cycles found in GFA. Removing them..." << std::endl;

    // Remove cycles
    std::unordered_map< size_t, size_t > newToOld;
    for(size_t i = 0; i < numNodes; i++){
        newToOld[i] = i;
    }

    adj.clear();
    adj.resize(numNodes);

    for(auto& sequence: intSequences){
        if(sequence.size() == 0){
            continue;
        }
        size_t currentNode = sequence[0];
        for(size_t i = 1; i < sequence.size(); i++){
            size_t nextNode = sequence[i];
            bool neighbourFound = false;
            for(auto neighbour: adj[currentNode]){
                if(newToOld[neighbour] == nextNode){
                    neighbourFound = true;
                    sequence[i] = neighbour;
                    nextNode = neighbour;
                    break;
                }
            }

            if(neighbourFound){
                currentNode = nextNode;
                continue;
            }

            // if(!neighbourFound){
            //     adj[currentNode].push_back(nextNode);
            // }

            std::vector< bool > visited(numNodes, false);
            if(!pathExists(nextNode, currentNode, visited)){
                adj[currentNode].push_back(nextNode);
                currentNode = nextNode;
            } else {
                // size_t j = 0;
                // for(; j < adj[currentNode].size(); j++){
                //     if(adj[currentNode][j] == nextNode){
                //         break;
                //     }
                // }
                // if(j < adj[currentNode].size()){
                //     adj[currentNode].erase(adj[currentNode].begin() + j);
                // } else {
                //     std::cerr << "Unexpectedly didn't find nextNode in currentNode's adjacency list!" << std::endl;
                //     exit(-1);
                // }
                adj[currentNode].push_back(numNodes);
                adj.push_back({});
                newToOld[numNodes] = newToOld[nextNode];
                nodeIdToString.push_back(nodeIdToString[newToOld[numNodes]]);
                intNodeToSequence.push_back(nodes[nodeIdToString[numNodes]]);
                sequence[i] = numNodes;
                currentNode = numNodes;
                numNodes++;
            }
        }
    }

    if(checkForCycles()){
        std::cerr << "FATAL: Cycle removal unsuccessful! Please report this." << std::endl;
        exit(-1);
    }

    // Assign sequences to integer node IDs


    // std::cout << "Cycle Check New: " << checkForCycles() << "\n";

    // std::cout << "Sequences:\n";
    // for(auto sequence: intSequences){
    //     for(auto u: sequence){
    //         std::cout << u << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // for(int i = 0; i < autoIncrId; i++){
    //     std::cout << i << ": ";
    //     for(auto u: adj[i]){
    //         std::cout << u << " ";
    //     }
    //     std::cout << std::endl;
    // }

}

bool PangenomeMAT2::GFAGraph::checkForCyclesHelper(size_t nodeId, std::vector< int >& color){
    color[nodeId] = 1;
    for(auto u: adj[nodeId]){
        if(color[u] == 1){
            return true;
        }
        if(color[u] == 0){
            if(checkForCyclesHelper(u, color)){
                return true;
            }
        }
    }
    color[nodeId] = 2;
    return false;
}

bool PangenomeMAT2::GFAGraph::checkForCycles(){

    std::vector< int > color(numNodes, 0);
    for(size_t i = 0; i < numNodes; i++){
        if(color[i] == 0){
            if(checkForCyclesHelper(i, color)){
                return true;
            }
        }
    }
    return false;
}

std::vector< std::vector< int64_t > > PangenomeMAT2::GFAGraph::getAlignedSequences(const std::vector< size_t >& topoArray){
    std::vector< std::vector< int64_t > > newSequences;

    for(auto sequence: intSequences){
        size_t p1 = 0, p2 = 0;
        std::vector< int64_t > newSequence;
        while(p1 < topoArray.size() && p2 < sequence.size()){
            if(topoArray[p1] == sequence[p2]){
                newSequence.push_back(sequence[p2]);
                p2++;
            } else {
                newSequence.push_back(-1);
            }
            p1++;
        }
        // Adding gaps towards the end
        while(newSequence.size() < topoArray.size()){
            newSequence.push_back(-1);
        }
        newSequences.push_back(newSequence);
    }
    return newSequences;
}

void PangenomeMAT2::GFAGraph::topologicalSortHelper(size_t nodeId, std::vector< size_t >& topoArray, std::vector< bool >& visited){
    visited[nodeId] = true;
    for(auto u: adj[nodeId]){
        if(!visited[u]){
            topologicalSortHelper(u, topoArray, visited);
        }
    }
    topoArray.push_back(nodeId);
}

std::vector< size_t > PangenomeMAT2::GFAGraph::getTopologicalSort(){
    std::vector< size_t > topoArray;
    std::vector< bool > visited(numNodes, false);

    for(size_t i = 0; i < numNodes; i++){
        if(!visited[i]){
            topologicalSortHelper(i, topoArray, visited);
        }
    }
    std::reverse(topoArray.begin(), topoArray.end());

    return topoArray;
}

PangenomeMAT2::Pangraph::Pangraph(Json::Value& pangraphData){
    // load paths
    for(size_t i = 0; i < pangraphData["paths"].size(); i++){
        Json::Value path = pangraphData["paths"][(int)i];
        for(size_t j = 0; j < path["blocks"].size(); j++){
            paths[path["name"].asString()].push_back(path["blocks"][(int)j]["id"].asString());
        }
    }

    // load blocks
    for(size_t i = 0; i < pangraphData["blocks"].size(); i++){
        std::string blockId = pangraphData["blocks"][(int)i]["id"].asString();
        std::string sequence = pangraphData["blocks"][(int)i]["sequence"].asString();
        std::transform(sequence.begin(), sequence.end(),sequence.begin(), ::toupper);
        stringIdToConsensusSeq[blockId] = sequence;
        // stringIdToGaps[blockId].push_back(std::make_pair);
        std::vector< std::string > gapMemberNames = pangraphData["blocks"][(int)i]["gaps"].getMemberNames();
        for(auto member: gapMemberNames){
            stringIdToGaps[blockId].push_back( std::make_pair( std::stoi(member), pangraphData["blocks"][(int)i]["gaps"][member].asInt() ) );
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["mutate"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["number"].asInt();
            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["mutate"][(int)j][1].size(); k++){
                std::string mutationString = pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                substitutions[blockId][seqName][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt(), mutationString) );
                // std::cout << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt() << " " << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString() << std::endl;
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["insert"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["number"].asInt();
            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["insert"][(int)j][1].size(); k++){
                std::string mutationString = pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                insertions[blockId][seqName][number].push_back( std::make_tuple( pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][0].asInt(), pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][1].asInt(), mutationString ) );
                // std::cout << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt() << " " << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString() << std::endl;
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["delete"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["number"].asInt();
            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["delete"][(int)j][1].size(); k++){
                deletions[blockId][seqName][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][0].asInt(), pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][1].asInt() ) );
                // std::cout << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt() << " " << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString() << std::endl;
            }
        }
        
    }

    // Auto increment ID to assign to nodes
    numNodes = 0;

    // Old string Node ID to new integer Node ID
    std::unordered_map< std::string, size_t > stringToNodeId;

    for(const auto& p: paths){
        for(const auto& block: p.second){
            if(stringToNodeId.find(block) == stringToNodeId.end()){
                stringToNodeId[block] = numNodes;
                intIdToStringId[numNodes] = block;
                numNodes++;
            }
        }
    }

    adj.resize(numNodes);
    std::set< std::pair< size_t, size_t > > edges;

    for(const auto& p: paths){
        for(size_t i = 1; i < p.second.size(); i++){
            size_t nId1 = stringToNodeId[p.second[i-1]];
            size_t nId2 = stringToNodeId[p.second[i]];
            if(edges.find( std::make_pair(nId1, nId2) ) == edges.end()){
                edges.insert(std::make_pair(nId1, nId2));
                adj[nId1].push_back(nId2);
            }
        }
    }

    for(const auto& p: paths){
        for(size_t i = 0; i < p.second.size(); i++){
            intSequences[p.first].push_back(stringToNodeId[p.second[i]]);
        }
    }

    if(!checkForCycles()){
        return;
    }

    std::cout << "Cycles found in Pangraph. Removing them..." << std::endl;

    // Remove cycles
    // std::unordered_map< size_t, size_t > newToOld;
    // for(size_t i = 0; i < numNodes; i++){
    //     newToOld[i] = i;
    // }
    adj.clear();
    adj.resize(numNodes);

    for(auto& sequence: intSequences){
        if(sequence.second.size() == 0){
            continue;
        }
        size_t currentNode = sequence.second[0];
        for(size_t i = 1; i < sequence.second.size(); i++){
            size_t nextNode = sequence.second[i];
            bool neighbourFound = false;
            for(auto neighbour: adj[currentNode]){
                if(intIdToStringId[neighbour] == intIdToStringId[nextNode]){
                    neighbourFound = true;
                    nextNode = neighbour;
                    sequence.second[i] = neighbour;
                    break;
                }
            }
            if(neighbourFound){
                currentNode = nextNode;
                continue;
            }
            std::vector< bool > visited(numNodes, false);
            if(!pathExists(nextNode, currentNode,visited)){
                adj[currentNode].push_back(nextNode);
                currentNode = nextNode;
            } else {
                adj[currentNode].push_back(numNodes);
                adj.push_back({});
                // newToOld[numNodes] = newToOld[nextNode];
                intIdToStringId[numNodes] = intIdToStringId[nextNode];
                sequence.second[i] = numNodes;
                currentNode = nextNode;
                numNodes++;
            }

        }
    }

    // for(auto sequence: intSequences){

    // }

    if(!checkForCycles()){
        std::cout << "Cycles Removed!" << std::endl;
        return;
    } else {
        std::cout << "Error: Cycle removal failed!" << std::endl;
        return;
    }

}

bool PangenomeMAT2::Pangraph::pathExists(size_t nId1, size_t nId2, std::vector< bool >& visited){
    if(nId1 == nId2){
        return true;
    }
    visited[nId1] = true;

    for(auto u: adj[nId1]){
        if(!visited[u] && pathExists(u, nId2, visited)){
            return true;
        }
    }
    return false;
}

bool PangenomeMAT2::Pangraph::checkForCyclesHelper(size_t nodeId, std::vector< int >& color){
    color[nodeId] = 1;
    for(auto u: adj[nodeId]){
        if(color[u] == 1){
            return true;
        }
        if(color[u] == 0){
            if(checkForCyclesHelper(u, color)){
                return true;
            }
        }
    }
    color[nodeId] = 2;
    return false;
}

bool PangenomeMAT2::Pangraph::checkForCycles(){
    std::vector< int > color(numNodes, 0);
    for(size_t i = 0; i < numNodes; i++){
        if(color[i] == 0){
            if(checkForCyclesHelper(i, color)){
                return true;
            }
        }
    }
    return false;
}

std::unordered_map< std::string, std::vector< int > > PangenomeMAT2::Pangraph::getAlignedSequences(const std::vector< size_t >& topoArray){
    std::unordered_map< std::string, std::vector< int > > alignedSequences;
    for(auto p: intSequences){
        size_t p1 = 0, p2 = 0;
        while(p1 < topoArray.size() && p2 < p.second.size()){
            if(topoArray[p1] == p.second[p2]){
                alignedSequences[p.first].push_back(topoArray[p1]);
                p2++;
            } else {
                alignedSequences[p.first].push_back(-1);
            }
            p1++;
        }
        while(alignedSequences[p.first].size() < topoArray.size()){
            alignedSequences[p.first].push_back(-1);
        }
    }
    return alignedSequences;
}

void PangenomeMAT2::Pangraph::topologicalSortHelper(size_t nodeId, std::vector< size_t >& topoArray, std::vector< bool >& visited){
    visited[nodeId] = true;
    for(auto u: adj[nodeId]){
        if(!visited[u]){
            topologicalSortHelper(u, topoArray, visited);
        }
    }
    topoArray.push_back(nodeId);
}

std::vector< size_t > PangenomeMAT2::Pangraph::getTopologicalSort(){
    std::vector< size_t > topoArray;
    std::vector< bool > visited(numNodes, false);

    for(size_t i = 0; i < numNodes; i++){
        if(!visited[i]){
            topologicalSortHelper(i, topoArray, visited);
        }
    }
    std::reverse(topoArray.begin(), topoArray.end());

    return topoArray;
}

PangenomeMAT2::TreeGroup::TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile){
    for(size_t i = 0; i < treeFiles.size(); i++){
        trees.emplace_back(treeFiles[i]);
    }
    // mutation file format: mutation type (H or R), tree_1 index, sequence_1 name, tree_2 index sequence_2 name, start_point_1, end_point_1, start_point_2, end_point_2, tree_3 index (child tree)
    std::string line;
    while(getline(mutationFile, line, '\n')){
        std::vector< std::string > tokens;
        stringSplit(line, ' ', tokens);
        char mutationType = tokens[0][0];
        size_t treeIndex1 = std::stoll(tokens[1]);
        std::string sequenceId1 = tokens[2];
        size_t treeIndex2 = std::stoll(tokens[3]);
        std::string sequenceId2 = tokens[4];
        size_t startPoint1 = std::stoll(tokens[5]);
        size_t endPoint1 = std::stoll(tokens[6]);
        size_t startPoint2 = std::stoll(tokens[7]);
        size_t endPoint2 = std::stoll(tokens[8]);
        size_t treeIndex3 = std::stoll(tokens[9]);

        sequence_t sequence1, sequence2;
        blockExists_t blockExists1, blockExists2;
        trees[treeIndex1].getSequenceFromReference(sequence1, blockExists1, sequenceId1);
        trees[treeIndex2].getSequenceFromReference(sequence2, blockExists2, sequenceId2);

        std::tuple< int,int,int,int > t_start1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(startPoint1, sequence1, blockExists1);
        std::tuple< int,int,int,int > t_end1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(endPoint1, sequence1, blockExists1);
        std::tuple< int,int,int,int > t_start2 = trees[treeIndex1].globalCoordinateToBlockCoordinate(startPoint2, sequence2, blockExists2);
        std::tuple< int,int,int,int > t_end2 = trees[treeIndex1].globalCoordinateToBlockCoordinate(endPoint2, sequence2, blockExists2);

        complexMutations.emplace_back(mutationType, treeIndex1, treeIndex2, treeIndex3, sequenceId1, sequenceId2, t_start1, t_end1, t_start2, t_end2);

    }
}

PangenomeMAT2::TreeGroup::TreeGroup(std::ifstream& fin){
    MATNew::treeGroup TG;
    if(!TG.ParseFromIstream(&fin)){
        throw std::invalid_argument("Could not read tree group from input file.");
    }
    
    for(int i = 0; i < TG.trees_size(); i++){
        trees.emplace_back(TG.trees(i));
    }
}

void PangenomeMAT2::TreeGroup::printFASTA(std::ofstream& fout){
    for(auto& tree: trees){
        tree.printFASTA(fout);
    }
}

void PangenomeMAT2::TreeGroup::writeToFile(std::ofstream& fout){
    MATNew::treeGroup treeGroupToWrite;

    for(auto& tree: trees){
        MATNew::tree treeToWrite;
        Node* node = tree.root;
        tree.getNodesPreorder(node, treeToWrite);
        std::string newick = tree.getNewickString(node);

        treeToWrite.set_newick(newick);
        std::map< std::vector< uint32_t >, std::vector< int64_t > > consensusSeqToBlockIds;
        
        for(auto block: tree.blocks){
            MATNew::block b;
            int64_t blockId;
            if(block.secondaryBlockId != -1){
                blockId = ((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId;
                b.set_blockid(blockId);
                b.set_blockgapexist(true);
            } else {
                blockId = ((int64_t)block.primaryBlockId << 32);
                b.set_blockid(blockId);
                b.set_blockgapexist(false);
            }
            b.set_chromosomename(block.chromosomeName);
            consensusSeqToBlockIds[block.consensusSeq].push_back(blockId);
            // for(auto n: block.consensusSeq){
            //     b.add_consensusseq(n);
            // }
            treeToWrite.add_blocks();
            *treeToWrite.mutable_blocks( treeToWrite.blocks_size() - 1 ) = b;
        }

        for(auto u: consensusSeqToBlockIds){
            MATNew::consensusSeqToBlockIds c;
            for(auto v: u.first){
                c.add_consensusseq(v);
            }
            for(auto v: u.second){
                c.add_blockid(v);
            }
            treeToWrite.add_consensusseqmap();
            *treeToWrite.mutable_consensusseqmap( treeToWrite.consensusseqmap_size() - 1 ) = c;
        }
        for(size_t i = 0; i < tree.gaps.size(); i++){
            MATNew::gapList gl;
            for(size_t j = 0; j < tree.gaps[i].nucPosition.size(); j++){
                gl.add_nucposition(tree.gaps[i].nucPosition[j]);
                gl.add_nucgaplength(tree.gaps[i].nucGapLength[j]);
            }
            if(tree.gaps[i].secondaryBlockId != -1){
                gl.set_blockid(((int64_t)tree.gaps[i].primaryBlockId << 32) + tree.gaps[i].secondaryBlockId);
                gl.set_blockgapexist(true);
            } else {
                gl.set_blockid(((int64_t)tree.gaps[i].primaryBlockId << 32));
                gl.set_blockgapexist(false);
            }
            treeToWrite.add_gaps();
            *treeToWrite.mutable_gaps( treeToWrite.gaps_size() - 1 ) = gl;
        }

        treeGroupToWrite.add_trees();
        *treeGroupToWrite.mutable_trees( treeGroupToWrite.trees_size() - 1 ) = treeToWrite;
    }

    for(auto cm: complexMutations){
        treeGroupToWrite.add_complexmutations();
        *treeGroupToWrite.mutable_complexmutations( treeGroupToWrite.complexmutations_size()-1 ) = cm.toProtobuf();
    }

    if(!treeGroupToWrite.SerializeToOstream(&fout)){
        std::cerr << "Failed to write to output file." << std::endl;
    }
}