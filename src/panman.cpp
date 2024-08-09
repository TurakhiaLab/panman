#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_invoke.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_map.h>
#include <ctime>
#include <iomanip>
#include <mutex>
#include <chrono>
#include <filesystem>
#include "chaining.cpp"
#include "rotation.cpp"
#include "fitchSankoff.cpp"
#include "summary.cpp"
#include "fasta.cpp"
#include "maf.cpp"
#include "vcf.cpp"
#include "subnet.cpp"
#include "gfa.cpp"
#include "annotate.cpp"
#include "reroot.cpp"
#include "aaTrans.cpp"

#include "panmanUtils.hpp"




char panmanUtils::getNucleotideFromCode(int code) {
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

char panmanUtils::getCodeFromNucleotide(char nuc) {
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

std::vector<int> getSankoffVector(char nuc) {
    std::vector<int> sankoffVector(5,SANKOFF_INF);
    switch(nuc) {
    case 'A':
        sankoffVector[1]=0;
    case 'C':
        sankoffVector[2]=0;
    case 'G':
        sankoffVector[3]=0;
    case 'T':
        sankoffVector[4]=0;
    case 'R':
        sankoffVector[1]=0;
        sankoffVector[4]=0;
    case 'Y':
        sankoffVector[2]=0;
        sankoffVector[4]=0;
    case 'S':
        sankoffVector[2]=0;
        sankoffVector[3]=0;
    case 'W':
        sankoffVector[1]=0;
        sankoffVector[4]=0;
    case 'K':
        sankoffVector[3]=0;
        sankoffVector[4]=0;
    case 'M':
        sankoffVector[1]=0;
        sankoffVector[2]=0;
    case 'B':
        sankoffVector[3]=0;
        sankoffVector[4]=0;
    case 'D':
        sankoffVector[1]=0;
        sankoffVector[3]=0;
        sankoffVector[4]=0;
    case 'H':
        sankoffVector[1]=0;
        sankoffVector[2]=0;
        sankoffVector[4]=0;
    case 'V':
        sankoffVector[1]=0;
        sankoffVector[2]=0;
        sankoffVector[3]=0;
    case 'N':
        sankoffVector[1]=0;
        sankoffVector[2]=0;
        sankoffVector[3]=0;
        sankoffVector[4]=0;
    default:
        sankoffVector[0]=0;
    }
    return sankoffVector;
}

// For reverse complement
char panmanUtils::getComplementCharacter(char nuc) {
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

std::string panmanUtils::getDate() {
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);
    std::string date;
    date += std::to_string(now->tm_year + 1900)
            + std::to_string(now->tm_mon + 1)
            +  std::to_string(now->tm_mday);
    return date;
}

panmanUtils::Node::Node(std::string id, float len) {
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

panmanUtils::Node::Node(std::string id, Node* par, float len) {
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

panmanUtils::Block::Block(size_t pBlockId, std::string seq) {
    primaryBlockId = pBlockId;
    secondaryBlockId = -1;
    for(size_t i = 0; i < seq.length(); i+=8) {
        uint32_t currentConsensusSeq = 0;
        for(size_t j = i; j < std::min(i+8, seq.length()); j++) {
            int code = panmanUtils::getCodeFromNucleotide(seq[j]);
            currentConsensusSeq^=(code << (4*(7-(j-i))));
        }
        consensusSeq.push_back(currentConsensusSeq);
    }
}

panmanUtils::Block::Block(int32_t pBlockId, int32_t sBlockId, const std::vector< uint32_t >& seq) {
    primaryBlockId = pBlockId;
    secondaryBlockId = sBlockId;
    consensusSeq = seq;
}

void panmanUtils::stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
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

panmanUtils::Node* panmanUtils::Tree::createTreeFromNewickString(std::string newickString) {
    newickString = panmanUtils::stripString(newickString);

    panmanUtils::Node* newTreeRoot = nullptr;

    std::vector<std::string> leaves;
    std::vector<size_t> numOpen;
    std::vector<size_t> numClose;
    std::vector<std::queue<float>> branchLen (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    stringSplit(newickString, ',', s1);

    for(size_t i = 0; i < s1.size(); i++) {
        if(s1[i].length() == 0) {
            std::cout << s1[i-2] << " " << s1[i-1] << " " << s1[i] << " " << s1[i+1] << " " << s1[i+2] << std::endl;
        }
    }

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

void panmanUtils::Tree::assignMutationsToNodes(Node* root, size_t& currentIndex,
        std::vector< panman::node >& nodes) {
    std::vector< panmanUtils::NucMut > storedNucMutation;

    for(int i = 0; i < nodes[currentIndex].mutations_size(); i++) {
        for(auto nucMut: nodes[currentIndex].mutations(i).nucmutation()) {
            storedNucMutation.push_back( panmanUtils::NucMut(nucMut,
                                         nodes[currentIndex].mutations(i).blockid(),
                                         nodes[currentIndex].mutations(i).blockgapexist()));
        }
    }

    std::vector< panmanUtils::BlockMut > storedBlockMutation;
    for(int i = 0; i < nodes[currentIndex].mutations_size(); i++) {
        panmanUtils::BlockMut tempBlockMut;
        if(nodes[currentIndex].mutations(i).blockmutexist()) {
            tempBlockMut.loadFromProtobuf(nodes[currentIndex].mutations(i));
            storedBlockMutation.push_back(tempBlockMut);
        }
    }

    for(int i = 0; i < nodes[currentIndex].annotations_size(); i++) {
        root->annotations.push_back(nodes[currentIndex].annotations(i));
        annotationsToNodes[nodes[currentIndex].annotations(i)].push_back(root->identifier);
    }

    root->nucMutation = storedNucMutation;
    root->blockMutation = storedBlockMutation;

    for(auto child: root->children) {
        currentIndex++;
        assignMutationsToNodes(child, currentIndex, nodes);
    }

}




bool panmanUtils::Tree::hasPolytomy(Node* node) {
    if(node->children.size() > 2) {
        return true;
    }
    for(auto u: node->children) {
        if(hasPolytomy(u)) {
            return true;
        }
    }
    return false;
}

panmanUtils::Tree::Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype,
                        std::string reference) {
    if(ftype == panmanUtils::FILE_TYPE::GFA) {
        std::map< std::string, std::string > nodes;
        std::map< std::string, std::vector< std::pair<std::string, bool> > > paths;
        std::string line;
        while(getline(fin, line, '\n')) {
            std::vector< std::string > separatedLine;
            stringSplit(line, '\t', separatedLine);
            if(separatedLine[0] == "S") {
                nodes[separatedLine[1]] = separatedLine[2];
            } else if(separatedLine[0] == "P") {
                std::vector< std::string > v;
                std::vector< std::pair< std::string, bool > > currentPath;
                stringSplit(separatedLine[2], ',', v);
                for(size_t i = 0; i < v.size(); i++) {
                    currentPath.push_back(std::make_pair(v[i].substr(0,v[i].size()-1),
                                                         (v[i][v[i].size()-1] == '+')?true:false));
                }
                paths[separatedLine[1]] = currentPath;
            }
        }
        std::vector< std::vector< std::pair< std::string, bool > > > stringSequences;
        std::vector< std::string > sequenceIds;
        for(auto p: paths) {
            sequenceIds.push_back(p.first);
            stringSequences.push_back(p.second);
        }

        GfaGraph g(sequenceIds, stringSequences, nodes);
        std::cout << "Graph without cycles created" << std::endl;

        std::vector< size_t > topoArray = g.getTopologicalSort();

        std::vector< std::vector< int64_t > > alignedSequences = g.getAlignedSequences(topoArray);
        std::vector< std::vector< int > > alignedStrandSequences = g
                .getAlignedStrandSequences(topoArray);


        std::string newickString;
        secondFin >> newickString;
        root = createTreeFromNewickString(newickString);

        std::unordered_map< std::string, std::vector< int64_t > > pathIdToSequence;
        std::unordered_map< std::string, std::vector< int > > pathIdToStrandSequence;
        for(size_t i = 0; i < g.pathIds.size(); i++) {
            pathIdToSequence[g.pathIds[i]] = alignedSequences[i];
            pathIdToStrandSequence[g.pathIds[i]] = alignedStrandSequences[i];
        }

        for(size_t i = 0; i < topoArray.size(); i++) {
            blocks.emplace_back(i, g.intNodeToSequence[topoArray[i]]);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string,
            std::pair< BlockMutationType, bool > > > globalMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;
            for(const auto& u: pathIdToSequence) {
                if(u.second[i] == -1) {
                    states[u.first] = 1;
                } else if(pathIdToStrandSequence[u.first][i]) {
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
    } else if(ftype == panmanUtils::FILE_TYPE::PANGRAPH) {
        std::string newickString;
        secondFin >> newickString;
        Json::Value pangraphData;
        fin >> pangraphData;

        auto start = std::chrono::high_resolution_clock::now();

        panmanUtils::Pangraph pg(pangraphData);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds timing = end -start;
        std::cout << "Psuedo Root Calculated in: " << timing.count() << " nanoseconds \n";

        circularSequences = pg.circularSequences;
        sequenceInverted = pg.sequenceInverted;
        rotationIndexes = pg.rotationIndexes;

        std::vector< size_t > topoArray = pg.getTopologicalSort();
        std::cout << "Length of Pseudo Root: " << topoArray.size() << endl;
        std::unordered_map< std::string, std::vector< int > >
        alignedSequences = pg.getAlignedSequences(topoArray);
        std::unordered_map< std::string, std::vector< int > >
        alignedStrandSequences = pg.getAlignedStrandSequences(topoArray);

        root = createTreeFromNewickString(newickString);

        // Check if tree is a polytomy to check if Sankoff algorithm needs to be applied
        bool polytomy = hasPolytomy(root);

        for(size_t i = 0; i < topoArray.size(); i++) {
            blocks.emplace_back(i, pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]]);
        }

        for(size_t i = 0; i < topoArray.size(); i++) {
            GapList g;
            g.primaryBlockId = i;
            g.secondaryBlockId = -1;
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++) {
                g.nucPosition.push_back(pg
                                        .stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first);
                g.nucGapLength.push_back(pg
                                         .stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second);
            }
            gaps.push_back(g);
        }

        tbb::concurrent_unordered_map< size_t, std::unordered_map< std::string,
            std::pair< BlockMutationType, bool > > > globalBlockMutations;

        std::cout << "Inferring mutations..." << std::endl;
        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            if(!polytomy) {
                // Apply Fitch's algorithm if not a Polytomy
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;

                int defaultState = -1;

                for(const auto& u: alignedSequences) {
                    if(reference.length()) {
                        if(u.first.find(reference) != std::string::npos) {
                            if(u.second[i] == -1) {
                                defaultState = 1;
                            } else if(alignedStrandSequences[u.first][i]) {
                                // forward strand
                                defaultState = 2;
                            } else {
                                // reverse strand
                                defaultState = 4;
                            }
                        }
                    }

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
                if(defaultState != -1) {
                    blockFitchBackwardPassNew(root, states, 1, defaultState);
                } else {
                    blockFitchBackwardPassNew(root, states, 1);
                }
                blockFitchAssignMutationsNew(root, states, mutations, 1);
                globalBlockMutations[i] = mutations;
            } else {
                // Apply Sankoff's algorithm if the tree is a Polytomy

                std::unordered_map< std::string, std::vector < int > > stateSets;
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;

                int defaultState = -1;

                for(const auto& u: alignedSequences) {
                    if(reference.length()) {
                        if(u.first.find(reference) != std::string::npos) {
                            if(u.second[i] == -1) {
                                defaultState = 0;
                            } else if(alignedStrandSequences[u.first][i]) {
                                // forward strand
                                defaultState = 1;
                            } else {
                                // reverse strand
                                defaultState = 2;
                            }
                        }
                    }

                    std::vector< int > currentState(3, SANKOFF_INF);
                    if(u.second[i] == -1) {
                        currentState[0] = 0;
                    } else if(alignedStrandSequences[u.first][i]) {
                        // forward strand
                        currentState[1] = 0;
                    } else {
                        // reverse strand
                        currentState[2] = 0;
                    }
                    stateSets[u.first] = currentState;
                }

                blockSankoffForwardPass(root, stateSets);
                if(defaultState != -1) {
                    blockSankoffBackwardPass(root, stateSets, states, 0, defaultState);
                } else {
                    blockSankoffBackwardPass(root, stateSets, states, 0);
                }
                blockSankoffAssignMutations(root, states, mutations, 0);
                globalBlockMutations[i] = mutations;
            }
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

        tbb::concurrent_unordered_map< std::string,
            std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        tbb::concurrent_unordered_map< std::string,
            std::vector< std::tuple< int,int,int,int,int,int > > > gapMutations;

        tbb::parallel_for((size_t)0, topoArray.size(), [&](size_t i) {
            std::string consensusSeq = pg.stringIdToConsensusSeq[pg.intIdToStringId[topoArray[i]]];
            std::vector< std::pair< char, std::vector< char > > > sequence(consensusSeq.size()+1,
            {'-', {}});
            for(size_t j = 0; j < consensusSeq.length(); j++) {
                sequence[j].first = consensusSeq[j];
            }
            for(size_t j = 0; j < pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]].size(); j++) {
                sequence[pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].first]
                .second.resize(pg.stringIdToGaps[pg.intIdToStringId[topoArray[i]]][j].second,
                               '-');
            }
            tbb::concurrent_unordered_map< std::string,
                std::vector< std::pair< char, std::vector< char > > > > individualSequences;

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
                    if(!polytomy) {
                        // Not a polytomy. Applying Fitch.
                        std::unordered_map< std::string, int > states;
                        std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;

                        int defaultState = -1;
                        for(const auto& u: individualSequences) {
                            if(reference.length()) {
                                if(u.first.find(reference) != std::string::npos) {
                                    if(u.second[j].second[k] != '-') {
                                        defaultState = (1 << getCodeFromNucleotide(u.second[j].second[k]));
                                    } else {
                                        defaultState = 1;
                                    }
                                }
                            }

                            if(u.second[j].second[k] != '-') {
                                states[u.first] = (1 << getCodeFromNucleotide(u.second[j].second[k]));
                            } else {
                                states[u.first] = 1;
                            }
                        }
                        nucFitchForwardPass(root, states);
                        if(defaultState != -1) {
                            nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].second[k])), defaultState);
                        } else {
                            nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                        }
                        nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].second[k])));
                        for(auto mutation: mutations) {
                            nodeMutexes[mutation.first].lock();
                            gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                            nodeMutexes[mutation.first].unlock();
                        }
                    } else {
                        // Since the topology is a polytomy, applying Sankoff
                        std::unordered_map< std::string, std::vector< int > > stateSets;
                        std::unordered_map< std::string, int > states;
                        std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;

                        int defaultState = -1;
                        for(const auto& u: individualSequences) {
                            if(reference.length()) {
                                if(u.first.find(reference) != std::string::npos) {
                                    if(u.second[j].second[k] != '-') {
                                        defaultState = getCodeFromNucleotide(u.second[j].second[k]);
                                    } else {
                                        defaultState = 0;
                                    }
                                }
                            }

                            std::vector< int > currentState(16, SANKOFF_INF);
                            if(u.second[j].second[k] != '-') {
                                currentState[getCodeFromNucleotide(u.second[j].second[k])] = 0;
                            } else {
                                currentState[0] = 0;
                            }
                            stateSets[u.first] = currentState;
                        }
                        nucSankoffForwardPass(root, stateSets);
                        if(defaultState != -1) {
                            nucSankoffBackwardPass(root, stateSets, states, getCodeFromNucleotide(sequence[j].second[k]), defaultState);
                        } else {
                            nucSankoffBackwardPass(root, stateSets, states, getCodeFromNucleotide(sequence[j].second[k]));
                        }
                        nucSankoffAssignMutations(root, states, mutations, getCodeFromNucleotide(sequence[j].second[k]));
                        for(auto mutation: mutations) {
                            nodeMutexes[mutation.first].lock();
                            gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, k, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                            nodeMutexes[mutation.first].unlock();
                        }
                    }
                });

                if(!polytomy) {
                    std::unordered_map< std::string, int > states;
                    std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;
                    int defaultState = -1;

                    for(const auto& u: individualSequences) {
                        if(u.first.find(reference) != std::string::npos) {
                            if(u.second[j].first != '-') {
                                defaultState = (1 << getCodeFromNucleotide(u.second[j].first));
                            } else {
                                defaultState = 1;
                            }
                        }

                        if(u.second[j].first != '-') {
                            states[u.first] = (1 << getCodeFromNucleotide(u.second[j].first));
                        } else {
                            states[u.first] = 1;
                        }
                    }
                    nucFitchForwardPass(root, states);
                    if(defaultState != -1) {
                        nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].first)), defaultState);
                    } else {
                        nucFitchBackwardPass(root, states, (1 << getCodeFromNucleotide(sequence[j].first)));
                    }
                    nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(sequence[j].first)));
                    for(auto mutation: mutations) {
                        nodeMutexes[mutation.first].lock();
                        nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                        nodeMutexes[mutation.first].unlock();
                    }
                } else {
                    // Since the topology is a polytomy, applying Sankoff
                    std::unordered_map< std::string, std::vector< int > > stateSets;
                    std::unordered_map< std::string, int > states;
                    std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;

                    int defaultState = -1;
                    for(const auto& u: individualSequences) {
                        if(reference.length()) {
                            if(u.first.find(reference) != std::string::npos) {
                                if(u.second[j].first != '-') {
                                    defaultState = getCodeFromNucleotide(u.second[j].first);
                                } else {
                                    defaultState = 0;
                                }
                            }
                        }

                        std::vector< int > currentState(16, SANKOFF_INF);
                        if(u.second[j].first != '-') {
                            currentState[getCodeFromNucleotide(u.second[j].first)] = 0;
                        } else {
                            currentState[0] = 0;
                        }
                        stateSets[u.first] = currentState;
                    }

                    nucSankoffForwardPass(root, stateSets);
                    if(defaultState != -1) {
                        nucSankoffBackwardPass(root, stateSets, states, getCodeFromNucleotide(sequence[j].first), defaultState);
                    } else {
                        nucSankoffBackwardPass(root, stateSets, states, getCodeFromNucleotide(sequence[j].first));
                    }
                    nucSankoffAssignMutations(root, states, mutations, getCodeFromNucleotide(sequence[j].first));
                    // for(const auto& u: individualSequences) {
                    //     // if(u.first != "Wuhan-Hu-1,")
                    //     // if((states[u.first] == 0 && u.second[j].first != '-') || (states[u.first] != 0 && getNucleotideFromCode(states[u.first]) != u.second[j].first)) {
                    //     //     std::cout << states[u.first] << " " << getNucleotideFromCode(states[u.first]) << " " << u.second[j].first << std::endl;
                    //     // }
                    //     if(u.first != "Wuhan-Hu-1,")
                    //     if(mutations.find(u.first) == mutations.end()) {
                    //         assert(states[u.first] == 0);
                    //         assert(u.second[j].first == '-');
                    //     } else {
                    //         assert(mutations[u.first].first != NucMutationType::ND);
                    //         assert(getNucleotideFromCode(states[u.first]) == mutations[u.first].second);
                    //         assert(u.second[j].first == mutations[u.first].second);
                    //     }
                    // }
                    // for(const auto& u: individualSequences) {
                    //     if(u.first == "Wuhan-Hu-1,"){
                    //         continue;
                    //     }
                    //     std::vector<Node*> path;
                    //     Node* current = allNodes[u.first];
                    //     while(current != root) {
                    //         path.push_back(current);
                    //         current = current->parent;
                    //     }
                    //     path.push_back(current);
                    //     char nuc = getCodeFromNucleotide(sequence[j].first);
                    //     assert(nuc != '-');
                    //     for(auto itr = path.rbegin(); itr != path.rend(); itr++) {
                    //         if(mutations.find((*itr)->identifier) != mutations.end())
                    //         nuc = mutations[(*itr)->identifier].second;
                    //     }
                    //     assert(nuc == u.second[j].first);
                    // }
                    for(auto mutation: mutations) {
                        nodeMutexes[mutation.first].lock();
                        nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, j, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                        nodeMutexes[mutation.first].unlock();
                    }
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

    } else if(ftype == panmanUtils::FILE_TYPE::MSA) {
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
        root->blockMutation.emplace_back(0, std::make_pair(BlockMutationType::BI, false));

        tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
        std::unordered_map< std::string, std::mutex > nodeMutexes;

        for(auto u: allNodes) {
            nodeMutexes[u.first];
        }

        tbb::parallel_for((size_t)0, consensusSeq.length(), [&](size_t i) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;
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
}

void panmanUtils::Tree::protoMATToTree(const panman::tree& mainTree) {
    // Create tree
    root = createTreeFromNewickString(mainTree.newick());

    std::map< std::pair<int32_t, int32_t>, std::vector< uint32_t > > blockIdToConsensusSeq;

    for(int i = 0; i < mainTree.consensusseqmap_size(); i++) {
        std::vector< uint32_t > seq;
        for(int j = 0; j < mainTree.consensusseqmap(i).consensusseq_size(); j++) {
            seq.push_back(mainTree.consensusseqmap(i).consensusseq(j));
        }
        for(int j = 0; j < mainTree.consensusseqmap(i).blockid_size(); j++) {
            std::pair< int32_t, int32_t > blockId;
            blockId.first = (mainTree.consensusseqmap(i).blockid(j) >> 32);
            if(mainTree.consensusseqmap(i).blockgapexist(j)) {
                blockId.second = (mainTree.consensusseqmap(i).blockid(j) & 0xFFFFFFFF);
            } else {
                blockId.second = -1;
            }
            blockIdToConsensusSeq[blockId] = seq;
        }
    }

    std::vector< panman::node > storedNodes;
    for(int i = 0; i < mainTree.nodes_size(); i++) {
        storedNodes.push_back(mainTree.nodes(i));
    }

    size_t initialIndex = 0;

    assignMutationsToNodes(root, initialIndex, storedNodes);

    // Block sequence
    for(auto u: blockIdToConsensusSeq) {
        blocks.emplace_back(u.first.first, u.first.second, u.second);
    }

    // Gap List
    for(int i = 0; i < mainTree.gaps_size(); i++) {
        panmanUtils::GapList tempGaps;
        tempGaps.primaryBlockId = (mainTree.gaps(i).blockid() >> 32);
        tempGaps.secondaryBlockId = (mainTree.gaps(i).blockgapexist() ? (mainTree.gaps(i).blockid() & 0xFFFF): -1);
        for(int j = 0; j < mainTree.gaps(i).nucposition_size(); j++) {
            tempGaps.nucPosition.push_back(mainTree.gaps(i).nucposition(j));
            tempGaps.nucGapLength.push_back(mainTree.gaps(i).nucgaplength(j));
        }
        gaps.push_back(tempGaps);
    }

    // Circular offsets
    for(int i = 0; i < mainTree.circularsequences_size(); i++) {
        circularSequences[mainTree.circularsequences(i).sequenceid()] = mainTree.circularsequences(i).offset();
    }

    // Rotation Indexes
    for(int i = 0; i < mainTree.rotationindexes_size(); i++) {
        rotationIndexes[mainTree.rotationindexes(i).sequenceid()] = mainTree
                .rotationindexes(i).blockoffset();
    }

    // Sequence inverted
    for(int i = 0; i < mainTree.sequencesinverted_size(); i++) {
        sequenceInverted[mainTree.sequencesinverted(i).sequenceid()] = mainTree
                .sequencesinverted(i).inverted();
    }

    // Block gap list
    for(int i = 0; i < mainTree.blockgaps().blockposition_size(); i++) {
        blockGaps.blockPosition.push_back(mainTree.blockgaps().blockposition(i));
        blockGaps.blockGapLength.push_back(mainTree.blockgaps().blockgaplength(i));
    }

}

panmanUtils::Tree::Tree(const panman::tree& mainTree) {
    protoMATToTree(mainTree);
}

panmanUtils::Tree::Tree(std::istream& fin, FILE_TYPE ftype) {

    if(ftype == panmanUtils::FILE_TYPE::PANMAT) {
        panman::tree mainTree;
        if(!mainTree.ParseFromIstream(&fin)) {
            throw std::invalid_argument("Could not read tree from input file.");
        }

        protoMATToTree(mainTree);
    }
}


void panmanUtils::Tree::printBfs(Node* node) {
    if(node == nullptr) {
        node = root;
    }

    // Traversal test
    std::queue<Node *> bfsQueue;
    size_t prevLev = 0;

    bfsQueue.push(node);

    while(!bfsQueue.empty()) {
        Node* current = bfsQueue.front();
        bfsQueue.pop();

        if(current->level != prevLev) {
            std::cout << '\n';
            prevLev = current->level;
        }
        std::cout << '(' << current->identifier << "," << current->branchLength << ") ";

        for(auto child: current->children) {
            bfsQueue.push(child);
        }
    }
    std::cout << '\n';
}




void panmanUtils::Tree::dfsExpansion(panmanUtils::Node* node,
                                     std::vector< panmanUtils::Node* >& vec) {
    vec.push_back(node);
    for(auto child: node->children) {
        dfsExpansion(child, vec);
    }
}

std::string panmanUtils::Tree::getNewickString(Node* node) {

    std::vector< panmanUtils::Node* > traversal;
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
            newick += branch_length_stack.top();
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    newick += ';';

    return newick;
}

// Merge parent node and child node into parent node
void panmanUtils::Tree::mergeNodes(panmanUtils::Node* par, panmanUtils::Node* chi) {

    par->identifier = chi->identifier;
    par->branchLength += chi->branchLength;
    par->children = chi->children;

    // For block mutations, we cancel out irrelevant mutations
    std::map< std::pair<int, int>, std::pair< panmanUtils::BlockMutationType, bool > > bidMutations;

    for(auto mutation: par->blockMutation) {
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        bool type = (mutation.blockMutInfo);
        bool inversion = (mutation.inversion);

        if(type == panmanUtils::BlockMutationType::BI) {
            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair( panmanUtils::BlockMutationType::BI, inversion );
        } else {
            if(bidMutations.find(std::make_pair(primaryBlockId, secondaryBlockId)) != bidMutations.end()) {
                if(bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].first == panmanUtils::BlockMutationType::BI) {
                    // If it was insertion earlier
                    if(inversion) {
                        // This means that the new mutation is an inversion. So, inverted the strand of the inserted block
                        bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = !bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second;
                    } else {
                        // Actually a deletion. So insertion and deletion cancel out
                        bidMutations.erase(std::make_pair(primaryBlockId, secondaryBlockId));
                    }
                } else {
                    // If previous mutation was an inversion
                    if(!inversion) {
                        // Actually a deletion. Remove inversion mutation and put deletion instead
                        bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = false;
                    }
                    // deletion followed by inversion doesn't make sense
                }
            } else {
                bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair(panmanUtils::BlockMutationType::BD, inversion);
            }
        }
    }

    for(auto mutation: chi->blockMutation) {
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int type = (mutation.blockMutInfo);
        bool inversion = (mutation.inversion);

        if(type == panmanUtils::BlockMutationType::BI) {
            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair(panmanUtils::BlockMutationType::BI, inversion);
        } else {
            if(bidMutations.find(std::make_pair(primaryBlockId, secondaryBlockId)) != bidMutations.end()) {
                if(bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].first == panmanUtils::BlockMutationType::BI) {
                    // If it was insertion earlier
                    if(inversion) {
                        // This means that the new mutation is an inversion. So, inverted the strand of the inserted block
                        bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = !bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second;
                    } else {
                        // Actually a deletion. So insertion and deletion cancel out
                        bidMutations.erase(std::make_pair(primaryBlockId, secondaryBlockId));
                    }
                } else {
                    // If previous mutation was an inversion
                    if(!inversion) {
                        // Actually a deletion. Remove inversion mutation and put deletion instead
                        bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = false;
                    }
                    // deletion followed by inversion doesn't make sense
                }
            } else {
                bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair(panmanUtils::BlockMutationType::BD, inversion);
            }
        }
    }

    std::vector< panmanUtils::BlockMut > newBlockMutation;
    for(auto mutation: bidMutations) {
        if(mutation.second.first == panmanUtils::BlockMutationType::BI) {
            panmanUtils::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = panmanUtils::BlockMutationType::BI;
            tempBlockMut.inversion = mutation.second.second;
            newBlockMutation.push_back( tempBlockMut );
        } else {
            panmanUtils::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = panmanUtils::BlockMutationType::BD;
            tempBlockMut.inversion = mutation.second.second;
            newBlockMutation.push_back( tempBlockMut );
        }
    }

    par->blockMutation = newBlockMutation;

    for(auto mutation: chi->nucMutation) {
        par->nucMutation.push_back(mutation);
    }

    delete chi;
}

// Replace old < type, nuc > pair with new < type, nuc > pair
std::pair< int, int > panmanUtils::Tree::replaceMutation(std::pair<int,int> oldMutation, std::pair<int, int> newMutation) {
    std::pair<int, int> ans = newMutation;
    if(oldMutation.first == newMutation.first) {
        ans = newMutation;
    } else if(oldMutation.first == panmanUtils::NucMutationType::NSNPS) {
        // Insertion after substitution (doesn't make sense but just in case)
        if(newMutation.first == panmanUtils::NucMutationType::NSNPI) {
            ans.first = panmanUtils::NucMutationType::NSNPS;
        } else if(newMutation.first == panmanUtils::NucMutationType::NSNPD) {
            ans = newMutation;
        }
    } else if(oldMutation.first == panmanUtils::NucMutationType::NSNPI) {
        if(newMutation.first == panmanUtils::NucMutationType::NSNPS) {
            ans.first = panmanUtils::NucMutationType::NSNPI;
        } else if(newMutation.first == panmanUtils::NucMutationType::NSNPD) {
            // Cancel out the two mutations if deletion after insertion
            ans = std::make_pair(404, 404);
        }
    } else if(oldMutation.first == panmanUtils::NucMutationType::NSNPD) {
        if(newMutation.first == panmanUtils::NucMutationType::NSNPI) {
            ans.first = panmanUtils::NucMutationType::NSNPS;
        } else if(newMutation.first == panmanUtils::NucMutationType::NSNPS) {
            // Substitution after deletion. Doesn't make sense but still
            ans.first = panmanUtils::NucMutationType::NSNPI;
        }
    }
    return ans;
}

bool panmanUtils::Tree::debugSimilarity(const std::vector< panmanUtils::NucMut > array1, const std::vector< panmanUtils::NucMut > array2) {
    std::map< std::tuple< int, int, int, int >, std::pair< int, int > > mutationRecords1, mutationRecords2;

    for(auto mutation: array1) {
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = ((mutation.mutInfo) >> 4);

        if(type >= 3) {
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type) {
        case panmanUtils::NucMutationType::NS:
            newType = panmanUtils::NucMutationType::NSNPS;
            break;
        case panmanUtils::NucMutationType::ND:
            newType = panmanUtils::NucMutationType::NSNPD;
            break;
        case panmanUtils::NucMutationType::NI:
            newType = panmanUtils::NucMutationType::NSNPI;
            break;
        }

        for(int i = 0; i < len; i++) {
            int newChar = (((mutation.nucs) >> (4*(5 - i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1) {
                if(mutationRecords1.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords1.end()) {
                    mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
                        mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords1.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords1.end()) {
                    mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
                        mutationRecords1[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords1.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    for(auto mutation: array2) {
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = ((mutation.mutInfo) >> 4);

        if(type >= 3) {
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type) {
        case panmanUtils::NucMutationType::NS:
            newType = panmanUtils::NucMutationType::NSNPS;
            break;
        case panmanUtils::NucMutationType::ND:
            newType = panmanUtils::NucMutationType::NSNPD;
            break;
        case panmanUtils::NucMutationType::NI:
            newType = panmanUtils::NucMutationType::NSNPI;
            break;
        }

        for(int i = 0; i < len; i++) {
            int newChar = (((mutation.nucs) >> (4*(5 - i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1) {
                if(mutationRecords2.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords2.end()) {
                    mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
                        mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords2.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords2.end()) {
                    mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
                        mutationRecords2[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                    } else {
                        mutationRecords2.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos ));
                    }
                }
            }
        }
    }

    std::vector< std::tuple< int, int, int, int, int, int > > mutationArray1, mutationArray2;
    for(auto u: mutationRecords1) {
        mutationArray1.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }
    for(auto u: mutationRecords2) {
        mutationArray2.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }

    if(mutationArray1.size() != mutationArray2.size()) {
        std::cout << "sizes don't match " << mutationArray1.size() << " " << mutationArray2.size() << std::endl;
        return false;
    }

    for(size_t i = 0; i < mutationArray1.size(); i++) {
        if(mutationArray1[i] != mutationArray2[i]) {
            std::cout << i << "th index doesn't match" << std::endl;
            return false;
        }
    }

    return true;
}

std::vector< panmanUtils::NucMut > panmanUtils::Tree::consolidateNucMutations(const std::vector< panmanUtils::NucMut >& nucMutation) {
    // primaryBid, secondaryBid, pos, gap_pos -> type, nuc
    std::map< std::tuple< int32_t, int32_t, int32_t, int32_t >, std::pair< int, int > > mutationRecords;
    for(auto mutation: nucMutation) {
        int primaryBlockId = mutation.primaryBlockId;
        int secondaryBlockId = mutation.secondaryBlockId;
        int pos = mutation.nucPosition;
        int gapPos = mutation.nucGapPosition;

        // I'm using int instead of NucMutationType because I want the 404 mutation too.
        int type = ((mutation.mutInfo) & 0x7);
        int len = (((mutation.mutInfo) >> 4));

        if(type >= 3) {
            len = 1;
        }

        // Replace variable length mutations into SNP. They will be combined later
        int newType = type;
        switch(type) {
        case panmanUtils::NucMutationType::NS:
            newType = panmanUtils::NucMutationType::NSNPS;
            break;
        case panmanUtils::NucMutationType::ND:
            newType = panmanUtils::NucMutationType::NSNPD;
            break;
        case panmanUtils::NucMutationType::NI:
            newType = panmanUtils::NucMutationType::NSNPI;
            break;
        }

        for(int i = 0; i < len; i++) {
            int newChar = (((mutation.nucs) >> (4*(5-i))) & 0xF);

            std::pair< int, int > newMutation = std::make_pair( newType, newChar );
            if(gapPos != -1) {
                if(mutationRecords.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )) == mutationRecords.end()) {
                    mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
                        mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i )] = newMutation;
                    } else {
                        mutationRecords.erase(std::make_tuple( primaryBlockId, secondaryBlockId, pos, gapPos + i ));
                    }
                }
            } else {
                if(mutationRecords.find(std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )) == mutationRecords.end()) {
                    mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )] = newMutation;
                } else {
                    std::pair< int, int > oldMutation = mutationRecords[std::make_tuple( primaryBlockId, secondaryBlockId, pos + i, gapPos )];
                    newMutation = replaceMutation(oldMutation, newMutation);
                    if(newMutation.first != 404) {
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
    for(auto u: mutationRecords) {
        mutationArray.push_back( std::make_tuple( std::get<0>(u.first), std::get<1>(u.first), std::get<2>(u.first), std::get<3>(u.first), u.second.first, u.second.second ) );
    }

    // mutation array is already sorted since mutationRecord was sorted
    std::vector< panmanUtils::NucMut > consolidatedMutationArray;

    for(size_t i = 0; i < mutationArray.size(); i++) {
        size_t j = i + 1;
        for(; j < std::min(i + 6, mutationArray.size()); j++) {
            if(std::get<3>(mutationArray[i]) != -1) {
                // gapPos exists
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && std::get<1>(mutationArray[i]) == std::get<1>(mutationArray[j]) && std::get<2>(mutationArray[i]) == std::get<2>(mutationArray[j])
                        && std::get<4>(mutationArray[i]) == std::get<4>(mutationArray[j]) && (size_t)(std::get<3>(mutationArray[j]) - std::get<3>(mutationArray[i])) == j - i)) {
                    break;
                }
            } else {
                if(!(std::get<0>(mutationArray[i]) == std::get<0>(mutationArray[j]) && std::get<1>(mutationArray[i]) == std::get<1>(mutationArray[j]) && (size_t)(std::get<2>(mutationArray[j]) - std::get<2>(mutationArray[i])) == j - i
                        && std::get<4>(mutationArray[i]) == std::get<4>(mutationArray[j]) && std::get<3>(mutationArray[j]) == std::get<3>(mutationArray[i]))) {
                    break;
                }
            }
        }

        if(j - i <= 1) {
            consolidatedMutationArray.push_back(panmanUtils::NucMut(mutationArray[i]));
            continue;
        }
        // combine mutations from i to j
        auto newMutation = panmanUtils::NucMut(mutationArray, i, j);

        consolidatedMutationArray.push_back(newMutation);

        i = j - 1;
    }

    return consolidatedMutationArray;
}


bool panmanUtils::Tree::panMATCoordinateGeq(const std::tuple< int, int, int, int >& coor1,
        const std::tuple< int, int, int, int >& coor2, bool strand) {

    if(coor1 == coor2) {
        return true;
    }

    if(std::get<0>(coor1) > std::get<0>(coor2)) {
        return true ^ (!strand);
    } else if (std::get<0>(coor1) < std::get<0>(coor2)) {
        return false ^ (!strand);
    }

    if(std::get<2>(coor1) > std::get<2>(coor2)) {
        return true ^ (!strand);
    } else if(std::get<2>(coor1) < std::get<2>(coor2)) {
        return false ^ (!strand);
    }

    if(std::get<3>(coor1) == -1) {
        // This means that the gap coordinate of the first coordinate is -1 but that of the second
        // coordinate is not
        return true ^ (!strand);
    } else if(std::get<3>(coor2) == -1) {
        // This means that the gap coordinate of the second coordinate is -1 but that of the first
        // coordinate is not
        return false ^ (!strand);
    }

    return (std::get<3>(coor1) > std::get<3>(coor2)) ^ (!strand);
}

bool panmanUtils::Tree::panMATCoordinateLeq(const std::tuple< int, int, int, int >& coor1,
        const std::tuple< int, int, int, int >& coor2, bool strand) {

    if(coor1 == coor2) {
        return true;
    }

    if(std::get<0>(coor1) < std::get<0>(coor2)) {
        return true ^ (!strand);
    } else if (std::get<0>(coor1) > std::get<0>(coor2)) {
        return false ^ (!strand);
    }

    if(std::get<2>(coor1) < std::get<2>(coor2)) {
        return true ^ (!strand);
    } else if(std::get<2>(coor1) > std::get<2>(coor2)) {
        return false ^ (!strand);
    }

    if(std::get<3>(coor1) == -1) {
        // This means that the gap coordinate of the first coordinate is -1 but that of the second
        // coordinate is not
        return false ^ (!strand);
    } else if(std::get<3>(coor2) == -1) {
        // This means that the gap coordinate of the second coordinate is -1 but that of the first
        // coordinate is not
        return true ^ (!strand);
    }

    return (std::get<3>(coor1) < std::get<3>(coor2)) ^ (!strand);
}

panmanUtils::Node* panmanUtils::Tree::extractPanMATSegmentHelper(panmanUtils::Node* node,
        const std::tuple< int, int, int, int >& start, const std::tuple< int, int, int, int >& end,
        const blockStrand_t& rootBlockStrand) {

    int startBlockID = std::get<0>(start);
    int endBlockID = std::get<0>(end);

    panmanUtils::Node* newNode = new panmanUtils::Node(node->identifier, node->branchLength);

    // Push all nucleotide mutations in range to the new node
    for(auto mutation: node->nucMutation) {
        std::tuple< int, int, int, int > mutationCoordinate = std::make_tuple(
                    mutation.primaryBlockId, -1, mutation.nucPosition, mutation.nucGapPosition);
        bool strand = rootBlockStrand[mutation.primaryBlockId].first;

        if(panMATCoordinateGeq(mutationCoordinate, start, strand)
                && panMATCoordinateLeq(mutationCoordinate, end, strand)) {
            if(mutation.primaryBlockId == startBlockID) {
                if(strand) {
                    mutation.nucPosition -= std::get<2>(start);
                }
            } else if(mutation.primaryBlockId == endBlockID) {
                if(!strand) {
                    mutation.nucPosition -= std::get<2>(end);
                }
            }

            mutation.primaryBlockId -= startBlockID;
            newNode->nucMutation.push_back(mutation);
        } else {
            int type = ((mutation.mutInfo) & 0x7);
            if(type < 3) {
                int len = ((mutation.mutInfo) >> 4);
                std::tuple< int, int, int, int > maxMutCoordinate = mutationCoordinate;
                if(mutation.nucGapPosition == -1) {
                    std::get<2>(maxMutCoordinate) += len-1;
                    if(panMATCoordinateGeq(maxMutCoordinate, start, strand)
                            && panMATCoordinateLeq(maxMutCoordinate, end, strand)) {

                        if(newNode->identifier == "node_1" || newNode->identifier=="node_7" || newNode->identifier =="node_758" || newNode->identifier=="node_760" || newNode->identifier=="node_761" || newNode->identifier=="node_810" || newNode->identifier=="node_811" || newNode->identifier=="USA/MI-CDC-ASC210606604/2021|OL969243.1|2021-12-03") {
                            std::cout << mutation.primaryBlockId << " " << mutation.nucPosition << " " << len << " " << mutation.nucGapPosition << std::endl;
                        }

                        if(mutation.primaryBlockId == startBlockID) {
                            if(strand) {
                                mutation.nucPosition -= std::get<2>(start);
                            }
                        } else if(mutation.primaryBlockId == endBlockID) {
                            if(!strand) {
                                mutation.nucPosition -= std::get<2>(end);
                            }
                        }

                        mutation.primaryBlockId -= startBlockID;

                        int diff = std::get<2>(start) - std::get<2>(mutationCoordinate);
                        mutation.nucPosition += diff;
                        len -= diff;
                        mutation.nucs <<= diff*4;
                        mutation.nucs &= 0xFFFFFF;
                        mutation.mutInfo = (len << 4);
                        mutation.mutInfo |= type;
                        std::cout << "NP " << diff << " " << len << " " << mutation.nucPosition << " " <<  (int)mutation.mutInfo << " " << mutation.nucs << std::endl;

                        newNode->nucMutation.push_back(mutation);
                    }
                } else {
                    std::get<3>(maxMutCoordinate) += len-1;
                    if(panMATCoordinateGeq(maxMutCoordinate, start, strand)
                            && panMATCoordinateLeq(maxMutCoordinate, end, strand)) {

                        if(newNode->identifier == "node_1" || newNode->identifier=="node_7" || newNode->identifier =="node_758" || newNode->identifier=="node_760" || newNode->identifier=="node_761" || newNode->identifier=="node_810" || newNode->identifier=="node_811" || newNode->identifier=="USA/MI-CDC-ASC210606604/2021|0L969243-1|2021-12-03") {
                            std::cout << mutation.primaryBlockId << " " << mutation.nucPosition << " " << len << " " << mutation.nucGapPosition << std::endl;
                        }


                        if(mutation.primaryBlockId == startBlockID) {
                            if(strand) {
                                mutation.nucPosition -= std::get<2>(start);
                            }
                        } else if(mutation.primaryBlockId == endBlockID) {
                            if(!strand) {
                                mutation.nucPosition -= std::get<2>(end);
                            }
                        }

                        mutation.primaryBlockId -= startBlockID;

                        int diff = std::get<3>(maxMutCoordinate) - std::get<3>(start);
                        mutation.nucGapPosition += diff;
                        len -= diff;
                        mutation.nucs <<= 24-4*(len);
                        mutation.nucs &= 0xFFFFFF;
                        mutation.mutInfo = (len << 4);
                        mutation.mutInfo &= type;

                        newNode->nucMutation.push_back(mutation);
                    }
                }
            }
        }
    }

    // Push all block mutations in range to the new node
    for(auto mutation: node->blockMutation) {
        if(mutation.primaryBlockId >= std::get<0>(start)
                && mutation.primaryBlockId <= std::get<0>(end)) {
            mutation.primaryBlockId -= startBlockID;
            newNode->blockMutation.push_back(mutation);
        }
    }

    newNode->children.resize(node->children.size(), nullptr);

    tbb::parallel_for((size_t)0, node->children.size(), [&](size_t i) {
        panmanUtils::Node* child = extractPanMATSegmentHelper(node->children[i], start, end,
                                   rootBlockStrand);
        newNode->children[i] = child;
        child->parent = newNode;
    });

    return newNode;

}

void panmanUtils::Tree::extractPanMATSegment(std::ostream& fout, int64_t start, int64_t end) {
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;

    // Extract Root Sequence
    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // Get PanMAT coordinates from global coordinates
    std::tuple< int, int, int, int > panMATStart = globalCoordinateToBlockCoordinate(start,
            rootSequence, rootBlockExists, rootBlockStrand);
    std::tuple< int, int, int, int > panMATEnd = globalCoordinateToBlockCoordinate(end,
            rootSequence, rootBlockExists, rootBlockStrand);

    std::cout << std::get<0>(panMATStart) << " " << std::get<2>(panMATStart) << " " << std::get<3>(panMATStart) << std::endl;

    panmanUtils::Node* newRoot = extractPanMATSegmentHelper(root, panMATStart, panMATEnd,
                                 rootBlockStrand);
    adjustLevels(newRoot);

    // New Block List
    std::vector< panmanUtils::Block > newBlocks;

    // New Gap List
    std::vector< panmanUtils::GapList > newGaps;

    // First block
    int firstBlockID = std::get<0>(panMATStart);
    int firstNucPosition = std::get<2>(panMATStart);
    int firstNucGapPosition = std::get<3>(panMATStart);

    // Last block
    int lastBlockID = std::get<0>(panMATEnd);
    int lastNucPosition = std::get<2>(panMATEnd);
    int lastNucGapPosition = std::get<3>(panMATEnd);

    // Create consensus sequence of first block
    std::string consensusSequence;
    for(size_t i = 0; i < blocks[firstBlockID].consensusSeq.size(); i++) {
        bool endFlag = false;
        for(size_t j = 0; j < 8; j++) {
            const int nucCode = (((blocks[firstBlockID].consensusSeq[i]) >> (4*(7 - j))) & 15);

            if(nucCode == 0) {
                endFlag = true;
                break;
            }
            const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

            consensusSequence += nucleotide;
        }

        if(endFlag) {
            break;
        }
    }

    if(rootBlockStrand[firstBlockID].first) {
        consensusSequence = consensusSequence.substr(firstNucPosition);
    } else {
        consensusSequence = consensusSequence.substr(0, firstNucPosition + 1);
    }

    newBlocks.emplace_back(0, consensusSequence);

    size_t currentBlockID = 1;
    for(size_t i = firstBlockID + 1; i < lastBlockID; i++) {
        newBlocks.push_back(blocks[i]);
        newBlocks[newBlocks.size() - 1].primaryBlockId = currentBlockID;
        currentBlockID++;
    }

    // Create consensus sequence of last block
    consensusSequence = "";
    for(size_t i = 0; i < blocks[lastBlockID].consensusSeq.size(); i++) {
        bool endFlag = false;
        for(size_t j = 0; j < 8; j++) {
            const int nucCode = (((blocks[lastBlockID].consensusSeq[i]) >> (4*(7 - j))) & 15);

            if(nucCode == 0) {
                endFlag = true;
                break;
            }
            const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

            consensusSequence += nucleotide;
        }

        if(endFlag) {
            break;
        }
    }

    if(rootBlockStrand[lastBlockID].first) {
        consensusSequence = consensusSequence.substr(0, lastNucPosition + 1);
    } else {
        consensusSequence = consensusSequence.substr(lastNucPosition);
    }
    newBlocks.emplace_back(currentBlockID, consensusSequence);

    for(size_t i = 0; i < gaps.size(); i++) {
        if(gaps[i].primaryBlockId > firstBlockID && gaps[i].primaryBlockId < lastBlockID) {
            GapList gl;
            gl.primaryBlockId = gaps[i].primaryBlockId - firstBlockID;
            gl.secondaryBlockId = -1;
            gl.nucPosition = gaps[i].nucPosition;
            gl.nucGapLength = gaps[i].nucGapLength;
            newGaps.push_back(gl);
        } else if(gaps[i].primaryBlockId == firstBlockID) {
            GapList gl;
            gl.primaryBlockId = gaps[i].primaryBlockId - firstBlockID;
            gl.secondaryBlockId = -1;
            for(int j = 0; j < gaps[i].nucPosition.size(); j++) {
                if(rootBlockStrand[firstBlockID].first) {
                    if(gaps[i].nucPosition[j] >= firstNucPosition) {
                        gl.nucPosition.push_back(gaps[i].nucPosition[j] - firstNucPosition);
                        gl.nucGapLength.push_back(gaps[i].nucGapLength[j]);
                    }
                } else {
                    if(gaps[i].nucPosition[j] <= firstNucPosition) {
                        gl.nucPosition.push_back(gaps[i].nucPosition[j]);
                        gl.nucGapLength.push_back(gaps[i].nucGapLength[j]);
                    }
                }
            }
            newGaps.push_back(gl);
        } else if(gaps[i].primaryBlockId == lastBlockID) {
            GapList gl;
            gl.primaryBlockId = gaps[i].primaryBlockId - firstBlockID;
            gl.secondaryBlockId = -1;
            for(int j = 0; j < gaps[i].nucPosition.size(); j++) {
                if(rootBlockStrand[lastBlockID].first) {
                    if(gaps[i].nucPosition[j] <= lastNucPosition) {
                        gl.nucPosition.push_back(gaps[i].nucPosition[j]);
                        gl.nucGapLength.push_back(gaps[i].nucGapLength[j]);
                    }
                } else {
                    if(gaps[i].nucPosition[j] >= lastNucPosition) {
                        gl.nucPosition.push_back(gaps[i].nucPosition[j] - lastNucPosition);
                        gl.nucGapLength.push_back(gaps[i].nucGapLength[j]);
                    }
                }
            }
            newGaps.push_back(gl);
        }
    }

    panman::tree treeToWrite;
    getNodesPreorder(newRoot, treeToWrite);

    std::string newick = getNewickString(newRoot);
    std::string newick2 = getNewickString(root);

    treeToWrite.set_newick(newick);

    std::map< std::vector< uint32_t >, std::vector< std::pair< int64_t, bool > > >
    consensusSeqToBlockIds;

    for(auto block: newBlocks) {
        int64_t blockId;
        bool blockGapExists = false;
        if(block.secondaryBlockId != -1) {
            blockId = ((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId;
            blockGapExists = true;
        } else {
            blockId = ((int64_t)block.primaryBlockId << 32);
        }
        consensusSeqToBlockIds[block.consensusSeq].push_back(
            std::make_pair(blockId, blockGapExists));
    }

    for(auto u: consensusSeqToBlockIds) {
        panman::consensusSeqToBlockIds c;
        for(auto v: u.first) {
            c.add_consensusseq(v);
        }
        for(auto v: u.second) {
            c.add_blockid(v.first);
            c.add_blockgapexist(v.second);
        }
        treeToWrite.add_consensusseqmap();
        *treeToWrite.mutable_consensusseqmap( treeToWrite.consensusseqmap_size() - 1 ) = c;
    }

    for(size_t i = 0; i < newGaps.size(); i++) {
        panman::gapList gl;
        for(size_t j = 0; j < newGaps[i].nucPosition.size(); j++) {
            gl.add_nucposition(newGaps[i].nucPosition[j]);
            gl.add_nucgaplength(newGaps[i].nucGapLength[j]);
        }
        gl.set_blockid(((int64_t)newGaps[i].primaryBlockId << 32));
        gl.set_blockgapexist(false);
        treeToWrite.add_gaps();
        *treeToWrite.mutable_gaps( treeToWrite.gaps_size() - 1 ) = gl;
    }

    if (!treeToWrite.SerializeToOstream(&fout)) {
        std::cerr << "Failed to write to output file." << std::endl;
    }


}

void panmanUtils::Tree::getNodesPreorder(panmanUtils::Node* root, panman::tree& treeToWrite) {

    panman::node n;
    std::map< std::pair< int32_t, int32_t >, std::pair< std::vector< panman::nucMut >, int > > blockToMutations;
    std::map< std::pair< int32_t, int32_t >, bool > blockToInversion;

    for(size_t i = 0; i < root->nucMutation.size(); i++) {
        const panmanUtils::NucMut& mutation = root->nucMutation[i];

        panman::nucMut nm;
        nm.set_nucposition(mutation.nucPosition);
        if(mutation.nucGapPosition != -1) {
            nm.set_nucgapposition(mutation.nucGapPosition);
            nm.set_nucgapexist(true);
        } else {
            nm.set_nucgapexist(false);
        }

        nm.set_mutinfo((((mutation.nucs) >> (24 - (mutation.mutInfo >> 4)*4)) << 8) + mutation.mutInfo);
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].first.push_back(nm);
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].second = 2;
    }

    for(size_t i = 0; i < root->blockMutation.size(); i++) {
        const panmanUtils::BlockMut& mutation = root->blockMutation[i];
        blockToMutations[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)].second = mutation.blockMutInfo;
        blockToInversion[std::make_pair(mutation.primaryBlockId, mutation.secondaryBlockId)] = mutation.inversion;
    }

    for(auto u: blockToMutations) {
        panman::mutation mutation;
        mutation.set_blockmutexist((u.second.second != 2));
        mutation.set_blockmutinfo(u.second.second);
        if(u.second.second != 2) {
            // block mutation exists
            mutation.set_blockinversion(blockToInversion[u.first]);
        } else {
            mutation.set_blockinversion(true);
        }

        int32_t primaryBlockId = u.first.first;
        int32_t secondaryBlockId = u.first.second;
        if(secondaryBlockId != -1) {
            mutation.set_blockid(((int64_t)primaryBlockId << 32) + secondaryBlockId);
            mutation.set_blockgapexist(true);
        } else {
            mutation.set_blockid(((int64_t)primaryBlockId << 32));
            mutation.set_blockgapexist(false);
        }
        for(auto v: u.second.first) {
            mutation.add_nucmutation();
            *mutation.mutable_nucmutation(mutation.nucmutation_size() - 1) = v;
        }
        n.add_mutations();
        *n.mutable_mutations(n.mutations_size() - 1) = mutation;
    }

    for(size_t i = 0; i < root->annotations.size(); i++) {
        n.add_annotations(root->annotations[i]);
    }

    treeToWrite.add_nodes();
    *treeToWrite.mutable_nodes( treeToWrite.nodes_size() - 1 ) = n;

    for(auto child: root->children) {
        getNodesPreorder(child, treeToWrite);
    }
}

void getNodesRootedAt(std::set< std::string >& nodeIds, panmanUtils::Node* node) {
    if(node == nullptr) {
        return;
    }

    nodeIds.insert(node->identifier);
    for(auto child: node->children) {
        getNodesRootedAt(nodeIds, child);
    }
}

// Write PanMAT to file
void panmanUtils::Tree::writeToFile(std::ostream& fout, panmanUtils::Node* node) {
    if(node == nullptr) {
        node = root;
    }
    // Get the list of nodes in subtree rooted at `node`, so only their information is recorded in
    // the written PanMAT - useful in the case of subtree extract
    std::set< std::string > nodeIds;
    getNodesRootedAt(nodeIds, node);

    panman::tree treeToWrite;
    getNodesPreorder(node, treeToWrite);

    std::string newick = getNewickString(node);

    treeToWrite.set_newick(newick);

    std::map< std::vector< uint32_t >, std::vector< std::pair< int64_t, bool > > >
    consensusSeqToBlockIds;

    for(auto block: blocks) {
        int64_t blockId;
        bool blockGapExists = false;
        if(block.secondaryBlockId != -1) {
            blockId = ((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId;
            blockGapExists = true;
        } else {
            blockId = ((int64_t)block.primaryBlockId << 32);
        }
        consensusSeqToBlockIds[block.consensusSeq].push_back(
            std::make_pair(blockId, blockGapExists));
    }

    for(auto u: consensusSeqToBlockIds) {
        panman::consensusSeqToBlockIds c;
        for(auto v: u.first) {
            c.add_consensusseq(v);
        }
        for(auto v: u.second) {
            c.add_blockid(v.first);
            c.add_blockgapexist(v.second);
        }
        treeToWrite.add_consensusseqmap();
        *treeToWrite.mutable_consensusseqmap( treeToWrite.consensusseqmap_size() - 1 ) = c;
    }

    for(size_t i = 0; i < gaps.size(); i++) {
        panman::gapList gl;
        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            gl.add_nucposition(gaps[i].nucPosition[j]);
            gl.add_nucgaplength(gaps[i].nucGapLength[j]);
        }
        if(gaps[i].secondaryBlockId != -1) {
            gl.set_blockid(((int64_t)gaps[i].primaryBlockId << 32) + gaps[i].secondaryBlockId);
            gl.set_blockgapexist(true);
        } else {
            gl.set_blockid(((int64_t)gaps[i].primaryBlockId << 32));
            gl.set_blockgapexist(false);
        }
        treeToWrite.add_gaps();
        *treeToWrite.mutable_gaps( treeToWrite.gaps_size() - 1 ) = gl;
    }

    for(auto u: circularSequences) {
        // Check if sequence is a part of the subtree being written
        if(nodeIds.find(u.first) == nodeIds.end()) {
            continue;
        }

        panman::circularOffset co;
        co.set_sequenceid(u.first);
        co.set_offset(u.second);
        treeToWrite.add_circularsequences();
        *treeToWrite.mutable_circularsequences(treeToWrite.circularsequences_size()-1) = co;
    }

    for(auto u: rotationIndexes) {
        // Check if sequence is a part of the subtree being written
        if(nodeIds.find(u.first) == nodeIds.end()) {
            continue;
        }

        panman::rotationIndex ri;
        ri.set_sequenceid(u.first);
        ri.set_blockoffset(u.second);
        treeToWrite.add_rotationindexes();
        *treeToWrite.mutable_rotationindexes(treeToWrite.rotationindexes_size()-1) = ri;
    }

    for(auto u: sequenceInverted) {
        // Check if sequence is a part of the subtree being written
        if(nodeIds.find(u.first) == nodeIds.end()) {
            continue;
        }

        panman::sequenceInverted si;
        si.set_sequenceid(u.first);
        si.set_inverted(u.second);
        treeToWrite.add_sequencesinverted();
        *treeToWrite.mutable_sequencesinverted(treeToWrite.sequencesinverted_size()-1) = si;
    }

    if (!treeToWrite.SerializeToOstream(&fout)) {
        std::cerr << "Failed to write to output file." << std::endl;
    }
}

void panmanUtils::Tree::getBlockSequenceFromReference(block_t& sequence, bool& blockExists, bool& blockStrand, std::string reference, int64_t primaryBlockId, int64_t secondaryBlockId) {
    Node* referenceNode = nullptr;

    for(auto u: allNodes) {
        // printf("%s\n", u.first);
        // std::cerr << u.first << std::endl;
        if(u.second->children.size() == 0 && u.first == reference) {
            referenceNode = u.second;
            break;
        }
    }

    if(referenceNode == nullptr) {
        std::cerr << "Error: Reference sequence with matching name not found!" << std::endl;
        return;
    }

    std::vector< panmanUtils::Node* > path;
    Node* it = referenceNode;

    while(it != root) {
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(root);

    // Get all blocks on the path
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(auto mutation: (*node)->blockMutation) {
            int pBlockId = mutation.primaryBlockId;
            int sBlockId = mutation.secondaryBlockId;

            if(pBlockId != primaryBlockId || sBlockId != secondaryBlockId) {
                continue;
            }

            int type = (mutation.blockMutInfo);
            bool inversion = mutation.inversion;

            if(type == panmanUtils::BlockMutationType::BI) {
                blockExists = true;

                // if insertion of inverted block takes place, the strand is backwards
                blockStrand = !inversion;
            } else {
                if(inversion) {
                    // This is not actually a deletion but an inversion
                    blockStrand = !blockStrand;
                } else {
                    blockExists = false;
                    blockStrand = true;
                }
            }
        }
    }

    if(!blockExists) {
        return;
    }

    for(size_t i = 0; i < blocks.size(); i++) {

        int32_t pBlockId = ((int32_t)blocks[i].primaryBlockId);
        int32_t sBlockId = ((int32_t)blocks[i].secondaryBlockId);

        if(pBlockId != primaryBlockId || sBlockId != secondaryBlockId) {
            continue;
        }

        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

                if(secondaryBlockId != -1) {
                    sequence.push_back({nucleotide, {}});
                } else {
                    sequence.push_back({nucleotide, {}});
                }
            }
            if(endFlag) {
                break;
            }
        }
        sequence.push_back({'x', {}});
    }

    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++) {
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);

        if(primaryBId != primaryBlockId || secondaryBId != secondaryBlockId) {
            continue;
        }

        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];

            sequence[pos].second.resize(len, '-');
        }
    }

    // Apply nucleotide mutations
    for(auto node = path.rbegin(); node != path.rend(); node++) {

        for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {

            int32_t pBlockId = (*node)->nucMutation[i].primaryBlockId;
            int32_t sBlockId = (*node)->nucMutation[i].secondaryBlockId;

            if(pBlockId != primaryBlockId || sBlockId != secondaryBlockId) {
                continue;
            }

            int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
            int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
            uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3) {

                int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                if(type == panmanUtils::NucMutationType::NS) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[nucPosition].second[nucGapPosition+j] = newVal;
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[nucPosition+j].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NI) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[nucPosition].second[nucGapPosition+j] = newVal;
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[nucPosition+j].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::ND) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            sequence[nucPosition].second[nucGapPosition+j] = '-';
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            sequence[nucPosition+j].first = '-';
                        }
                    }
                }
            } else {
                if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(nucGapPosition != -1) {
                        sequence[nucPosition].second[nucGapPosition] = newVal;
                    } else {
                        sequence[nucPosition].first = newVal;
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPI) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(nucGapPosition != -1) {
                        sequence[nucPosition].second[nucGapPosition] = newVal;
                    } else {
                        sequence[nucPosition].first = newVal;
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPD) {
                    if(nucGapPosition != -1) {
                        sequence[nucPosition].second[nucGapPosition] = '-';
                    } else {
                        sequence[nucPosition].first = '-';
                    }
                }
            }
        }
    }
}

void panmanUtils::Tree::printMutations(std::ostream& fout) {

    // Get reference sequence
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;

    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // sequence_t st;
    // blockExists_t bt;
    // blockStrand_t bst;
    // getSequenceFromReference(st, bt, bst, "USA/CA-CDC-QDX21497008/2021|MW666944.1|2021-01-27");

    // std::cout << st[55].first[132].second[40] << std::endl;


    tbb::concurrent_map< std::tuple< int, int, int >, size_t > panMATCoordinateToGlobal;
    tbb::concurrent_map< std::tuple< int, int, int >, char > rootCurrentCharacter;

    tbb::concurrent_unordered_map< std::string,
        std::vector< std::tuple< char, size_t, char, char, bool > > > nodeMutations;

    tbb::concurrent_unordered_map< size_t, bool > rootPresentBlocks;

    tbb::concurrent_map< std::tuple< int, int, int >, bool > isGapCoordinate;

    // convert PanMAT coordinate to global reference coordinate
    size_t rootCtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            rootPresentBlocks[i] = true;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            // if(rootCtr == 240) {
                            //     std::cout << rootSequence[i].first[j].first << " " << i << " " << j << " " << k << std::endl;
                            // }
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                }
            }
        } else {
            rootPresentBlocks[i] = false;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                }
            }
        }
    }

    // tbb::concurrent_map< std::tuple< std::string, int, int, int >, char > seqChar;

    nodeMutations[root->identifier];
    int countNodeID = 0;
    tbb::parallel_for_each(allNodes, [&](auto u) {
        // for (auto &u: allNodes) {
        std::map< std::tuple< std::string, int, int, int >, char > seqChar;
        std::cout << u.first << "\t" << "\t" << countNodeID++;
        sequence_t st;
        blockExists_t bt;
        blockStrand_t bst;
        getSequenceFromReference(st, bt, bst, u.first);

        printf("\t%s", "seqChar");

        for(size_t i = 0; i < st.size(); i++) {
            if(bst[i].first) {
                for(size_t j = 0; j < st[i].first.size(); j++) {
                    for(size_t k = 0; k < st[i].first[j].second.size(); k++) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x' && bt[i].first) {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = st[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                }
            } else {
                for(size_t j = st[i].first.size() - 1; j + 1 > 0; j--) {
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x') {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = rootSequence[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                    for(size_t k = st[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                }
            }
        }
        // });
        // }


        // Compute mutations for each of the other sequences
        // tbb::parallel_for_each(allNodes, [&](auto u) {
        // countNodeID = 0;
        // for (auto &u: allNodes) {
        // std::cout << "Iteration 1: " << countNodeID++ << std::endl;

        if(u.first == root->identifier) {
            return;
        }

        tbb::concurrent_map< std::tuple< int, int, int >, char > currentCharacter = rootCurrentCharacter;
        tbb::concurrent_unordered_map< size_t, bool > presentBlocks = rootPresentBlocks;

        nodeMutations[u.first];

        Node* it = u.second;
        std::vector< panmanUtils::Node* > pathBlock;
        std::vector< panmanUtils::Node* > path;

        // while(it != root) {
        path.push_back(it);
        // it = it->parent;
        // }

        while(it != root) {
            pathBlock.push_back(it);
            it = it->parent;
        }

        std::vector< std::pair< size_t, std::tuple< char, size_t, char, char, bool > > > currentNodeMutations;

        printf("\t%s", "blockMut");
        for(auto node = pathBlock.rbegin(); node != pathBlock.rend(); node++) {
            for(auto mutation: (*node)->blockMutation) {
                int32_t primaryBlockId = mutation.primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                bool type = mutation.blockMutInfo;
                bool inversion = mutation.inversion;
                if(type == 1) {
                    // insertion
                    presentBlocks[primaryBlockId] = true;
                    if(inversion) {
                        std::cout << "INVERTED BLOCK FOUND" << std::endl;
                    }
                } else {
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        std::cout << "INVERSION FOUND" << std::endl;
                    } else {
                        // Actually a deletion
                        presentBlocks[primaryBlockId] = false;
                    }
                }
            }
        }

        printf("\t%s", "nucMut");
        for(auto node = path.rend()-1; node != path.rend(); node++) {
            for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {
                int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
                int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
                uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';

                if(type < 3) {
                    int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                    if(type == panmanUtils::NucMutationType::NS) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)];
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition+j)];
                                    if(oldVal == '-' || oldVal == 'x') {
                                        continue;
                                    }
                                    if(node == path.rend()-1)
                                        currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition+j, -1)];
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)];
                                    // if(u.first == "Denmark/DCGC-504971/2022|OX187739.1|2022-04-28" && panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)] == 185) {
                                    //     std::cout << oldVal << " " << newVal << " " << primaryBlockId << " " << nucPosition+j << " " << -1 << std::endl;
                                    // }
                                    if(oldVal == '-' || oldVal == 'x') {
                                        // std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                                        continue;
                                    }
                                    // if(oldVal == newVal) {
                                    //     std::cout << primaryBlockId << " " << nucPosition << " " << nucGapPosition+j << " " << oldVal << " " << newVal << std::endl;
                                    // }
                                    if(node == path.rend()-1)
                                        currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, -1)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)] = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));

                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition + j)];
                                if(node == path.rend()-1) {
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                                }
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition + j, nucGapPosition)];
                                if(node == path.rend()-1) {
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));
                                }
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(nucGapPosition != -1) {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)] = newVal;
                    } else {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, -1)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, -1)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)] = newVal;
                    }
                }
            }
        }

        printf("\t%s\n", "Adding");
        for(auto mut: currentNodeMutations) {
            if(presentBlocks.find(mut.first) != presentBlocks.end()) {
                nodeMutations[u.first].push_back(mut.second);
            }
        }
    });
    // }

    printf("%s\n", "Printting mutations");

    for(auto& u: nodeMutations) {
        // print all substitutions first
        fout << "Substitutions:\t";
        // fout << u.first << '\t';
        for(auto v: u.second) {
            if(std::get<0>(v) == 'S') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<2>(v) << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Insertions:\t";
        // fout << u.first << '\t';
        // print insertions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'I') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Deletions:\t";
        // fout << u.first << '\t';
        // print deletions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'D') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<2>(v);
            }
        }
        fout << '\n';
    }

}



void panmanUtils::Tree::printNodePaths(std::ostream& fout) {

    for (auto &u: allNodes) {
        Node* it = u.second;

        while(it != root) {
            std::cout << it->identifier << "\t";
            it = it->parent;
            if (it != root) std::cout << "<\t";
            else std::cout << "\n";
        }

    }

}


void panmanUtils::Tree::printMutationsNew(std::ostream& fout) {

    // Get reference sequence
    sequence_t rootSequence;
    blockExists_t rootBlockExists;
    blockStrand_t rootBlockStrand;
    getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand, root->identifier);

    // sequence_t st;
    // blockExists_t bt;
    // blockStrand_t bst;
    // getSequenceFromReference(st, bt, bst, "USA/CA-CDC-QDX21497008/2021|MW666944.1|2021-01-27");

    // std::cout << st[55].first[132].second[40] << std::endl;


    tbb::concurrent_map< std::tuple< int, int, int >, size_t > panMATCoordinateToGlobal;
    tbb::concurrent_map< std::tuple< int, int, int >, char > rootCurrentCharacter;

    tbb::concurrent_unordered_map< std::string,
        std::vector< std::tuple< char, size_t, char, char, bool > > > nodeMutations;

    tbb::concurrent_unordered_map< size_t, bool > rootPresentBlocks;

    tbb::concurrent_map< std::tuple< int, int, int >, bool > isGapCoordinate;

    // convert PanMAT coordinate to global reference coordinate
    size_t rootCtr = 0;
    for(size_t i = 0; i < rootSequence.size(); i++) {
        if(rootBlockExists[i].first) {
            rootPresentBlocks[i] = true;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    if(rootSequence[i].first[j].first != '-' && rootSequence[i].first[j].first != 'x') {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = rootSequence[i].first[j].first;
                        isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                        rootCtr++;
                    } else {
                        rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,-1)] = true;
                    }
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        if(rootSequence[i].first[j].second[k] != '-' && rootSequence[i].first[j].second[k] != 'x') {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = rootSequence[i].first[j].second[k];
                            // if(rootCtr == 240) {
                            //     std::cout << rootSequence[i].first[j].first << " " << i << " " << j << " " << k << std::endl;
                            // }
                            isGapCoordinate[std::make_tuple(i,j,k)] = false;
                            rootCtr++;
                        } else {
                            rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                            isGapCoordinate[std::make_tuple(i,j,k)] = true;
                        }
                    }
                }
            }
        } else {
            rootPresentBlocks[i] = false;
            if(rootBlockStrand[i].first) {
                for(size_t j = 0; j < rootSequence[i].first.size(); j++) {
                    for(size_t k = 0; k < rootSequence[i].first[j].second.size(); k++) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                }
            } else {
                for(size_t j = rootSequence[i].first.size() - 1; j + 1 > 0; j--) {
                    panMATCoordinateToGlobal[std::make_tuple(i,j,-1)] = rootCtr;
                    rootCurrentCharacter[std::make_tuple(i,j,-1)] = '-';
                    isGapCoordinate[std::make_tuple(i,j,-1)] = false;
                    for(size_t k = rootSequence[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        panMATCoordinateToGlobal[std::make_tuple(i,j,k)] = rootCtr;
                        rootCurrentCharacter[std::make_tuple(i,j,k)] = '-';
                        isGapCoordinate[std::make_tuple(i,j,k)] = false;
                    }
                }
            }
        }
    }

    tbb::concurrent_map< std::tuple< std::string, int, int, int >, char > seqChar;

    tbb::parallel_for_each(allNodes, [&](auto u) {
        sequence_t st;
        blockExists_t bt;
        blockStrand_t bst;
        getSequenceFromReference(st, bt, bst, u.first);

        for(size_t i = 0; i < st.size(); i++) {
            if(bst[i].first) {
                for(size_t j = 0; j < st[i].first.size(); j++) {
                    for(size_t k = 0; k < st[i].first[j].second.size(); k++) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x' && bt[i].first) {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = st[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                }
            } else {
                for(size_t j = st[i].first.size() - 1; j + 1 > 0; j--) {
                    if(st[i].first[j].first != '-' && st[i].first[j].first != 'x') {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = rootSequence[i].first[j].first;
                    } else {
                        seqChar[std::make_tuple(u.first,i,j,-1)] = '-';
                    }
                    for(size_t k = st[i].first[j].second.size() - 1; k + 1 > 0; k--) {
                        if(st[i].first[j].second[k] != '-' && st[i].first[j].second[k] != 'x') {
                            seqChar[std::make_tuple(u.first,i,j,k)] = st[i].first[j].second[k];
                        } else {
                            seqChar[std::make_tuple(u.first,i,j,k)] = '-';
                        }
                    }
                }
            }
        }
    });

    nodeMutations[root->identifier];

    // Compute mutations for each of the other sequences
    tbb::parallel_for_each(allNodes, [&](auto u) {
        if(u.first == root->identifier) {
            return;
        }

        tbb::concurrent_map< std::tuple< int, int, int >, char > currentCharacter = rootCurrentCharacter;
        tbb::concurrent_unordered_map< size_t, bool > presentBlocks = rootPresentBlocks;

        nodeMutations[u.first];

        Node* it = u.second;
        std::vector< panmanUtils::Node* > path;

        while(it != root) {
            path.push_back(it);
            it = it->parent;
        }

        std::vector< std::pair< size_t, std::tuple< char, size_t, char, char, bool > > > currentNodeMutations;

        for(auto node = path.rbegin(); node != path.rend(); node++) {
            for(auto mutation: (*node)->blockMutation) {
                int32_t primaryBlockId = mutation.primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                bool type = mutation.blockMutInfo;
                bool inversion = mutation.inversion;
                if(type == 1) {
                    // insertion
                    presentBlocks[primaryBlockId] = true;
                    if(inversion) {
                        std::cout << "INVERTED BLOCK FOUND" << std::endl;
                    }
                } else {
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        std::cout << "INVERSION FOUND" << std::endl;
                    } else {
                        // Actually a deletion
                        presentBlocks[primaryBlockId] = false;
                    }
                }
            }
        }

        for(auto node = path.rend()-1; node != path.rend(); node++) {
            for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {
                int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
                if(rootPresentBlocks.find(primaryBlockId) == rootPresentBlocks.end()) {
                    continue;
                }

                int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
                int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
                uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';

                if(type < 3) {
                    int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                    if(type == panmanUtils::NucMutationType::NS) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)];
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition+j)];
                                    if(oldVal == '-' || oldVal == 'x') {
                                        continue;
                                    }
                                    // if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(presentBlocks[primaryBlockId]) {
                                    char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition+j, -1)];
                                    // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)];
                                    // if(u.first == "Denmark/DCGC-504971/2022|OX187739.1|2022-04-28" && panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)] == 185) {
                                    //     std::cout << oldVal << " " << newVal << " " << primaryBlockId << " " << nucPosition+j << " " << -1 << std::endl;
                                    // }
                                    if(oldVal == '-' || oldVal == 'x') {
                                        // std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                                        continue;
                                    }
                                    // if(oldVal == newVal) {
                                    //     std::cout << primaryBlockId << " " << nucPosition << " " << nucGapPosition+j << " " << oldVal << " " << newVal << std::endl;
                                    // }
                                    // if(node == path.rend()-1)
                                    currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, -1)])));
                                    // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, newVal));
                                }
                                currentCharacter[std::make_tuple(primaryBlockId, nucPosition + j, -1)] = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                // if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                // if(node == path.rend()-1)
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('I', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], '-', newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));

                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition + j)];
                                // if(node == path.rend()-1){
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition + j)])));
                                // }
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition + j, nucGapPosition)];
                                // if(node == path.rend()-1){
                                currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('D', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition + j, -1)], oldVal, '-', isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition + j, -1)])));
                                // }
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(nucGapPosition != -1) {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, nucGapPosition)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            // if(node == path.rend()-1)
                            currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)] = newVal;
                    } else {
                        if(presentBlocks[primaryBlockId]) {
                            char oldVal = seqChar[std::make_tuple((*node)->parent->identifier,primaryBlockId, nucPosition, -1)];
                            // char oldVal = currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)];
                            if(oldVal == '-' || oldVal == 'x') {
                                std::cout << "NOT ACTUALLY A SUBSTITUTION" << std::endl;
                            }
                            // if(node == path.rend()-1)
                            currentNodeMutations.push_back(std::make_pair(primaryBlockId, std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal, isGapCoordinate[std::make_tuple(primaryBlockId, nucPosition, -1)])));
                            // nodeMutations[u.first].push_back(std::make_tuple('S', panMATCoordinateToGlobal[std::make_tuple(primaryBlockId, nucPosition, -1)], oldVal, newVal));
                        }
                        currentCharacter[std::make_tuple(primaryBlockId, nucPosition, -1)] = newVal;
                    }
                }
            }
        }

        for(auto mut: currentNodeMutations) {
            if(presentBlocks.find(mut.first) != presentBlocks.end()) {
                nodeMutations[u.first].push_back(mut.second);
            }
        }
    });

    for(auto& u: nodeMutations) {
        // print all substitutions first
        fout << "Substitutions:\t";
        fout << u.first << '\t';
        for(auto v: u.second) {
            if(std::get<0>(v) == 'S') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<2>(v) << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Insertions:\t";
        fout << u.first << '\t';
        // print insertions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'I') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<3>(v);
            }
        }
        fout << '\n';

        fout << "Deletions:\t";
        fout << u.first << '\t';
        // print deletions
        for(auto v: u.second) {
            if(std::get<0>(v) == 'D') {
                fout << " > " << (std::get<4>(v) ? "g" : "") << std::get<1>(v)+1 << std::get<2>(v);
            }
        }
        fout << '\n';
    }

}

void panmanUtils::Tree::getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists, blockStrand_t& blockStrand, std::string reference, bool rotateSequence, int* rotIndex) {
    Node* referenceNode = nullptr;

    for(auto u: allNodes) {
        // printf("%s\n",u.first);
        // std::cerr << u.first << std::endl;
        if(u.first == reference) {
            referenceNode = u.second;
            break;
        }
    }

    // printf(reference)

    if(referenceNode == nullptr) {
        std::cerr << "Error: Reference sequence with matching name not found: " << reference << std::endl;
        return;
    }

    std::vector< panmanUtils::Node* > path;
    Node* it = referenceNode;

    while(it != root) {
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(root);

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    sequence.resize(blocks.size() + 1);
    blockExists.resize(blocks.size() + 1, {false, {}});
    blockStrand.resize(blocks.size() + 1, {true, {}});

    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        blockStrand[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], true);
    }

    int32_t maxBlockId = 0;

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
                const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

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

    // Assigning nucleotide gaps
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

    // Get all blocks on the path
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(auto mutation: (*node)->blockMutation) {
            int primaryBlockId = mutation.primaryBlockId;
            int secondaryBlockId = mutation.secondaryBlockId;
            int type = (mutation.blockMutInfo);
            bool inversion = mutation.inversion;

            if(type == panmanUtils::BlockMutationType::BI) {
                if(secondaryBlockId != -1) {
                    blockExists[primaryBlockId].second[secondaryBlockId] = true;

                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
                } else {
                    blockExists[primaryBlockId].first = true;

                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId].first = !inversion;
                }
            } else {
                if(inversion) {
                    // This is not actually a deletion but an inversion
                    if(secondaryBlockId != -1) {
                        blockStrand[primaryBlockId].second[secondaryBlockId] = !blockStrand[primaryBlockId].second[secondaryBlockId];
                    } else {
                        blockStrand[primaryBlockId].first = !blockStrand[primaryBlockId].first;
                    }
                } else {
                    // Actually a deletion
                    if(secondaryBlockId != -1) {
                        blockExists[primaryBlockId].second[secondaryBlockId] = false;
                        blockStrand[primaryBlockId].second[secondaryBlockId] = true;
                    } else {
                        blockExists[primaryBlockId].first = false;
                        blockStrand[primaryBlockId].first = true;
                    }
                }
            }
        }
    }

    // Apply nucleotide mutations
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {

            int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = (*node)->nucMutation[i].secondaryBlockId;

            if(secondaryBlockId != -1) {
                if(!blockExists[primaryBlockId].second[secondaryBlockId]) {
                    continue;
                }
            } else {
                if(!blockExists[primaryBlockId].first) {
                    continue;
                }
            }

            int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
            int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
            uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3) {

                int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                if(type == panmanUtils::NucMutationType::NS) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NI) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::ND) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            }
                        }
                    }
                }
            } else {
                if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPI) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPD) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    if(rotateSequence) {
        if(rotationIndexes.find(reference) != rotationIndexes.end() && rotationIndexes[reference] != 0) {
            int ctr = -1, rotInd = 0;
            for(size_t i = 0; i < blockExists.size(); i++) {
                if(blockExists[i].first) {
                    ctr++;
                }
                if(ctr == rotationIndexes[reference]) {
                    rotInd = i;
                    break;
                }
            }
            if(rotIndex != nullptr) {
                *rotIndex = rotInd;
            }
            rotate(sequence.begin(), sequence.begin() + rotInd, sequence.end());
            rotate(blockExists.begin(), blockExists.begin() + rotInd, blockExists.end());
            rotate(blockStrand.begin(), blockStrand.begin() + rotInd, blockStrand.end());
        }

        if(sequenceInverted.find(reference) != sequenceInverted.end() && sequenceInverted[reference]) {
            reverse(sequence.begin(), sequence.end());
            reverse(blockExists.begin(), blockExists.end());
            reverse(blockStrand.begin(), blockStrand.end());
        }
    }
}

std::string panmanUtils::Tree::getStringFromReference(std::string reference, bool aligned,  bool incorporateInversions) {

    Node* referenceNode = nullptr;

    for(auto u: allNodes) {
        // printf("%s\n", u.first);
        // std::cerr << u.first << std::endl;
        if(u.first == reference) {
            referenceNode = u.second;
            break;
        }
    }

    if(referenceNode == nullptr) {
        return "Error: Reference sequence with matching name not found!";
    }

    std::vector< panmanUtils::Node* > path;
    Node* it = referenceNode;

    while(it != root) {
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(root);

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

    // Create block consensus sequences
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
                const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

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

    // Assigning nucleotide gaps
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

    // Get all blocks on the path
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(auto mutation: (*node)->blockMutation) {
            int primaryBlockId = mutation.primaryBlockId;
            int secondaryBlockId = mutation.secondaryBlockId;
            int type = (mutation.blockMutInfo);
            bool inversion = mutation.inversion;

            if(type == panmanUtils::BlockMutationType::BI) {
                if(secondaryBlockId != -1) {
                    blockExists[primaryBlockId].second[secondaryBlockId] = true;

                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
                } else {
                    blockExists[primaryBlockId].first = true;

                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId].first = !inversion;
                }
            } else {
                if(inversion) {
                    // This is not actually a deletion but an inversion
                    if(secondaryBlockId != -1) {
                        blockStrand[primaryBlockId].second[secondaryBlockId] = !blockStrand[primaryBlockId].second[secondaryBlockId];
                    } else {
                        blockStrand[primaryBlockId].first = !blockStrand[primaryBlockId].first;
                    }
                } else {
                    // Actually a deletion
                    if(secondaryBlockId != -1) {
                        blockExists[primaryBlockId].second[secondaryBlockId] = false;
                        blockStrand[primaryBlockId].second[secondaryBlockId] = true;
                    } else {
                        blockExists[primaryBlockId].first = false;
                        blockStrand[primaryBlockId].first = true;
                    }
                }
            }
        }
    }

    // Apply nucleotide mutations
    for(auto node = path.rbegin(); node != path.rend(); node++) {

        for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {

            int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = (*node)->nucMutation[i].secondaryBlockId;

            if(secondaryBlockId != -1) {
                if(!blockExists[primaryBlockId].second[secondaryBlockId]) {
                    continue;
                }
            } else {
                if(!blockExists[primaryBlockId].first) {
                    continue;
                }
            }

            int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
            int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
            uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3) {

                int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                if(type == panmanUtils::NucMutationType::NS) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NI) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::ND) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            }

                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            }
                        }
                    }
                }
            } else {
                if(type == panmanUtils::NucMutationType::NSNPS) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPI) {
                    newVal = panmanUtils::getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == panmanUtils::NucMutationType::NSNPD) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    if(!aligned && rotationIndexes.find(reference) != rotationIndexes.end() && rotationIndexes[reference] != 0) {
        int ctr = -1, rotInd = 0;
        for(size_t i = 0; i < blockExists.size(); i++) {
            if(blockExists[i].first) {
                ctr++;
            }
            if(ctr == rotationIndexes[reference]) {
                rotInd = i;
                break;
            }
        }
        rotate(sequence.begin(), sequence.begin() + rotInd, sequence.end());
        rotate(blockExists.begin(), blockExists.begin() + rotInd, blockExists.end());
        rotate(blockStrand.begin(), blockStrand.begin() + rotInd, blockStrand.end());
    }

    if(sequenceInverted.find(reference) != sequenceInverted.end() && sequenceInverted[reference]) {
        reverse(sequence.begin(), sequence.end());
        reverse(blockExists.begin(), blockExists.end());
        reverse(blockStrand.begin(), blockStrand.end());
    }

    std::string sequenceString;
    for(size_t i = 0; i < sequence.size(); i++) {
        // Iterate through gap blocks - CURRENTLY NOT BEING USED
        for(size_t j = 0; j < sequence[i].second.size(); j++) {
            if(blockExists[i].second[j]) {
                if(blockStrand[i].second[j]) {
                    // If forward strand
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                            if(sequence[i].second[j][k].second[w] == 'x' || sequence[i].second[j][k].second[w] == '-' ) {
                                if(aligned) {
                                    sequenceString+='-';
                                }
                            } else {
                                sequenceString += sequence[i].second[j][k].second[w];
                            }
                        }
                        if(sequence[i].second[j][k].first == 'x' || sequence[i].second[j][k].first == '-') {
                            if(aligned) {
                                sequenceString+='-';
                            }
                        } else {
                            sequenceString += sequence[i].second[j][k].first;
                        }
                    }
                } else {
                    for(size_t k = sequence[i].second[j].size()-1; k + 1 > 0; k--) {
                        // If reverse strand
                        if(sequence[i].second[j][k].first == 'x' || sequence[i].second[j][k].first == '-') {
                            if(aligned) {
                                sequenceString+='-';
                            }
                        } else {
                            sequenceString += sequence[i].second[j][k].first;
                        }
                        for(size_t w = sequence[i].second[j][k].second.size() - 1; w + 1 > 0; w--) {
                            if(sequence[i].second[j][k].second[w] == 'x' || sequence[i].second[j][k].second[w] == '-') {
                                if(aligned) {
                                    sequenceString+='-';
                                }
                            } else {
                                sequenceString += sequence[i].second[j][k].second[w];
                            }
                        }
                    }
                }
            } else {
                if(aligned) {
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                            sequenceString+='-';
                        }
                        sequenceString+='-';
                    }
                }
            }
        }

        // Main block
        if(blockExists[i].first) {
            if(blockStrand[i].first) {
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                        if(sequence[i].first[j].second[k] == 'x' || sequence[i].first[j].second[k] == '-' ) {
                            // This shouldn't be possible but I'm still keeping it since it doesn't hurt
                            if(aligned) {
                                sequenceString += '-';
                            }
                        } else {
                            sequenceString += sequence[i].first[j].second[k];
                        }
                    }
                    if(sequence[i].first[j].first == 'x' || sequence[i].first[j].first == '-') {
                        if(aligned) {
                            sequenceString += '-';
                        }
                    } else {
                        sequenceString += sequence[i].first[j].first;
                    }
                }
            } else {
                // If reverse strand
                for(size_t j = sequence[i].first.size()-1; j + 1 > 0; j--) {
                    if(sequence[i].first[j].first == 'x' || sequence[i].first[j].first == '-') {
                        if(aligned) {
                            sequenceString += '-';
                        }
                    } else {
                        sequenceString += getComplementCharacter(sequence[i].first[j].first);
                    }
                    for(size_t k = sequence[i].first[j].second.size() - 1; k+1 > 0; k--) {
                        if(sequence[i].first[j].second[k] == 'x' || sequence[i].first[j].second[k] == '-' ) {
                            // This shouldn't be possible but I'm still keeping it since it doesn't hurt
                            if(aligned) {
                                sequenceString += '-';
                            }
                        } else {
                            sequenceString += getComplementCharacter(sequence[i].first[j].second[k]);
                        }
                    }
                }
            }
        } else {
            if(aligned) {
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                        sequenceString+='-';
                    }
                    sequenceString+='-';
                }
            }
        }
    }

    int offset = 0;
    if(!aligned && circularSequences.find(reference) != circularSequences.end()) {
        offset = circularSequences[reference];
    }
    if(offset == 0) {
        return sequenceString;
    } else {
        return sequenceString.substr(offset) + sequenceString.substr(0,offset);
    }

}

std::string panmanUtils::stripString(std::string s) {
    while(s.length() && s[s.length() - 1] == ' ') {
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++) {
        if(s[i] != ' ') {
            return s.substr(i);
        }
    }
    return s;
}

std::string panmanUtils::stripGaps(const std::string sequenceString) {
    std::string result;
    for(auto u: sequenceString) {
        if(u != '-' && u != 'x') {
            result+=u;
        }
    }
    return result;
}

bool panmanUtils::Tree::verifyVCFFile(std::ifstream& fin) {

    for(auto u: allNodes) {
        if(u.second->children.size() == 0) {
            std::cout << u.first << std::endl;
            fin.clear();
            fin.seekg(0);
            if(getSequenceFromVCF(u.first, fin) != getStringFromReference(u.first, false)) {
                return false;
            }
            std::cout << u.first << std::endl;
        }
    }

    return true;
}

void panmanUtils::Tree::vcfToFASTA(std::ifstream& fin, std::ofstream& fout) {
    for(auto u: allNodes) {
        if(u.second->children.size() == 0) {
            fin.clear();
            fin.seekg(0);
            std::string sequenceString = getSequenceFromVCF(u.first, fin);
            fout << '>' << u.first << '\n';
            for(size_t i = 0; i < sequenceString.size(); i+=70) {
                fout << sequenceString.substr(i,std::min(70, (int)sequenceString.length() - (int)i)) << '\n';
            }
        }
    }
}

std::string panmanUtils::Tree::getSequenceFromVCF(std::string sequenceId, std::ifstream& fin) {
    std::string line;

    // get reference line
    for(int i = 0; i < 4; i++) {
        std::getline(fin, line);
    }

    if(line.substr(0,12) != "##reference=") {
        std::cout << "Incorrect line format: " << line << std::endl;
        return "";
    }

    std::string referenceSequenceId = line.substr(12);

    std::string referenceSequence = stripGaps(getStringFromReference(referenceSequenceId));


    if(sequenceId == referenceSequenceId) {
        return referenceSequence;
    }

    // column headers
    std::getline(fin, line);

    std::vector< std::string > columnWords;
    std::string word;

    for(size_t i = 0; i < line.size(); i++) {
        if(line[i] != ' ' && line[i]!='\t') {
            word += line[i];
        } else {
            if(word.length()) {
                columnWords.push_back(word);
                word = "";
            }
        }
    }
    if(word.length()) {
        columnWords.push_back(word);
    }

    int sequenceIndex = -1;

    for(size_t i = 9; i < columnWords.size(); i++) {
        if(columnWords[i] == sequenceId) {
            sequenceIndex = i;
            break;
        }
    }

    if(sequenceIndex == -1) {
        std::cout << "sequence not found! "<< sequenceId << std::endl;
        return "";
    }

    // To account for insertions
    std::vector< std::pair< char, std::vector< char > > > alteredSequence;
    for(auto u: referenceSequence) {
        alteredSequence.push_back({u, {}});
    }
    alteredSequence.push_back({'-', {}});

    while(getline(fin, line)) {
        std::vector< std::string > words;
        std::string word;

        for(size_t i = 0; i < line.size(); i++) {
            if(line[i] != ' ' && line[i] != '\t') {
                word += line[i];
            } else {
                if(word.length()) {
                    words.push_back(word);
                    word = "";
                }
            }
        }

        if(word.length()) {
            words.push_back(word);
            word="";
        }

        int choice = std::stoll(words[sequenceIndex]);
        if(choice == 0) {
            continue;
        }

        choice--;

        int position = std::stoll(words[1]);

        std::string ref = words[3];

        std::string altStrings = words[4];

        std::string currentAlt;
        std::vector< std::string > altChoices;

        for(auto u: altStrings) {
            if(u != ',') {
                currentAlt += u;
            } else {
                if(currentAlt.length()) {
                    altChoices.push_back(currentAlt);
                    currentAlt = "";
                }
            }
        }

        if(currentAlt.length()) {
            altChoices.push_back(currentAlt);
            currentAlt = "";
        }

        std::string alt = altChoices[choice];

        if(ref != ".") {
            int len = ref.length();
            for(int i = position; i < position + len; i++) {
                alteredSequence[i].first = '-';
            }
        }


        if(alt != ".") {
            if(alt.length() && alteredSequence[position].second.size()) {
                std::cout << "VCF Error: alternate sequence already exists at position " << position <<"!" << std::endl;
                std::cout << sequenceId << " " << referenceSequenceId << std::endl;
            }
            for(size_t i = 0; i < alt.length(); i++) {
                alteredSequence[position].second.push_back(alt[i]);
            }
        }

    }

    std::string finalSequence;
    for(size_t i = 0; i < alteredSequence.size(); i++) {
        for(size_t j = 0; j < alteredSequence[i].second.size(); j++) {
            if(alteredSequence[i].second[j] != '-') {
                finalSequence += alteredSequence[i].second[j];
            }
        }
        if(alteredSequence[i].first != '-') {
            finalSequence += alteredSequence[i].first;
        }

    }

    // std::string alteredSequenceOriginal = stripGaps(getStringFromReference(sequenceId));

    return finalSequence;

}

void panmanUtils::Tree::printFASTAParallel(std::ofstream& fout, bool aligned) {

    std::mutex fastaMutex;
    size_t lineSize = 70;

    tbb::parallel_for_each(allNodes, [&](auto n) {
        if(n.second->children.size() == 0) {
            std::string sequence;
            sequence = getStringFromReference(n.first, aligned);

            fastaMutex.lock();
            fout << '>' << n.first << '\n';
            for(size_t i = 0; i < sequence.size(); i+=lineSize) {
                fout << sequence.substr(i, std::min(lineSize, sequence.size() - i)) << '\n';
            }
            fastaMutex.unlock();
        }

    });
}


void panmanUtils::Tree::printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout) {
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
            stringSplit(separatedLine[2], ',', paths[separatedLine[1]]);
        }
    }
    for(auto p: paths) {
        fout << ">" << p.first << "\n";
        std::string sequence;
        for(auto s: p.second) {
            char strand = s[s.length()-1];
            s.pop_back();
            if(strand == '+') {
                sequence += nodes[s];
            } else {
                for (std::string::reverse_iterator rit=nodes[s].rbegin(); rit!=nodes[s].rend(); ++rit) {
                    sequence += getComplementCharacter(*rit);
                }
            }
        }
        for(size_t i = 0; i < sequence.size(); i+=70) {
            fout << sequence.substr(i,std::min((size_t)70, sequence.size() - i)) << '\n';
        }
    }
}

int32_t panmanUtils::Tree::getUnalignedGlobalCoordinate(int32_t primaryBlockId,
        int32_t secondaryBlockId, int32_t pos, int32_t gapPos, const sequence_t& sequence,
        const blockExists_t& blockExists, const blockStrand_t& blockStrand, int circularOffset) {

    // std::cout << "P " << primaryBlockId << " " << secondaryBlockId << " " << pos << " " << gapPos << " " << sequence[primaryBlockId].first[pos].first << std::endl;

    int ctr = 0;
    int ans = -1;
    int len = 0;
    for(size_t i = 0; i < blockExists.size(); i++) {
        if(!blockExists[i].first) {
            continue;
        }
        if(blockStrand[i].first) {
            for(size_t k = 0; k < sequence[i].first.size(); k++) {
                for(size_t w = 0; w < sequence[i].first[k].second.size(); w++) {
                    if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x') {
                        if((int)i == primaryBlockId && secondaryBlockId == -1 && (int)k == pos && (int)w == gapPos) {
                            ans = ctr;
                        }
                        if(ans==-1) {
                            ctr++;
                        }
                        len++;
                    }
                }
                if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                    if((int)i == primaryBlockId && secondaryBlockId == -1 && (int)k == pos && gapPos == -1) {
                        ans = ctr;
                    }
                    if(ans==-1) {
                        ctr++;
                    }
                    len++;
                }
            }
        } else {
            for(size_t k = sequence[i].first.size() - 1; k + 1 > 0; k--) {
                if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                    if((int)i == primaryBlockId && secondaryBlockId == -1 && (int)k == pos && gapPos == -1) {
                        ans = ctr;
                    }
                    if(ans==-1) {
                        ctr++;
                    }
                    len++;
                }
                for(size_t w = sequence[i].first[k].second.size() - 1; w+1 > 0; w--) {
                    if(sequence[i].first[k].second[w] != '-'
                            && sequence[i].first[k].second[w] != 'x') {
                        if((int)i == primaryBlockId && secondaryBlockId == -1 && (int)k == pos
                                && (int)w == gapPos) {
                            ans = ctr;
                        }
                        if(ans==-1) {
                            ctr++;
                        }
                        len++;
                    }
                }
            }
        }
    }

    // std::cout << "ANS: " << ans << " " << circularOffset << std::endl;
    ans -= circularOffset;
    if(ans < 0) {
        ans += len;
    }
    return ans;
}

std::tuple< int, int, int, int > panmanUtils::Tree::globalCoordinateToBlockCoordinate(
    int64_t globalCoordinate, const sequence_t& sequence, const blockExists_t& blockExists,
    const blockStrand_t& blockStrand, int64_t circularOffset) {

    // Computing length of sequence
    int len = 0;
    for(size_t i = 0; i < blockExists.size(); i++) {
        if(!blockExists[i].first) {
            continue;
        }
        for(size_t k = 0; k < sequence[i].first.size(); k++) {
            for(size_t w = 0; w < sequence[i].first[k].second.size(); w++) {
                if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x') {
                    len++;
                }
            }
            if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                len++;
            }
        }
    }


    // Adjusting for circular offset
    if(circularOffset + globalCoordinate < len) {
        globalCoordinate += circularOffset;
    } else {
        globalCoordinate = globalCoordinate + circularOffset - len;
    }

    int ctr = 0;
    for(size_t i = 0; i < blockExists.size(); i++) {
        if(!blockExists[i].first) {
            continue;
        }
        if(blockStrand[i].first) {
            for(size_t k = 0; k < sequence[i].first.size(); k++) {
                for(size_t w = 0; w < sequence[i].first[k].second.size(); w++) {
                    if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x') {
                        if(ctr == globalCoordinate) {
                            return std::make_tuple(i, -1, k, w);
                        }
                        ctr++;
                    }
                }
                if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                    if(ctr == globalCoordinate) {
                        return std::make_tuple(i, -1, k, -1);
                    }
                    ctr++;
                }
            }
        } else {
            for(size_t k = sequence[i].first.size() - 1; k + 1 > 0; k--) {
                if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                    if(ctr == globalCoordinate) {
                        return std::make_tuple(i, -1, k, -1);
                    }
                    ctr++;
                }
                for(size_t w = sequence[i].first[k].second.size() - 1; w + 1 > 0; w--) {
                    if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x') {
                        if(ctr == globalCoordinate) {
                            return std::make_tuple(i, -1, k, w);
                        }
                        ctr++;
                    }
                }
            }
        }
    }
    return std::make_tuple(-1,-1,-1,-1);
}

void panmanUtils::Tree::adjustLevels(Node* node) {
    if(node->parent == nullptr) {
        node->level = 1;
    } else {
        node->level = node->parent->level + 1;
    }
    for(auto u: node->children) {
        adjustLevels(u);
    }
}

panmanUtils::Node* panmanUtils::Tree::transformHelper(Node* node) {
    if(node == root) {
        if(node->children.size() > 1) {
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
    for(size_t i = 0; i < par->children.size(); i++) {
        if(par->children[i]->identifier == node->identifier) {
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

void panmanUtils::Tree::transform(Node* node) {
    Node* par = node->parent;
    if(par == nullptr) {
        // already root
        return;
    }
    if(par == root) {
        // Parent already root. The root will contain the same sequence as the node
        node->branchLength = 0;
        return;
    }

    // remove node from parent's children
    for(size_t i = 0; i < par->children.size(); i++) {
        if(par->children[i]->identifier == node->identifier) {
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

panmanUtils::Tree::Tree(Node* newRoot, const std::vector< Block >& b,
                        const std::vector< GapList >& g, std::unordered_map< std::string, int >& cs,
                        std::unordered_map< std::string, int >& ri,
                        std::unordered_map< std::string, bool >& si,
                        const BlockGapList& bgl) {

    std::set< std::string > nodeIds;
    getNodesRootedAt(nodeIds, newRoot);

    root = newRoot;
    blocks = b;
    gaps = g;

    // Adjusting allNodes and other node data
    std::queue< Node* > q;
    q.push(root);
    while(!q.empty()) {
        Node* current = q.front();
        q.pop();
        allNodes[current->identifier] = current;
        if(cs.find(current->identifier) != cs.end()) {
            circularSequences[current->identifier] = cs[current->identifier];
        }
        if(ri.find(current->identifier) != ri.end()) {
            rotationIndexes[current->identifier] = ri[current->identifier];
        }
        if(si.find(current->identifier) != si.end()) {
            sequenceInverted[current->identifier] = si[current->identifier];
        }
        for(auto child: current->children) {
            q.push(child);
        }
    }

    blockGaps = bgl;
}

std::pair< panmanUtils::Tree, panmanUtils::Tree > panmanUtils::Tree::splitByComplexMutations(const std::string& nodeId3) {

    Node* newRoot = allNodes[nodeId3];

    // getting all prior mutations
    // Block mutations
    std::vector< panmanUtils::Node* > path;
    panmanUtils::Node* currentNode = newRoot;
    while(currentNode != nullptr) {
        path.push_back(currentNode);
        currentNode = currentNode->parent;
    }

    // For block mutations, we cancel out irrelevant mutations
    std::map< std::pair<int, int>, std::pair< panmanUtils::BlockMutationType, bool > > bidMutations;
    std::vector< panmanUtils::BlockMut > newBlockMutation;
    std::vector< panmanUtils::NucMut > newNucMutation;

    for(auto itr = path.rbegin(); itr!= path.rend(); itr++) {
        for(auto mutation: (*itr)->blockMutation) {
            int primaryBlockId = mutation.primaryBlockId;
            int secondaryBlockId = mutation.secondaryBlockId;
            bool type = (mutation.blockMutInfo);
            bool inversion = (mutation.inversion);

            if(type == panmanUtils::BlockMutationType::BI) {
                bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair( panmanUtils::BlockMutationType::BI, inversion );
            } else {
                if(bidMutations.find(std::make_pair(primaryBlockId, secondaryBlockId)) != bidMutations.end()) {
                    if(bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].first == panmanUtils::BlockMutationType::BI) {
                        // If it was insertion earlier
                        if(inversion) {
                            // This means that the new mutation is an inversion. So, inverted the strand of the inserted block
                            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = !bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second;
                        } else {
                            // Actually a deletion. So insertion and deletion cancel out
                            bidMutations.erase(std::make_pair(primaryBlockId, secondaryBlockId));
                        }
                    } else {
                        // If previous mutation was an inversion
                        if(!inversion) {
                            // Actually a deletion. Remove inversion mutation and put deletion instead
                            bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)].second = false;
                        }
                        // deletion followed by inversion doesn't make sense
                    }
                } else {
                    bidMutations[std::make_pair(primaryBlockId, secondaryBlockId)] = std::make_pair(panmanUtils::BlockMutationType::BD, inversion);
                }
            }
        }
        for(auto mutation: (*itr)->nucMutation) {
            newNucMutation.push_back( mutation );
        }
    }

    for(auto mutation: bidMutations) {
        if(mutation.second.first == panmanUtils::BlockMutationType::BI) {
            panmanUtils::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = panmanUtils::BlockMutationType::BI;
            tempBlockMut.inversion = mutation.second.second;
            newBlockMutation.push_back( tempBlockMut );
        } else {
            panmanUtils::BlockMut tempBlockMut;
            tempBlockMut.primaryBlockId = mutation.first.first;
            tempBlockMut.secondaryBlockId = mutation.first.second;
            tempBlockMut.blockMutInfo = panmanUtils::BlockMutationType::BD;
            tempBlockMut.inversion = mutation.second.second;
            newBlockMutation.push_back( tempBlockMut );
        }
    }

    // Finally assigning Block Mutations
    newRoot->blockMutation = newBlockMutation;

    // Assigning Nuc Mutations
    newRoot->nucMutation = consolidateNucMutations(newNucMutation);

    // Adjusting parent and child pointers
    if(newRoot->parent != nullptr) {
        for(size_t i = 0; i < newRoot->parent->children.size(); i++) {
            if(newRoot->parent->children[i]->identifier == nodeId3) {
                newRoot->parent->children.erase(newRoot->parent->children.begin() + i);
                break;
            }
        }
    }
    newRoot->parent = nullptr;

    // Creating a whole new tree out of the child
    Tree childTree(newRoot, blocks, gaps, circularSequences, rotationIndexes, sequenceInverted,
                   blockGaps);

    // Removing child tree's nodes from current tree
    std::queue< Node* > q;
    q.push(newRoot);
    while(!q.empty()) {
        Node* current = q.front();
        q.pop();
        allNodes.erase(current->identifier);
        circularSequences.erase(current->identifier);
        rotationIndexes.erase(current->identifier);
        sequenceInverted.erase(current->identifier);
        for(auto child: current->children) {
            q.push(child);
        }
    }

    return std::make_pair(*this, childTree);

}

panmanUtils::GfaGraph::GfaGraph(const std::vector< std::string >& pathNames, const std::vector< std::vector< std::pair< std::string, bool > > >& sequences, std::map< std::string, std::string >& nodes) {
    pathIds = pathNames;

    // New integer Node ID to old string Node ID
    std::vector< std::string > nodeIdToString;

    // Auto increment ID to assign to nodes
    numNodes = 0;

    intSequences.resize(sequences.size());
    strandPaths.resize(sequences.size());

    std::unordered_map<int,std::string> intToString; // Locally stored -> Later on mapped to intIdToStringId
    std::vector<std::string> consensus = {};
    std::vector<int> intSequenceConsensus= {};

    for(size_t seqCount = 0; seqCount < sequences.size(); seqCount++) {
        if (seqCount == 0) {
            // Load first sequence path
            for(const auto& block: sequences[seqCount]) {
                consensus.push_back(block.first);
                // sample_base.push_back(block);
                intToString[numNodes] = block.first;
                intSequences[0].push_back(numNodes);
                strandPaths[seqCount].push_back(block.second);
                intSequenceConsensus.push_back(numNodes);
                numNodes++;
            }
        } else {
            std::vector<int> intSequenceSample;
            std::vector<int> intSequenceConsensusNew;
            std::vector<std::string> sample;
            std::vector<std::string> consensusNew;

            for(const auto& block: sequences[seqCount]) {
                sample.push_back(block.first);
                strandPaths[seqCount].push_back(block.second);
            }

            chain_align (consensus,
                         sample,
                         intSequenceConsensus,
                         intSequenceSample,
                         numNodes,
                         consensusNew,
                         intSequenceConsensusNew,
                         intToString);
            for (auto &b: intSequenceSample) {
                intSequences[seqCount].push_back(b);
            }
            consensus.clear();
            intSequenceConsensus.clear();
            for (auto &b: consensusNew) {
                consensus.push_back(b);
            }

            for (auto &b: intSequenceConsensusNew) {
                intSequenceConsensus.push_back(b);
            }
        }
        std::cout << seqCount << " " << intSequenceConsensus.size() << endl;
    }

    // re-assigning IDs in fixed order
    size_t reorder = 0;
    std::unordered_map<int, int> orderMap = {};
    for (auto &i: intSequenceConsensus) {
        orderMap[i] = reorder;
        nodeIdToString.push_back(intToString[i]);
        topoSortedIntSequences.push_back(reorder);
        reorder++;
    }

    for (auto &m: intSequences) {
        for (auto &s: m) {
            s = orderMap[s];
        }
    }

    for(size_t i = 0; i < reorder; i++) {
        intNodeToSequence.push_back(nodes[nodeIdToString[i]]);
    }
}

std::vector< std::vector< int64_t > > panmanUtils::GfaGraph::getAlignedSequences(const std::vector< size_t >& topoArray) {
    std::vector< std::vector< int64_t > > newSequences;

    for(auto sequence: intSequences) {
        size_t p1 = 0, p2 = 0;
        std::vector< int64_t > newSequence;
        while(p1 < topoArray.size() && p2 < sequence.size()) {
            if(topoArray[p1] == sequence[p2]) {
                newSequence.push_back(sequence[p2]);
                p2++;
            } else {
                newSequence.push_back(-1);
            }
            p1++;
        }
        // Adding gaps towards the end
        while(newSequence.size() < topoArray.size()) {
            newSequence.push_back(-1);
        }
        newSequences.push_back(newSequence);
    }
    return newSequences;
}

std::vector< std::vector< int > > panmanUtils::GfaGraph::getAlignedStrandSequences(
    const std::vector< size_t >& topoArray) {
    std::vector< std::vector< int > > alignedStrandSequences;

    int currentSequence = 0;
    for(auto sequence: intSequences) {
        size_t p1 = 0, p2 = 0;
        std::vector< int > strandSequence;
        while(p1 < topoArray.size() && p2 < sequence.size()) {
            if(topoArray[p1] == sequence[p2]) {
                strandSequence.push_back(strandPaths[currentSequence][p2]);
                p2++;
            } else {
                strandSequence.push_back(-1);
            }
            p1++;
        }
        // Adding gaps towards the end
        while(strandSequence.size() < topoArray.size()) {
            strandSequence.push_back(-1);
        }
        alignedStrandSequences.push_back(strandSequence);
        currentSequence++;
    }
    return alignedStrandSequences;
}

std::vector< size_t > panmanUtils::GfaGraph::getTopologicalSort() {
    return topoSortedIntSequences;
}


panmanUtils::Pangraph::Pangraph(Json::Value& pangraphData) {
    // load paths
    bool circular=false;
    for(size_t i = 0; i < pangraphData["paths"].size(); i++) {
        Json::Value path = pangraphData["paths"][(int)i];
        for(size_t j = 0; j < path["blocks"].size(); j++) {
            paths[path["name"].asString()].push_back(path["blocks"][(int)j]["id"].asString());
            strandPaths[path["name"].asString()].push_back(path["blocks"][(int)j]["strand"].asBool());
        }
        if(path["circular"].asBool() == true) {
            circular = true;
            circularSequences[path["name"].asString()] = -path["offset"].asInt();
        }
    }

    std::unordered_map<std::string, int> blockSizeMap;

    // load blocks
    for(size_t i = 0; i < pangraphData["blocks"].size(); i++) {
        std::string blockId = pangraphData["blocks"][(int)i]["id"].asString();
        std::string sequence = pangraphData["blocks"][(int)i]["sequence"].asString();
        std::transform(sequence.begin(), sequence.end(),sequence.begin(), ::toupper);
        stringIdToConsensusSeq[blockId] = sequence;
        blockSizeMap[blockId] = sequence.size();
        std::vector< std::string > gapMemberNames = pangraphData["blocks"][(int)i]["gaps"].getMemberNames();
        for(auto member: gapMemberNames) {
            stringIdToGaps[blockId].push_back( std::make_pair( std::stoi(member), pangraphData["blocks"][(int)i]["gaps"][member].asInt() ) );
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["mutate"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["mutate"][(int)j][1].size(); k++) {
                std::string mutationString = pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                substitutions[blockId][seqName][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt(), mutationString) );
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
    std::unordered_map<std::string, std::vector<string>> test;
    if (circular) {

        std::vector<std::string> sample_base = {};
        int seq_count = 0;
        std::string sample_base_string;

        std::vector<std::string> sample_new = {};
        for(const auto& p: paths) {
            // std::cout << p.first << "\n";

            test[p.first] = p.second;
            if (seq_count == 0) {
                std::unordered_map< std::string, size_t > baseBlockNumber;
                sequenceInverted[p.first] = false;
                rotationIndexes[p.first] = 0;

                for(const auto& block: p.second) {
                    blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                    baseBlockNumber[block]++;
                    sample_base.push_back(block);
                }
            } else {
                // Assigning block numbers
                std::unordered_map< std::string, size_t > baseBlockNumber;
                for(const auto& block: p.second) {
                    blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                    baseBlockNumber[block]++;
                }

                std::vector<std::string> sample_dumy = {};
                sample_new.clear();
                for(const auto& block: p.second) {
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

            for(const auto& block: p.second) {
                blockNumbers[p.first].push_back(baseBlockNumber[block]+1);
                baseBlockNumber[block]++;
            }
        }
    }

    // Auto increment ID to assign to nodes
    numNodes = 0;

    std::unordered_map<int,std::string> intToString; // Locally stored -> Later on mapped to intIdToStringId
    int seqCount = 0; // Current Sequence ID
    std::vector<std::string> consensus = {};
    std::vector<std::string> sample = {};
    std::vector<std::string> consensus_new = {};
    std::vector<int> intSequenceConsensus= {};
    std::vector<int> intSequenceSample= {};
    std::vector<int> intSequenceConsensusNew= {};

    std::cout << "Resolving rearrangements and duplications..." << std::endl;
    for(const auto& p: paths) {
        if (seqCount == 0) { // Load first sequence path
            for(const auto& block: p.second) {
                consensus.push_back(block);
                // sample_base.push_back(block);
                intToString[numNodes] = block;
                intSequences[p.first].push_back(numNodes);
                intSequenceConsensus.push_back(numNodes);
                numNodes++;
            }
        } else {
            intSequenceSample.clear();
            intSequenceConsensusNew.clear();
            sample.clear();
            consensus_new.clear();

            for(const auto& block: p.second) {
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
                         intSequenceConsensusNew,
                         intToString);
            for (auto &b: intSequenceSample) {
                intSequences[p.first].push_back(b);
            }
            consensus.clear();
            intSequenceConsensus.clear();
            for (auto &b: consensus_new) {
                consensus.push_back(b);
            }

            for (auto &b: intSequenceConsensusNew) {
                intSequenceConsensus.push_back(b);
            }

        }
        seqCount++;
        // std::cout << seqCount << " " << intSequenceConsensusNew.size() << endl;
    }

    // re-assigning IDs in fixed order
    int reorder = 0;
    std::unordered_map<int,int> order_map = {};
    for (auto &i: intSequenceConsensus) {
        order_map[i] = reorder;
        intIdToStringId[reorder] = intToString[i];
        topoSortedIntSequences.push_back(reorder);
        reorder++;
    }

    for (auto &m: intSequences) {
        for (auto &s: m.second) {
            s = order_map[s];
        }
    }
}

std::unordered_map< std::string,std::vector< int > > panmanUtils::Pangraph::getAlignedStrandSequences(const std::vector< size_t >& topoArray) {
    std::unordered_map< std::string, std::vector< int > > alignedStrandSequences;
    for(auto p: intSequences) {
        size_t p1 = 0, p2 = 0;
        while(p1 < topoArray.size() && p2 < p.second.size()) {
            if(topoArray[p1] == p.second[p2]) {
                alignedStrandSequences[p.first].push_back(strandPaths[p.first][p2]);
                p2++;
            } else {
                alignedStrandSequences[p.first].push_back(-1);
            }
            p1++;
        }
        while(alignedStrandSequences[p.first].size() < topoArray.size()) {
            alignedStrandSequences[p.first].push_back(-1);
        }
    }
    return alignedStrandSequences;
}

std::unordered_map< std::string, std::vector< int > > panmanUtils::Pangraph::getAlignedSequences(const std::vector< size_t >& topoArray) {
    std::unordered_map< std::string, std::vector< int > > alignedSequences;
    for(auto p: intSequences) {
        size_t p1 = 0, p2 = 0;
        while(p1 < topoArray.size() && p2 < p.second.size()) {
            if(topoArray[p1] == p.second[p2]) {
                alignedSequences[p.first].push_back(topoArray[p1]);
                p2++;
            } else {
                alignedSequences[p.first].push_back(-1);
            }
            p1++;
        }
        while(alignedSequences[p.first].size() < topoArray.size()) {
            alignedSequences[p.first].push_back(-1);
        }
    }
    return alignedSequences;
}

std::vector< size_t > panmanUtils::Pangraph::getTopologicalSort() {
    std::vector< size_t > topoArray;

    for (auto &i: topoSortedIntSequences) {
        topoArray.push_back(i);
    }

    return topoArray;
}

panmanUtils::TreeGroup::TreeGroup(std::vector< Tree* >& tg) {
    for (auto& t: tg) {
        trees.push_back(*t);
    }
}

panmanUtils::TreeGroup::TreeGroup(std::vector< Tree* >& tg, std::ifstream& mutationFile) {

    for (auto& t: tg) {
        trees.push_back(*t);
    }

    // mutation file format: mutation type (H or R), tree_1 index, sequence_1 name, tree_2 index, sequence_2 name, start_point_1, end_point_1, start_point_2, end_point_2, tree_3 index (child tree), sequence_3 (child sequence) name
    std::string line;
    while(getline(mutationFile, line, '\n')) {
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
        std::string sequenceId3 = tokens[10];
        bool splitOccurred = false;

        if(treeIndex3 == treeIndex1 && treeIndex3 == treeIndex2) {
            // If all three sequences are from the same tree, split this tree
            std::cout << "Performing Split" << std::endl;
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex1] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        } else if(treeIndex3 == treeIndex1) {
            // If child belongs to one parent's tree, split this tree
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex1] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        } else if(treeIndex3 == treeIndex2) {
            // If child belongs to one parent's tree, split this tree
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex2] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        }

        sequence_t sequence1, sequence2;
        blockExists_t blockExists1, blockExists2;
        blockStrand_t blockStrand1, blockStrand2;
        trees[treeIndex1].getSequenceFromReference(sequence1, blockExists1, blockStrand1, sequenceId1, true);
        trees[treeIndex2].getSequenceFromReference(sequence2, blockExists2, blockStrand2, sequenceId2, true);

        int64_t co1 = 0, co2 = 0;
        if(trees[treeIndex1].circularSequences.find(sequenceId1) != trees[treeIndex1].circularSequences.end()) {
            co1 = trees[treeIndex1].circularSequences[sequenceId1];
        }
        if(trees[treeIndex2].circularSequences.find(sequenceId2) != trees[treeIndex2].circularSequences.end()) {
            co2 = trees[treeIndex2].circularSequences[sequenceId2];
        }

        if(!splitOccurred) {
            trees[treeIndex3].reroot(sequenceId3);
        }

        std::tuple< int,int,int,int > t_start1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(startPoint1, sequence1, blockExists1, blockStrand1, co1);
        std::tuple< int,int,int,int > t_end1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(endPoint1, sequence1, blockExists1, blockStrand1, co1);
        std::tuple< int,int,int,int > t_start2 = trees[treeIndex2].globalCoordinateToBlockCoordinate(startPoint2, sequence2, blockExists2, blockStrand2, co2);
        std::tuple< int,int,int,int > t_end2 = trees[treeIndex2].globalCoordinateToBlockCoordinate(endPoint2, sequence2, blockExists2, blockStrand2, co2);

        complexMutations.emplace_back(mutationType, treeIndex1, treeIndex2, treeIndex3, sequenceId1, sequenceId2, sequenceId3, t_start1, t_end1, t_start2, t_end2);
    }
}


panmanUtils::TreeGroup::TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile) {
    for(size_t i = 0; i < treeFiles.size(); i++) {
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
        inPMATBuffer.push(boost::iostreams::gzip_decompressor());
        inPMATBuffer.push(treeFiles[i]);
        std::istream inputStream(&inPMATBuffer);

        trees.emplace_back(inputStream);
    }

    // mutation file format: mutation type (H or R), tree_1 index, sequence_1 name, tree_2 index, sequence_2 name, start_point_1, end_point_1, start_point_2, end_point_2, tree_3 index (child tree), sequence_3 (child sequence) name
    std::string line;
    while(getline(mutationFile, line, '\n')) {
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
        std::string sequenceId3 = tokens[10];
        bool splitOccurred = false;

        if(treeIndex3 == treeIndex1 && treeIndex3 == treeIndex2) {
            // If all three sequences are from the same tree, split this tree
            std::cout << "Performing Split" << std::endl;
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex1] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        } else if(treeIndex3 == treeIndex1) {
            // If child belongs to one parent's tree, split this tree
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex1] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        } else if(treeIndex3 == treeIndex2) {
            // If child belongs to one parent's tree, split this tree
            std::pair< panmanUtils::Tree, panmanUtils::Tree > parentAndChild = trees[treeIndex1].splitByComplexMutations(sequenceId3);
            splitOccurred = true;
            trees[treeIndex2] = parentAndChild.first;
            trees.push_back(parentAndChild.second);
            treeIndex3 = trees.size()-1;
        }

        sequence_t sequence1, sequence2;
        blockExists_t blockExists1, blockExists2;
        blockStrand_t blockStrand1, blockStrand2;
        trees[treeIndex1].getSequenceFromReference(sequence1, blockExists1, blockStrand1, sequenceId1, true);
        trees[treeIndex2].getSequenceFromReference(sequence2, blockExists2, blockStrand2, sequenceId2, true);

        int64_t co1 = 0, co2 = 0;
        if(trees[treeIndex1].circularSequences.find(sequenceId1) != trees[treeIndex1].circularSequences.end()) {
            co1 = trees[treeIndex1].circularSequences[sequenceId1];
        }
        if(trees[treeIndex2].circularSequences.find(sequenceId2) != trees[treeIndex2].circularSequences.end()) {
            co2 = trees[treeIndex2].circularSequences[sequenceId2];
        }

        if(!splitOccurred) {
            trees[treeIndex3].reroot(sequenceId3);
        }

        std::tuple< int,int,int,int > t_start1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(startPoint1, sequence1, blockExists1, blockStrand1, co1);
        std::tuple< int,int,int,int > t_end1 = trees[treeIndex1].globalCoordinateToBlockCoordinate(endPoint1, sequence1, blockExists1, blockStrand1, co1);
        std::tuple< int,int,int,int > t_start2 = trees[treeIndex2].globalCoordinateToBlockCoordinate(startPoint2, sequence2, blockExists2, blockStrand2, co2);
        std::tuple< int,int,int,int > t_end2 = trees[treeIndex2].globalCoordinateToBlockCoordinate(endPoint2, sequence2, blockExists2, blockStrand2, co2);

        complexMutations.emplace_back(mutationType, treeIndex1, treeIndex2, treeIndex3, sequenceId1, sequenceId2, sequenceId3, t_start1, t_end1, t_start2, t_end2);
    }
}

panmanUtils::TreeGroup::TreeGroup(std::istream& fin) {
    panman::treeGroup TG;

    if(!TG.ParseFromIstream(&fin)) {
        throw std::invalid_argument("Could not read tree group from input file.");
    }

    for(int i = 0; i < TG.trees_size(); i++) {
        trees.emplace_back(TG.trees(i));
    }
    for(int i = 0; i < TG.complexmutations_size(); i++) {
        complexMutations.emplace_back(TG.complexmutations(i));
    }
}

void panmanUtils::TreeGroup::printFASTA(std::ofstream& fout) {
    for(auto& tree: trees) {
        tree.printFASTA(fout);
    }
}

void panmanUtils::TreeGroup::writeToFile(std::ostream& fout) {
    panman::treeGroup treeGroupToWrite;

    for(auto& tree: trees) {
        panman::tree treeToWrite;
        Node* node = tree.root;

        tree.getNodesPreorder(node, treeToWrite);
        std::string newick = tree.getNewickString(node);

        treeToWrite.set_newick(newick);
        std::map< std::vector< uint32_t >, std::vector< std::pair< int64_t, bool > > >
        consensusSeqToBlockIds;

        for(auto block: tree.blocks) {
            int64_t blockId;
            bool blockGapExists = false;
            if(block.secondaryBlockId != -1) {
                blockId = ((int64_t)block.primaryBlockId << 32) + block.secondaryBlockId;
                blockGapExists = true;
            } else {
                blockId = ((int64_t)block.primaryBlockId << 32);
            }
            consensusSeqToBlockIds[block.consensusSeq].push_back(
                std::make_pair(blockId, blockGapExists));
        }

        for(auto u: consensusSeqToBlockIds) {
            panman::consensusSeqToBlockIds c;
            for(auto v: u.first) {
                c.add_consensusseq(v);
            }
            for(auto v: u.second) {
                c.add_blockid(v.first);
                c.add_blockgapexist(v.second);
            }
            treeToWrite.add_consensusseqmap();
            *treeToWrite.mutable_consensusseqmap( treeToWrite.consensusseqmap_size() - 1 ) = c;
        }
        for(size_t i = 0; i < tree.gaps.size(); i++) {
            panman::gapList gl;
            for(size_t j = 0; j < tree.gaps[i].nucPosition.size(); j++) {
                gl.add_nucposition(tree.gaps[i].nucPosition[j]);
                gl.add_nucgaplength(tree.gaps[i].nucGapLength[j]);
            }
            if(tree.gaps[i].secondaryBlockId != -1) {
                gl.set_blockid(((int64_t)tree.gaps[i].primaryBlockId << 32) + tree.gaps[i]
                               .secondaryBlockId);
                gl.set_blockgapexist(true);
            } else {
                gl.set_blockid(((int64_t)tree.gaps[i].primaryBlockId << 32));
                gl.set_blockgapexist(false);
            }
            treeToWrite.add_gaps();
            *treeToWrite.mutable_gaps( treeToWrite.gaps_size() - 1 ) = gl;
        }
        for(auto u: tree.circularSequences) {
            panman::circularOffset co;
            co.set_sequenceid(u.first);
            co.set_offset(u.second);
            treeToWrite.add_circularsequences();
            *treeToWrite.mutable_circularsequences(treeToWrite.circularsequences_size()-1) = co;
        }

        for(auto u: tree.rotationIndexes) {
            panman::rotationIndex ri;
            ri.set_sequenceid(u.first);
            ri.set_blockoffset(u.second);
            treeToWrite.add_rotationindexes();
            *treeToWrite.mutable_rotationindexes(treeToWrite.rotationindexes_size()-1) = ri;
        }

        for(auto u: tree.sequenceInverted) {
            panman::sequenceInverted si;
            si.set_sequenceid(u.first);
            si.set_inverted(u.second);
            treeToWrite.add_sequencesinverted();
            *treeToWrite.mutable_sequencesinverted(treeToWrite.sequencesinverted_size()-1) = si;
        }

        treeGroupToWrite.add_trees();
        *treeGroupToWrite.mutable_trees( treeGroupToWrite.trees_size() - 1 ) = treeToWrite;
    }

    for(auto cm: complexMutations) {
        treeGroupToWrite.add_complexmutations();
        *treeGroupToWrite.mutable_complexmutations(treeGroupToWrite
                .complexmutations_size()-1) = cm.toProtobuf();
    }

    if(!treeGroupToWrite.SerializeToOstream(&fout)) {
        std::cerr << "Failed to write to output file." << std::endl;
    }
}

void panmanUtils::TreeGroup::printComplexMutations(std::ostream& fout) {
    for(const auto& u: complexMutations) {
        sequence_t s1, s2;
        blockExists_t b1, b2;
        blockStrand_t str1, str2;

        // Circular Offsets
        int co1 = 0, co2 = 0;

        trees[u.treeIndex1].getSequenceFromReference(s1, b1, str1, u.sequenceId1, true);
        trees[u.treeIndex2].getSequenceFromReference(s2, b2, str2, u.sequenceId2, true);

        if(trees[u.treeIndex1].circularSequences.find(u.sequenceId1) != trees[u.treeIndex1].circularSequences.end()) {
            co1 = trees[u.treeIndex1].circularSequences[u.sequenceId1];
        }

        if(trees[u.treeIndex2].circularSequences.find(u.sequenceId2) != trees[u.treeIndex2].circularSequences.end()) {
            co2 = trees[u.treeIndex2].circularSequences[u.sequenceId2];
        }

        fout << u.mutationType
             << " " << u.treeIndex1
             << " " << u.sequenceId1
             << " " << u.treeIndex2
             << " " << u.sequenceId2
             << " " << trees[u.treeIndex1].getUnalignedGlobalCoordinate(u.primaryBlockIdStart1,
                     u.secondaryBlockIdStart1, u.nucPositionStart1, u.nucGapPositionStart1, s1, b1,
                     str1, co1)
             << " " << trees[u.treeIndex1].getUnalignedGlobalCoordinate(u.primaryBlockIdEnd1,
                     u.secondaryBlockIdEnd1, u.nucPositionEnd1, u.nucGapPositionEnd1, s1, b1,
                     str1, co1)
             << " " << trees[u.treeIndex2].getUnalignedGlobalCoordinate(u.primaryBlockIdStart2,
                     u.secondaryBlockIdStart2, u.nucPositionStart2, u.nucGapPositionStart2, s2, b2,
                     str2, co2)
             << " " << trees[u.treeIndex2].getUnalignedGlobalCoordinate(u.primaryBlockIdEnd2,
                     u.secondaryBlockIdEnd2, u.nucPositionEnd2, u.nucGapPositionEnd2, s2, b2,
                     str2, co2)
             << " " << u.treeIndex3 << " " << u.sequenceId3 << "\n";
    }
}
