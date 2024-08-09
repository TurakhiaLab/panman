#include "panmanUtils.hpp"


void panmanUtils::Tree::reroot(std::string sequenceName) {
    if(allNodes.find(sequenceName) == allNodes.end()) {
        std::cerr << "Sequence with name " << sequenceName << " not found!" << std::endl;
        return;
    }
    panmanUtils::Node* newRoot = allNodes[sequenceName];
    if(newRoot->children.size()) {
        std::cerr << "Node with id " << sequenceName << " is not a tip!" << std::endl;
        return;
    }
    sequence_t sequence;
    blockExists_t blockExists;
    blockStrand_t blockStrand;

    getSequenceFromReference(sequence, blockExists, blockStrand, sequenceName);

    std::unordered_map<std::string, sequence_t> nodeIdToSequence;
    std::unordered_map<std::string, blockExists_t> nodeIdToBlockExists;
    std::unordered_map<std::string, blockStrand_t> nodeIdToBlockStrand;

    for(const auto& u: allNodes) {
        if(u.second->children.size() == 0) {
            sequence_t tempSequence;
            blockExists_t tempBlockExists;
            blockExists_t tempBlockStrand;

            getSequenceFromReference(tempSequence, tempBlockExists, tempBlockStrand, u.first);
            nodeIdToSequence[u.first] = tempSequence;
            nodeIdToBlockExists[u.first] = tempBlockExists;
            nodeIdToBlockStrand[u.first] = tempBlockStrand;
        }
    }

    // Transform tree topology
    transform(newRoot);

    std::cout << "Transformation complete!" << std::endl;

    // clear previous block mutations
    for(auto u: allNodes) {
        u.second->blockMutation.clear();
    }

    std::unordered_map< std::string, std::mutex > nodeMutexes;

    for(auto u: allNodes) {
        nodeMutexes[u.first];
    }

    // make new block mutations
    tbb::parallel_for((size_t)0, blockExists.size(), [&](size_t i) {
        tbb::parallel_for((size_t)0, blockExists[i].second.size(), [&](size_t j) {
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;
            // block gaps
            for(const auto& u: nodeIdToBlockExists) {
                if(!u.second[i].second[j]) {
                    // doesn't exist
                    states[u.first] = 1;
                } else if(nodeIdToBlockStrand[u.first][i].second[j]) {
                    // forward strand
                    states[u.first] = 2;
                } else {
                    // reverse strand
                    states[u.first] = 4;
                }
            }
            int defaultState;
            if(!blockExists[i].second[j]) {
                defaultState = 1;
            } else if(blockStrand[i].second[j]) {
                defaultState = 2;
            } else {
                defaultState = 4;
            }

            blockFitchForwardPassNew(root, states);
            blockFitchBackwardPassNew(root, states, 1, defaultState);
            blockFitchAssignMutationsNew(root, states, mutations, 1);

            for(auto u: mutations) {
                nodeMutexes[u.first].lock();
                allNodes[u.first]->blockMutation.emplace_back(i, u.second, j);
                nodeMutexes[u.first].unlock();
            }
        });
        std::unordered_map< std::string, int > states;
        std::unordered_map< std::string, std::pair< BlockMutationType, bool > > mutations;
        // main block
        for(const auto& u: nodeIdToBlockExists) {
            if(!u.second[i].first) {
                // doesn't exist
                states[u.first] = 1;
            } else if(nodeIdToBlockStrand[u.first][i].first) {
                // forward strand
                states[u.first] = 2;
            } else {
                // reverse strand
                states[u.first] = 4;
            }
        }
        int defaultState;
        if(!blockExists[i].first) {
            defaultState = 1;
        } else if(blockStrand[i].first) {
            defaultState = 2;
        } else {
            defaultState = 4;
        }

        blockFitchForwardPassNew(root, states);
        blockFitchBackwardPassNew(root, states, 1, defaultState);
        blockFitchAssignMutationsNew(root, states, mutations, 1);
        for(auto u: mutations) {
            nodeMutexes[u.first].lock();
            allNodes[u.first]->blockMutation.emplace_back(i, u.second, -1);
            nodeMutexes[u.first].unlock();
        }
    });

    std::cout << "Block Mutations Ready!" << std::endl;

    // clear previous nuc mutations
    for(auto u: allNodes) {
        u.second->nucMutation.clear();
    }

    tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > nonGapMutations;
    tbb::concurrent_unordered_map< std::string, std::vector< std::tuple< int,int,int,int,int,int > > > gapMutations;

    tbb::parallel_for((size_t)0, sequence.size(), [&](size_t i) {
        // main block
        std::vector< char > consensusSeq;
        bool endFlag = false;
        for(size_t t1 = 0; t1 < blocks.size(); t1++) {
            if(blocks[t1].primaryBlockId == (int)i && blocks[t1].secondaryBlockId == -1) {
                for(size_t t2 = 0; t2 < blocks[t1].consensusSeq.size(); t2++) {
                    for(size_t t3 = 0; t3 < 8; t3++) {
                        const int nucCode = (((blocks[t1].consensusSeq[t2]) >> (4*(7 - t3))) & 15);
                        if(nucCode == 0) {
                            endFlag = true;
                            break;
                        }
                        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

                        consensusSeq.push_back(nucleotide);
                    }
                    if(endFlag) {
                        break;
                    }
                }
                endFlag = true;
                break;
            }
        }

        if(!endFlag) {
            std::cerr << "FATAL: Block with id " << i << " " << -1 << " not found!" << std::endl;
            exit(-1);
        }

        consensusSeq.push_back('-');
        if(consensusSeq.size() != sequence[i].first.size()) {
            std::cerr << "FATAL: consenusSeq length doesn't match with sequence length" << std::endl;
            exit(-1);
        }

        tbb::parallel_for((size_t)0, sequence[i].first.size(), [&](size_t k) {
            tbb::parallel_for((size_t)0, sequence[i].first[k].second.size(), [&](size_t w) {
                // gap nuc
                std::unordered_map< std::string, int > states;
                std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;
                for(const auto& u: nodeIdToSequence) {
                    if(u.second[i].first[k].second[w] != '-' && u.second[i].first[k].second[w] != 'x') {
                        states[u.first] = (1 << getCodeFromNucleotide(u.second[i].first[k].second[w]));
                    }  else {
                        states[u.first] = 1;
                    }
                }
                int nucleotideCode = 1;
                if(sequence[i].first[k].second[w] != '-' && sequence[i].first[k].second[w] != 'x') {
                    nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].first[k].second[w]));
                }

                nucFitchForwardPass(root, states);
                nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
                nucFitchAssignMutations(root, states, mutations, 1);

                for(auto mutation: mutations) {
                    nodeMutexes[mutation.first].lock();
                    gapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, k, w, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
                    nodeMutexes[mutation.first].unlock();
                }
            });
            // main nuc
            std::unordered_map< std::string, int > states;
            std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > > mutations;

            for(const auto& u: nodeIdToSequence) {
                if(u.second[i].first[k].first != '-' && u.second[i].first[k].first != 'x') {
                    states[u.first] = (1 << getCodeFromNucleotide(u.second[i].first[k].first));
                }  else {
                    states[u.first] = 1;
                }
            }
            int nucleotideCode = 1;
            if(sequence[i].first[k].first != '-' && sequence[i].first[k].first != 'x') {
                nucleotideCode = (1 << getCodeFromNucleotide(sequence[i].first[k].first));
            }

            nucFitchForwardPass(root, states);
            nucFitchBackwardPass(root, states, nucleotideCode, nucleotideCode);
            nucFitchAssignMutations(root, states, mutations, (1 << getCodeFromNucleotide(consensusSeq[k])));

            for(auto mutation: mutations) {
                nodeMutexes[mutation.first].lock();
                nonGapMutations[mutation.first].push_back(std::make_tuple((int)i, -1, k, -1, mutation.second.first, getCodeFromNucleotide(mutation.second.second)));
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
            if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<1>(u.second[i]) != std::get<1>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])) {
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
            if(i - currentStart == 6 || std::get<0>(u.second[i]) != std::get<0>(u.second[i-1]) || std::get<1>(u.second[i]) != std::get<1>(u.second[i-1]) || std::get<2>(u.second[i]) != std::get<2>(u.second[i-1]) || std::get<3>(u.second[i]) != std::get<3>(u.second[i-1])+1 || std::get<4>(u.second[i]) != std::get<4>(u.second[i-1])) {
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
