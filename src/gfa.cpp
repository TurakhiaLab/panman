#include "panmanUtils.hpp"

void panmanUtils::Tree::convertToGFA(std::ostream& fout) {
    // First we check if there are any nucleotide mutations. If there are no nuc mutations, we can
    // simply construct a GFA of blocks
    bool nucMutationFlag = false;
    for(auto u: allNodes) {
        if(u.second->nucMutation.size() != 0) {
            nucMutationFlag = true;
        }
    }

    if(!nucMutationFlag) {
        // get nodes
        std::map<std::pair<int32_t, int32_t>, std::string> nodes;
        for(auto block: blocks) {
            int64_t primaryBlockId = block.primaryBlockId;
            int64_t secondaryBlockId = block.secondaryBlockId;
            std::string sequenceString;
            for(auto u: block.consensusSeq) {
                for(size_t k = 0; k < 8; k++) {
                    const int nucCode = (((u) >> (4*(7 - k))) & 15);
                    if(nucCode == 0) {
                        break;
                    }
                    const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
                    sequenceString += nucleotide;
                }
            }
            nodes[std::make_pair(primaryBlockId, secondaryBlockId)] = sequenceString;
        }

        // block presense map
        std::vector< std::pair< bool, std::vector< bool > > > blockExistsGlobal(blocks.size() + 1, {false, {}});
        // Assigning block gaps
        for(size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
            blockExistsGlobal[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        }

        tbb::concurrent_unordered_set< std::pair< std::pair<int32_t, int32_t>, std::pair<int32_t, int32_t> > > edges;
        tbb::concurrent_unordered_map< std::string, std::vector< std::pair<int32_t, int32_t> > > paths;

        // get all paths
        tbb::parallel_for_each(allNodes, [&](auto u) {
            if(u.second->children.size()) {
                return;
            }

            auto blockExists = blockExistsGlobal;

            std::vector< panmanUtils::Node* > path;

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

                    if(type == panmanUtils::BlockMutationType::BI) {
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
        size_t node_len = 32;
        size_t autoIncrId = 0;
        std::map< std::pair< std::tuple< int, size_t, size_t >, std::string >,
            std::pair< size_t, bool > > allSequenceNodes;
        std::mutex allSequenceNodeMutex;
        tbb::concurrent_unordered_map< std::string, std::vector< size_t > > paths;
        tbb::concurrent_unordered_map< std::string, std::vector< bool > > strandPaths;

        // for(const auto& u: allNodes) {
        tbb::parallel_for_each(allNodes, [&](const auto& u) {
            if(u.second->children.size() != 0) {
                return;
                // continue;
            }

            sequence_t sequence;
            blockExists_t blockExists;
            blockStrand_t blockStrand;
            getSequenceFromReference(sequence, blockExists, blockStrand, u.first);

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
                                if(currentSequence.length() == node_len) {
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
                            if(currentSequence.length() == node_len) {
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
                            currentStart = std::make_tuple(-1*(int)i-1,j,-1);
                            if(currentSequence.length() == node_len) {
                                currentSequence = stripGaps(currentSequence);
                                if(currentSequence.length()) {
                                    // Since the GFA stores the strand parameter, the reverse
                                    // complement will be computed anyway
                                    std::reverse(currentSequence.begin(), currentSequence.end());
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
                                currentStart = std::make_tuple(-1*(int)i,j,k);
                                if(currentSequence.length() == node_len) {
                                    currentSequence = stripGaps(currentSequence);
                                    if(currentSequence.length()) {
                                        // Since the GFA stores the strand parameter, the reverse
                                        // complement will be computed anyway
                                        std::reverse(currentSequence.begin(), currentSequence.end());
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
        });
        // }

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
                if(u.first.second) {
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

        // Remove duplicate nodes
        std::map< std::string, size_t > sequenceToId;
        size_t ctr = 1;
        for(const auto& u: finalNodes) {
            sequenceToId[u.second] = ctr++;
        }

        std::unordered_map< size_t, size_t > oldToNew;
        for(const auto& u: finalNodes) {
            if(sequenceToId.find(u.second) != sequenceToId.end()) {
                oldToNew[u.first.first] = sequenceToId[u.second];
            }
            // oldToNew[u.first.first] = ctr++;
        }

        std::set< std::pair< pair< size_t, bool >, pair< size_t, bool > > > edges;

        for(const auto& u: G) {
            if(finalNodes.find(u.first) == finalNodes.end()) {
                continue;
            }
            for(auto edge: u.second) {
                if(finalNodes.find(edge.second) == finalNodes.end()) {
                    continue;
                }
                edges.insert(std::make_pair(std::make_pair(oldToNew[u.first.first],
                                            u.first.second), std::make_pair(oldToNew[edge.second.first],
                                                    edge.second.second)));
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
            p.second = newPath;
            strandPaths[p.first] = newStrandPath;
        }

        // for(auto& p: paths) {
        //     fout << ">" << p.first << "\n";
        //     for(int i = 0; i < p.second.size(); i++) {
        //         fout << finalNodes[std::make_pair(p.second[i], strandPaths[p.first][i])];
        //     }
        //     fout << "\n";
        // }

        // convert node IDs to consecutive node IDs
        size_t currentID = 1;
        std::unordered_map< size_t, size_t > sequentialIds;
        for(auto u: finalNodes) {
            if(sequentialIds.find(oldToNew[u.first.first]) == sequentialIds.end()) {
                sequentialIds[oldToNew[u.first.first]] = currentID;
                currentID++;
            }
        }

        fout << "H\tVN:Z:1.1\n";
        std::unordered_map< size_t, bool > alreadyPrinted;
        for(auto u: finalNodes) {
            if(!alreadyPrinted[oldToNew[u.first.first]]) {
                alreadyPrinted[oldToNew[u.first.first]] = true;
                fout << "S\t" << sequentialIds[oldToNew[u.first.first]] << "\t" << u.second << "\n";
            }
        }

        for(auto u: edges) {
            fout << "L\t" << sequentialIds[u.first.first] << "\t" << (u.first.second? "+":"-")
                 << "\t" << sequentialIds[u.second.first] << "\t" << (u.second.second? "+":"-")
                 << "\t0M\n";
        }

        for(auto u: paths) {
            fout << "P\t" << u.first << "\t";
            for(size_t i = 0; i < u.second.size(); i++) {
                fout << sequentialIds[oldToNew[u.second[i]]] << (strandPaths[u.first][i]?"+":"-");
                if(i != u.second.size() - 1) {
                    fout << ",";
                }
            }
            fout << "\t*\n";
        }
    }
}

