
#include "panmanUtils.hpp"

void panmanUtils::Tree::generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout) {
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

void panmanUtils::Tree::printMAF(std::ostream& fout) {
    std::vector< std::string > sequenceNames;

    std::map< std::vector< uint32_t >, std::vector< std::pair< int,int > > > blocksWithSameSequences;

    for(auto b: blocks) {
        blocksWithSameSequences[b.consensusSeq].push_back(std::make_pair(b.primaryBlockId, b.secondaryBlockId));
    }

    // sequence name, primary bid, secondary bid -> actual sequence offset
    std::map< std::tuple< std::string, int32_t, int32_t >, int32_t > sequenceBlockToStartPoint;
    std::map< std::string, int32_t > sequenceLengths;

    int ctr2=0;

    for(auto u: allNodes) {
        if(u.second->children.size() == 0) {
            sequenceNames.push_back(u.first);
            sequence_t sequence;
            blockExists_t blockExists;
            blockStrand_t blockStrand;
            int rotIndex = 0;
            getSequenceFromReference(sequence, blockExists, blockStrand, u.first, true, &rotIndex);

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
                if(blockExists[i].first) {
                    sequenceBlockToStartPoint[std::make_tuple(u.first, blockIds[i], -1)] = ctr;
                    for(size_t j = 0; j < sequence[i].first.size(); j++) {
                        for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                            if(sequence[i].first[j].second[k] != '-' && sequence[i].first[j].second[k] != 'x') {
                                ctr++;
                            }
                            ctr2++;
                        }
                        if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                            ctr++;
                        }
                        ctr2++;
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

    std::atomic<int> ctr3 = 0;

    fout << "##maf version=1\n";

    for(auto& common: blocksWithSameSequences) {
        fout << "a\n";
        for(auto& b: common.second) {
            tbb::concurrent_unordered_map< std::string, std::pair< std::pair< int, int >, std::pair< std::string, bool > > > sequenceIdToSequence;

            int primaryBlockId = b.first;
            int secondaryBlockId = b.second;

            tbb::parallel_for_each(sequenceNames, [&](auto u) {
                if(sequenceBlockToStartPoint.find(std::make_tuple(u, primaryBlockId, -1)) == sequenceBlockToStartPoint.end()) {
                    return;
                }

                block_t sequence;
                bool blockExists = false, blockStrand = true;

                getBlockSequenceFromReference(sequence, blockExists, blockStrand, u, primaryBlockId,
                                              secondaryBlockId);

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
                ctr3+=stringSequence.length();
                sequenceIdToSequence[u] = std::make_pair(std::make_pair(primaryBlockId, secondaryBlockId), std::make_pair(stringSequence, blockStrand));
            });

            for(auto u: sequenceIdToSequence) {
                fout << "s\t" << u.first << "\t"<< sequenceBlockToStartPoint[std::make_tuple(u.first, u.second.first.first, u.second.first.second)] <<"\t" << stripGaps(u.second.second.first).length() << "\t" << (u.second.second.second? "+\t":"-\t") << sequenceLengths[u.first] << "\t" << u.second.second.first << "\n";
            }
        }
        fout << "\n";
    }
}
