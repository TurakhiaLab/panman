#include "panmanUtils.hpp"

void panmanUtils::printSequenceLines(const sequence_t& sequence,\
                                     const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize, bool aligned, std::ostream& fout, int offset, bool debug) {

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

std::string panmanUtils::printSequenceLinesNew(const std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::unordered_map<int, int>& blockLengths,
                          const std::vector<bool>& blockExists, 
                          const std::vector<bool>& blockStrand, size_t lineSize, bool aligned, int offset, bool debug) {

    // String that stores the sequence to be printed
    std::string line;

    for(size_t i = 0; i < blockExists.size(); i++) {
        // Non-gap block - the only type being used currently
        if(blockExists[i]) {
            // line += ">" + std::to_string(i) + "\n";
            // If forward strand
            if(blockStrand[i]) {
                // Iterate through main nucs
                for(size_t j = 0; j < sequence[i].size(); j++) {
                    // Gap nucs
                    for(size_t k = 0; k < sequence[i][j].second.size(); k++) {
                        if(sequence[i][j].second[k] != '-') {
                            line += sequence[i][j].second[k];
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                    // Main nuc
                    if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        line += sequence[i][j].first;
                    } else if(aligned && sequence[i][j].first != 'x') {
                        line += '-';
                    }
                }
            } else {
                // If reverse strand, iterate backwards
                for(size_t j = sequence[i].size()-1; j+1 > 0; j--) {
                    // Main nuc first since we are iterating in reverse direction
                    if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        line += getComplementCharacter(sequence[i][j].first);
                    } else if(aligned  && sequence[i][j].first != 'x') {
                        line += '-';
                    }

                    // Gap nucs
                    for(size_t k = sequence[i][j].second.size()-1; k+1 > 0; k--) {
                        if(sequence[i][j].second[k] != '-') {
                            line += getComplementCharacter(sequence[i][j].second[k]);
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                }
            }
        } else if(aligned) {
            // If aligned sequence is required, print gaps instead if block does not exist
            for(size_t j = 0; j < blockLengths[i]; j++) {
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

    // std::cout << line.size() << std::endl;

    // From offset to end
    std::string currentLine = "";
    int c=0;
    for(size_t i = ctr; i < line.length(); i++) {
        currentLine += line[i];
        c++;
        if(c == lineSize) {
            currentLine += '\n';
            c=0;
        }
    }
    // From beginning to offset
    for(size_t i = 0; i < ctr; i++) {
        currentLine += line[i];
        c++;
        if(c == lineSize) {
            currentLine += '\n';
            c=0;
        }
    }

    return currentLine;

}

void panmanUtils::printSubsequenceLines(const sequence_t& sequence,
                                     const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize, 
                                     const std::tuple<int, int, int, int>& panMATStart, 
                                     const std::tuple<int, int, int, int>& panMATEnd,
                                     bool aligned, std::ostream& fout, int offset, bool debug) {

    int primaryBlockIdStart = std::get<0>(panMATStart);
    int secondaryBlockIdStart = std::get<1>(panMATStart);
    int posStart = std::get<2>(panMATStart);
    int gapPosStart = std::get<3>(panMATStart);

    int primaryBlockIdEnd = std::get<0>(panMATEnd);
    int secondaryBlockIdEnd = std::get<1>(panMATEnd);
    int posEnd = std::get<2>(panMATEnd);
    int gapPosEnd = std::get<3>(panMATEnd);

    // String that stores the sequence to be printed
    std::string line;

    for(size_t i = primaryBlockIdStart; i < primaryBlockIdEnd-primaryBlockIdStart+1; i++) {
        
        // Non-gap block - the only type being used currently
        if(blockExists[i].first) {
            // If forward strand
            if(blockStrand[i].first) {
                // Iterate through main nucs
                size_t nucStart = (i==primaryBlockIdStart)? posStart: 0;
                size_t nucEnd = (i==primaryBlockIdEnd)? posEnd + 1: sequence[i].first.size();
                for(size_t j = nucStart; j < nucEnd; j++) {
                    // Gap nucs
                    size_t nucGapStart = (i==primaryBlockIdStart && j == posStart)? gapPosStart: 0;
                    size_t nucGapEnd = (i==primaryBlockIdEnd && j == posEnd)? gapPosEnd + 1: sequence[i].first[j].second.size();
                    for(size_t k = nucGapStart; k < sequence[i].first[j].second.size(); k++) {
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
                size_t nucStart = (i==primaryBlockIdStart)? posStart: 0;
                size_t nucEnd = (i==primaryBlockIdEnd)? posEnd: sequence[i].first.size() - 1;
                for(size_t j = nucEnd; j+1 > nucStart; j--) {
                    // Main nuc first since we are iterating in reverse direction
                    if(sequence[i].first[j].first != '-' && sequence[i].first[j].first != 'x') {
                        line += getComplementCharacter(sequence[i].first[j].first);
                    } else if(aligned) {
                        line += '-';
                    }

                    // Gap nucs
                    size_t nucGapStart = (i==primaryBlockIdStart && j == posStart)? gapPosStart: 0;
                    size_t nucGapEnd = (i==primaryBlockIdEnd && j == posEnd)? gapPosEnd: sequence[i].first[j].second.size() - 1;
                    for(size_t k = nucGapEnd; k+1 > nucGapStart; k--) {
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
            size_t nucStart = (i==primaryBlockIdStart)? posStart: 0;
            size_t nucEnd = (i==primaryBlockIdEnd)? posEnd + 1: sequence[i].first.size();
            for(size_t j = nucStart; j < nucEnd; j++) {
                for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                    line+='-';
                }
                line+='-';
            }
        }

    }

    std::cout << line << std::endl;
    // size_t ctr = 0;

    // if(offset != 0) {
    //     for(size_t i = 0; i < line.length(); i++) {
    //         if(line[i] != '-') {
    //             if(ctr == (size_t)offset) {
    //                 // mark starting point
    //                 ctr = i;
    //                 break;
    //             }
    //             ctr++;
    //         }
    //     }
    // }

    // // std::cout << line << std::endl;
    // std::string currentLine = "";
    // bool reachedEnd = false;
    // int newStart = (line.size()-1-ctr >= start)? ctr+start: start-line.size()-1-ctr;
    // int newEnd = (line.size()-1-ctr >= end)? ctr+end: end-line.size()-1-ctr;

    // // std::cout << newStart << " " << newEnd << " " << ctr << " " << start << " " << end << std::endl;
    // if (newStart > newEnd) {
    //     currentLine += line.substr(newStart, line.size()-newStart);
    //     currentLine += line.substr(0, newEnd+1);
    // } else {
    //     currentLine += line.substr(newStart, newEnd-newStart+1);
    // }
    // fout << currentLine << std::endl;
}

std::pair<std::vector<std::string>, std::vector<int>> panmanUtils::printSequenceLinesNewer(const std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::unordered_map<int, int>& blockLengths,
                          const std::vector<bool>& blockExists, 
                          const std::vector<bool>& blockStrand, size_t lineSize, bool aligned, int offset, bool debug) {

    // String that stores the sequence to be printed
    std::vector<std::string> lines;
    std::vector<int> blockLens;

    for(size_t i = 0; i < blockExists.size(); i++) {
        // Non-gap block - the only type being used currently
        if(blockExists[i]) {
            std::string line="";
            // If forward strand
            if(blockStrand[i]) {
                // Iterate through main nucs
                for(size_t j = 0; j < sequence[i].size(); j++) {
                    // Gap nucs
                    for(size_t k = 0; k < sequence[i][j].second.size(); k++) {
                        if(sequence[i][j].second[k] != '-') {
                            line += sequence[i][j].second[k];
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                    // Main nuc
                    if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        line += sequence[i][j].first;
                    } else if(aligned && sequence[i][j].first != 'x') {
                        line += '-';
                    }
                }
            } else {
                // If reverse strand, iterate backwards
                for(size_t j = sequence[i].size()-1; j+1 > 0; j--) {
                    // Main nuc first since we are iterating in reverse direction
                    if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        line += getComplementCharacter(sequence[i][j].first);
                    } else if(aligned  && sequence[i][j].first != 'x') {
                        line += '-';
                    }

                    // Gap nucs
                    for(size_t k = sequence[i][j].second.size()-1; k+1 > 0; k--) {
                        if(sequence[i][j].second[k] != '-') {
                            line += getComplementCharacter(sequence[i][j].second[k]);
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                }
            }

            blockLens.push_back(line.size());
            lines.push_back(line);

        } else if(aligned) {
            // If aligned sequence is required, print gaps instead if block does not exist
            blockLens.push_back(blockLengths[i]);
            lines.push_back("-");
        }
    }

    return std::make_pair(lines, blockLens);

}

// Depth first traversal FASTA writer
void panmanUtils::Tree::printFASTAHelper(panmanUtils::Node* root, sequence_t& sequence,
        blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {
        
    // Apply mutations
    // std::cout << root->identifier << " " << std::get<0>(panMATStart) << " " << std::get<1>(panMATStart) << " " << std::get<2>(panMATStart) << " " << std::get<3>(panMATStart) <<std::endl;
    // std::cout << root->identifier << " " << std::get<0>(panMATEnd) << " " << std::get<1>(panMATEnd) << " " << std::get<2>(panMATEnd) << " " << std::get<3>(panMATEnd) <<std::endl;

    // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
    std::vector< std::tuple< int32_t, int32_t, bool, bool, bool, bool > > blockMutationInfo;

    // Block Mutations
    for(auto mutation: root->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        if (secondaryBlockId != -1) {
            std::cout << "Error: Block Secondary ID is not -1" << std::endl;
            exit(0);
        }

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
                    // std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
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

        // }

        
    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < root->nucMutation.size(); i++) {
        int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

        // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
        int32_t nucPosition = root->nucMutation[i].nucPosition;
        int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
        uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            // Either S, I or D
            int len = ((root->nucMutation[i].mutInfo) >> 4);

            if(primaryBlockId >= sequence.size()) {
                std::cout << primaryBlockId << " " << sequence.size() << std::endl;
            }

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
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
        } else {
            if(type == panmanUtils::NucMutationType::NSNPS) {
                // SNP Substitution
                newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
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
            } else if(type == panmanUtils::NucMutationType::NSNPI) {
                // SNP Insertion
                newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
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
            } else if(type == panmanUtils::NucMutationType::NSNPD) {
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
    // }

    if(root->children.size() == 0 || rootSeq) {
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
            // std::cout << "rotating" << std::endl;
            rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
            rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
            rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
        }

        if(sequenceInverted.find(root->identifier) != sequenceInverted.end() && sequenceInverted[root->identifier]) {
            // std::cout << "inverting" << std::endl;
            reverse(sequencePrint.begin(), sequencePrint.end());
            reverse(blockExistsPrint.begin(), blockExistsPrint.end());
            reverse(blockStrandPrint.begin(), blockStrandPrint.end());
        }
        if (allIndex) {
            // bool* checkA;
            // bool* checkB;
            // *checkA = false;
            // *checkB = false;
            // int startCoordinate = getUnalignedGlobalCoordinate(std::get<0>(panMATStart),
            //                                            std::get<1>(panMATStart),
            //                                            std::get<2>(panMATStart),
            //                                            std::get<3>(panMATStart),
            //                                            sequencePrint,
            //                                            blockExistsPrint,
            //                                            blockStrandPrint,
            //                                            circularSequences[root->identifier],
            //                                            checkA
            //                                         );
    
            // int endCoordinate = getUnalignedGlobalCoordinate(std::get<0>(panMATEnd),
            //                                                 std::get<1>(panMATEnd),
            //                                                 std::get<2>(panMATEnd),
            //                                                 std::get<3>(panMATEnd),
            //                                                 sequencePrint,
            //                                                 blockExistsPrint,
            //                                                 blockStrandPrint,
            //                                                 circularSequences[root->identifier],
            //                                                 checkB
            //                                             );
            
            // if (checkA) {
            //     startCoordinate = -1;
            // }
            // if (checkB) {
            //     endCoordinate = -1;
            // }
            // std::cout << root->identifier << " " << startCoordinate << " " << endCoordinate << " offsets " << circularSequences[root->identifier] << " " << offset << std::endl;
            // std::cout << "printFASTA start" << std::get<0>(panMATStart) << " " << std::get<1>(panMATStart) << " " << std::get<2>(panMATStart) << " " << std::get<3>(panMATStart) << std::endl;
            // std::cout << "printFASTA end" << std::get<0>(panMATEnd) << " " << std::get<1>(panMATEnd) << " " << std::get<2>(panMATEnd) << " " << std::get<3>(panMATEnd) << std::endl;
            panmanUtils::printSubsequenceLines(sequencePrint, blockExistsPrint, blockStrandPrint, 70, panMATStart, panMATEnd, aligned, fout, offset);
        } else {
            panmanUtils::printSequenceLines(sequencePrint, blockExistsPrint, blockStrandPrint, 70, aligned, fout, offset);
        }
    } else {

        // DFS on children
        for(panmanUtils::Node* child: root->children) {
            printFASTAHelper(child, sequence, blockExists, blockStrand, fout, aligned, rootSeq, panMATStart, panMATEnd, allIndex);

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

    // std::cout << "Done iteration for node: " << root->identifier << std::endl; 

}

void panmanUtils::Tree::printSingleNodeHelper(std::vector<panmanUtils::Node*> &nodeList, int nodeListIndex, sequence_t& sequence,
        blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd) {
    
    panmanUtils::Node* node = nodeList[nodeListIndex--];

    // Apply mutations
    // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
    std::vector< std::tuple< int32_t, int32_t, bool, bool, bool, bool > > blockMutationInfo;

    // Block Mutations
    for(auto mutation: node->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        int32_t secondaryBlockId = mutation.secondaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        if (secondaryBlockId != -1) {
            std::cout << "Error: Block Secondary ID is not -1" << std::endl;
            exit(0);
        }

        // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
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
                    // std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
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

        // }

        
    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int32_t, int, int, char, char > > mutationInfo;

    // Nuc mutations
    for(size_t i = 0; i < node->nucMutation.size(); i++) {
        int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

        // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            // Either S, I or D
            int len = ((node->nucMutation[i].mutInfo) >> 4);

            if(primaryBlockId >= sequence.size()) {
                std::cout << primaryBlockId << " " << sequence.size() << std::endl;
            }

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                if(secondaryBlockId != -1) {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }

                    }
                } else {
                    if(nucGapPosition != -1) {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j];
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        }
                    } else {
                        for(int j = 0; j < len; j++) {
                            char oldVal = sequence[primaryBlockId].first[nucPosition+j].first;
                            newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                            sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
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
        } else {
            if(type == panmanUtils::NucMutationType::NSNPS) {
                // SNP Substitution
                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
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
            } else if(type == panmanUtils::NucMutationType::NSNPI) {
                // SNP Insertion
                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
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
            } else if(type == panmanUtils::NucMutationType::NSNPD) {
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
    // }

    if(nodeListIndex < 0) {
        // Print sequence
        fout << '>' << node->identifier << std::endl;

        int offset = 0;
        if(!aligned && circularSequences.find(node->identifier) != circularSequences.end()) {
            // If MSA is to be printed, offset doesn't matter
            offset = circularSequences[node->identifier];
        }
        sequence_t sequencePrint = sequence;
        blockExists_t blockExistsPrint = blockExists;
        blockStrand_t blockStrandPrint = blockStrand;

        if(rotationIndexes.find(node->identifier) != rotationIndexes.end() && rotationIndexes[node->identifier] != 0) {
            int ctr = -1, rotInd = 0;
            for(size_t i = 0; i < blockExistsPrint.size(); i++) {
                if(blockExistsPrint[i].first) {
                    ctr++;
                }
                if(ctr == rotationIndexes[node->identifier]) {
                    rotInd = i;
                    break;
                }
            }
            rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
            rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
            rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
        }

        if(sequenceInverted.find(node->identifier) != sequenceInverted.end() && sequenceInverted[node->identifier]) {
            reverse(sequencePrint.begin(), sequencePrint.end());
            reverse(blockExistsPrint.begin(), blockExistsPrint.end());
            reverse(blockStrandPrint.begin(), blockStrandPrint.end());
        }

        panmanUtils::printSubsequenceLines(sequencePrint, blockExistsPrint, blockStrandPrint, 70, panMATStart, panMATEnd, aligned, fout, offset);
    } else {
            printSingleNodeHelper(nodeList, nodeListIndex, sequence, blockExists, blockStrand, fout, aligned, rootSeq, panMATStart, panMATEnd);
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

void panmanUtils::Tree::printFASTA(std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {
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
    
    printFASTAHelper(root, sequence, blockExists, blockStrand, fout, aligned, rootSeq, panMATStart, panMATEnd, allIndex);

}

void panmanUtils::Tree::printSingleNode(std::ostream& fout, const sequence_t& sequenceRef,
                                         const blockExists_t& blockExistsRef, const blockStrand_t& blockStrandRef,
                                         std::string nodeIdentifier, std::tuple< int, int, int, int >& panMATStart, std::tuple< int, int, int, int >& panMATEnd) {
    // List nodes from root to nodeIdentifier
    std::vector<Node*> nodeList;
    Node* newNode = allNodes[nodeIdentifier];
    nodeList.push_back(newNode);
    while (newNode->parent != nullptr) {
        newNode = newNode->parent;
        nodeList.push_back(newNode);
    }

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
    // Run traversal on nodeList to extract sequences
    printSingleNodeHelper(nodeList, (nodeList.size()-1), sequence, blockExists, blockStrand, fout, false, false, panMATStart, panMATEnd);

}

void panmanUtils::Tree::printFASTAParallel(std::ostream& fout, bool aligned) {

    std::mutex fastaMutex;
    size_t lineSize = 70;
    std::cout << tbb::this_task_arena::max_concurrency() << std::endl;

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

// void panmanUtils::Tree::printFASTAHelperNew(panmanUtils::Node* root, 
//                             std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
//                             std::vector<bool>& blockExists, 
//                             std::vector<bool>& blockStrand, std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {
        
//     // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
//     std::vector< std::tuple< int32_t, bool, bool, bool, bool > > blockMutationInfo;

//     // Block Mutations
//     for(auto mutation: root->blockMutation) {
//         int32_t primaryBlockId = mutation.primaryBlockId;
//         bool type = mutation.blockMutInfo;
//         bool inversion = mutation.inversion;

//         if(type == 1) {
//             // insertion

//             bool oldStrand;
//             bool oldMut;
//             oldStrand = blockStrand[primaryBlockId];
//             oldMut = blockExists[primaryBlockId];
//             blockExists[primaryBlockId] = true;

//             // if insertion of inverted block takes place, the strand is backwards
//             blockStrand[primaryBlockId] = !inversion;
//             blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, true, !inversion) );


//         } else {
//             bool oldMut;
//             bool oldStrand;
//             if(inversion) {
//                 // This means that this is not a deletion, but instead an inversion
//                 oldStrand = blockStrand[primaryBlockId];
//                 oldMut = blockExists[primaryBlockId];
//                 blockStrand[primaryBlockId] = !oldStrand;
                
//                 if(oldMut != true) {
//                     std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
//                 }
//                 blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
//             } else {
//                 // Actually a deletion
//                 oldStrand = blockStrand[primaryBlockId];
//                 oldMut = blockExists[primaryBlockId];
//                 blockExists[primaryBlockId] = false;

//                 // resetting strand to true during deletion
//                 blockStrand[primaryBlockId] = true;
//             }
//             blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, false, true) );

//         }

//         // }

        
//     }

//     // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
//     std::vector< std::tuple< int32_t, int, int, char, char > > mutationInfo;

//     // Nuc mutations
//     for(size_t i = 0; i < root->nucMutation.size(); i++) {
//         int32_t primaryBlockId = root->nucMutation[i].primaryBlockId;
//         int32_t secondaryBlockId = root->nucMutation[i].secondaryBlockId;

//         // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
//         int32_t nucPosition = root->nucMutation[i].nucPosition;
//         int32_t nucGapPosition = root->nucMutation[i].nucGapPosition;
//         uint32_t type = (root->nucMutation[i].mutInfo & 0x7);
//         char newVal = '-';

//         if(type < 3) {
//             // Either S, I or D
//             int len = ((root->nucMutation[i].mutInfo) >> 4);

//             if(primaryBlockId >= sequence.size()) {
//                 std::cout << primaryBlockId << " " << sequence.size() << std::endl;
//             }

//             if(type == panmanUtils::NucMutationType::NS) {
//                 // Substitution
//                 if(nucGapPosition != -1) {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                         newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
//                     }
//                 } else {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                         newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                         sequence[primaryBlockId][nucPosition+j].first = newVal;
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
//                     }
//                 }
//             } else if(type == panmanUtils::NucMutationType::NI) {
//                 // Insertion
//                 if(nucGapPosition != -1) {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                         newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
//                     }
//                 } else {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                         newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                         sequence[primaryBlockId][nucPosition+j].first = newVal;
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
//                     }
//                 }
//             } else if(type == panmanUtils::NucMutationType::ND) {
//                 // Deletion
//                 if(nucGapPosition != -1) {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
//                     }
//                 } else {
//                     for(int j = 0; j < len; j++) {
//                         char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                         sequence[primaryBlockId][nucPosition+j].first = '-';
//                         mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
//                     }
//                 }
//             }
//         } else {
//             if(type == panmanUtils::NucMutationType::NSNPS) {
//                 // SNP Substitution
//                 newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
//                 if(nucGapPosition != -1) {
//                     char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                     sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
//                 } else {
//                     char oldVal = sequence[primaryBlockId][nucPosition].first;
//                     sequence[primaryBlockId][nucPosition].first = newVal;
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
//                 }
//             } else if(type == panmanUtils::NucMutationType::NSNPI) {
//                 // SNP Insertion
//                 newVal = panmanUtils::getNucleotideFromCode(((root->nucMutation[i].nucs) >> 20) & 0xF);
//                 if(nucGapPosition != -1) {
//                     char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                     sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
//                 } else {
//                     char oldVal = sequence[primaryBlockId][nucPosition].first;
//                     sequence[primaryBlockId][nucPosition].first = newVal;
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
//                 }
//             } else if(type == panmanUtils::NucMutationType::NSNPD) {
//                 // SNP Deletion
//                 if(nucGapPosition != -1) {
//                     char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                     sequence[primaryBlockId][nucPosition].second[nucGapPosition] = '-';
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
//                 } else {
//                     char oldVal = sequence[primaryBlockId][nucPosition].first;
//                     sequence[primaryBlockId][nucPosition].first = '-';
//                     mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
//                 }
//             }
//         }
//     }
//     // }

//     if(root->children.size() == 0 || rootSeq) {
//         // Print sequence

//         fout << '>' << root->identifier << std::endl;

//         int offset = 0;
//         if(!aligned && circularSequences.find(root->identifier) != circularSequences.end()) {
//             // If MSA is to be printed, offset doesn't matter
//             offset = circularSequences[root->identifier];
//         }
//         std::vector< std::vector< std::pair< char, std::vector< char > > > > sequencePrint = sequence;
//         std::vector< bool > blockExistsPrint = blockExists;
//         std::vector< bool > blockStrandPrint = blockStrand;

//         if(rotationIndexes.find(root->identifier) != rotationIndexes.end() && rotationIndexes[root->identifier] != 0) {
//             int ctr = -1, rotInd = 0;
//             for(size_t i = 0; i < blockExistsPrint.size(); i++) {
//                 if(blockExistsPrint[i]) {
//                     ctr++;
//                 }
//                 if(ctr == rotationIndexes[root->identifier]) {
//                     rotInd = i;
//                     break;
//                 }
//             }
//             // std::cout << "rotating" << std::endl;
//             rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
//             rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
//             rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
//         }

//         if(sequenceInverted.find(root->identifier) != sequenceInverted.end() && sequenceInverted[root->identifier]) {
//             // std::cout << "inverting" << std::endl;
//             reverse(sequencePrint.begin(), sequencePrint.end());
//             reverse(blockExistsPrint.begin(), blockExistsPrint.end());
//             reverse(blockStrandPrint.begin(), blockStrandPrint.end());
//         }
        
//          panmanUtils::printSequenceLinesNew(sequencePrint, blockExistsPrint, blockStrandPrint, 70, aligned, fout, offset);
//     } else {

//         // DFS on children
//         for(panmanUtils::Node* child: root->children) {
//             printFASTAHelperNew(child, sequence, blockExists, blockStrand, fout, aligned, rootSeq, panMATStart, panMATEnd, allIndex);

//         }
//     }


//     // Undo block mutations when current node and its subtree have been processed
//     for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
//         auto mutation = *it;
//         blockExists[std::get<0>(mutation)] = std::get<1>(mutation);
//         blockStrand[std::get<0>(mutation)] = std::get<2>(mutation);

//     }

//     // Undo nuc mutations when current node and its subtree have been processed
//     for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
//         auto mutation = *it;
//         if(std::get<2>(mutation) != -1) {
//             sequence[std::get<0>(mutation)][std::get<1>(mutation)].second[std::get<2>(mutation)] = std::get<3>(mutation);
//         } else {
//             sequence[std::get<0>(mutation)][std::get<1>(mutation)].first = std::get<3>(mutation);
//         }

//     }

//     // std::cout << "Done iteration for node: " << root->identifier << std::endl; 
// }


// void panmanUtils::Tree::printFASTANew(std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {
//     // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
//     std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence(blocks.size() + 1);
//     std::vector< bool >  blockExists(blocks.size() + 1, false, {});
//     std::vector< bool >  blockStrand(blocks.size() + 1, true, {});


//     int32_t maxBlockId = 0;

//     // Create consensus sequence of blocks
//     for(size_t i = 0; i < blocks.size(); i++) {

//         int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);

//         maxBlockId = std::max(maxBlockId, primaryBlockId);

//         for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
//             bool endFlag = false;
//             for(size_t k = 0; k < 8; k++) {
//                 const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

//                 if(nucCode == 0) {
//                     endFlag = true;
//                     break;
//                 }
//                 const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

//                 sequence[primaryBlockId].push_back({nucleotide, {}});
//             }

//             if(endFlag) {
//                 break;
//             }
//         }

//         // End character to incorporate for gaps at the end
//         sequence[primaryBlockId].push_back({'x', {}});
//     }

//     sequence.resize(maxBlockId + 1);
//     blockExists.resize(maxBlockId + 1);
//     blockStrand.resize(maxBlockId + 1);

//     // Assigning nucleotide gaps in blocks
//     for(size_t i = 0; i < gaps.size(); i++) {
//         int32_t primaryBId = (gaps[i].primaryBlockId);
//         int32_t secondaryBId = (gaps[i].secondaryBlockId);

//         for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
//             int len = gaps[i].nucGapLength[j];
//             int pos = gaps[i].nucPosition[j];
//             sequence[primaryBId][pos].second.resize(len, '-');
//         }
//     }

//     // Run depth first traversal to extract sequences
    
//     printFASTAHelperNew(root, sequence, blockExists, blockStrand, fout, aligned, rootSeq, panMATStart, panMATEnd, allIndex);

// }

void getNodesFromTipToRoot(panmanUtils::Node* node, std::vector<panmanUtils::Node*> &nodesFromTipToRoot){
    panmanUtils::Node* currentNode = node;
    while(currentNode->parent != nullptr){
        nodesFromTipToRoot.push_back(currentNode);
        currentNode = currentNode->parent;
    }
    nodesFromTipToRoot.push_back(currentNode);

    reverse(nodesFromTipToRoot.begin(), nodesFromTipToRoot.end());

    return;
}

void getBlockSequence(std::vector<panmanUtils::Node*> &nodesFromTipToRoot, 
                      std::vector< bool >& blockExists){
    // panmanUtils::Node* node;
    for (auto node: nodesFromTipToRoot){ 
        // node = nodesFromTipToRoot[i];
        for(auto mutation: node->blockMutation) {
            int32_t primaryBlockId = mutation.primaryBlockId;
            bool type = mutation.blockMutInfo;
            bool inversion = mutation.inversion;
            if(type == 1) {
                // insertion
                blockExists[primaryBlockId] = true;
            } else {
                // deletion
                if(!inversion) {
                    blockExists[primaryBlockId] = false;
                }
            }
        }
    }

}

std::string panmanUtils::Tree::printFASTAUltraFastHelper(
                            const std::vector<bool>& blockSequence,
                            std::unordered_map<int, int>& blockLengths,
                            const std::vector<panmanUtils::Node*>& nodesFromTipToRoot, 
                            std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                            std::vector<bool>& blockExists, 
                            std::vector<bool>& blockStrand, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {

    
    for (auto node: nodesFromTipToRoot){
        // Block Mutations
        for(auto mutation: node->blockMutation) {
            int32_t primaryBlockId = mutation.primaryBlockId;
            bool type = mutation.blockMutInfo;
            bool inversion = mutation.inversion;
            if (blockSequence[primaryBlockId]) {
                if(type == 1) {
                    // insertion
                    bool oldStrand;
                    bool oldMut;
                    oldStrand = blockStrand[primaryBlockId];
                    oldMut = blockExists[primaryBlockId];
                    blockExists[primaryBlockId] = true;
                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId] = !inversion;
                } else {
                    bool oldMut;
                    bool oldStrand;
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        oldStrand = blockStrand[primaryBlockId];
                        oldMut = blockExists[primaryBlockId];
                        blockStrand[primaryBlockId] = !oldStrand;
                        if(oldMut != true) {
                            // std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
                        }
                    } else {
                        // Actually a deletion
                        oldStrand = blockStrand[primaryBlockId];
                        oldMut = blockExists[primaryBlockId];
                        blockExists[primaryBlockId] = false;
                        // resetting strand to true during deletion
                        blockStrand[primaryBlockId] = true;
                    }
                }
            }
        }

        // Nuc mutations
        for(size_t i = 0; i < node->nucMutation.size(); i++) {
            int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

            if (blockSequence[primaryBlockId]) {
            // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
                int32_t nucPosition = node->nucMutation[i].nucPosition;
                int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
                uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';

                if(type < 3) {
                    // Either S, I or D
                    int len = ((node->nucMutation[i].mutInfo) >> 4);

                    if(primaryBlockId >= sequence.size()) {
                        std::cout << primaryBlockId << " " << sequence.size() << std::endl;
                    }

                    if(type == panmanUtils::NucMutationType::NS) {
                        // Substitution
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition+j].first = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        // Insertion
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition+j].first = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        // Deletion
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                sequence[primaryBlockId][nucPosition+j].first = '-';
                            }
                        }
                    }
                } else {
                    if(type == panmanUtils::NucMutationType::NSNPS) {
                        // SNP Substitution
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = newVal;
                        }
                    } else if(type == panmanUtils::NucMutationType::NSNPI) {
                        // SNP Insertion
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = newVal;
                        }
                    } else if(type == panmanUtils::NucMutationType::NSNPD) {
                        // SNP Deletion
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    
    // Store sequence
    panmanUtils::Node* tipNode = nodesFromTipToRoot[nodesFromTipToRoot.size()-1];
    std::string line="";
    line += '>' + tipNode->identifier + '\n';
    

    int offset = 0;
    if(!aligned && circularSequences.find(tipNode->identifier) != circularSequences.end()) {
        // If MSA is to be printed, offset doesn't matter
        offset = circularSequences[tipNode->identifier];
    }
    std::vector< std::vector< std::pair< char, std::vector< char > > > > sequencePrint = sequence;
    std::vector< bool > blockExistsPrint = blockExists;
    std::vector< bool > blockStrandPrint = blockStrand;

    if(rotationIndexes.find(tipNode->identifier) != rotationIndexes.end() && rotationIndexes[tipNode->identifier] != 0) {
        int ctr = -1, rotInd = 0;
        for(size_t i = 0; i < blockExistsPrint.size(); i++) {
            if(blockExistsPrint[i]) {
                ctr++;
            }
            if(ctr == rotationIndexes[tipNode->identifier]) {
                rotInd = i;
                break;
            }
        }
        // std::cout << "rotating" << std::endl;
        rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
        rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
        rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
    }

    if(sequenceInverted.find(tipNode->identifier) != sequenceInverted.end() && sequenceInverted[tipNode->identifier]) {
        // std::cout << "inverting" << std::endl;
        reverse(sequencePrint.begin(), sequencePrint.end());
        reverse(blockExistsPrint.begin(), blockExistsPrint.end());
        reverse(blockStrandPrint.begin(), blockStrandPrint.end());
    }
    
    line += panmanUtils::printSequenceLinesNew(sequencePrint, blockLengths, blockExistsPrint, blockStrandPrint, 70, aligned, offset, false);
    return line;
}

void panmanUtils::Tree::printFASTAUltraFast(std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {

    std::unordered_map< std::string, std::mutex > nodeMutexes;
    for(auto u: allNodes) {
        if (u.second->children.size() == 0) {
            nodeMutexes[u.first];
        }
    }

    std::mutex printMutex;
	int counting=0;
    // for (auto &keyValue: allNodes) {
    tbb::parallel_for_each(allNodes.begin(), allNodes.end(), [&](const std::pair<std::string, panmanUtils::Node*>& keyValue) {
        panmanUtils::Node* node = keyValue.second;
        // Create a stringstream for each thread to avoid race conditions on fout
        if (node->children.size() != 0) {
            return;
        }

        // Get block sequnece of the Tip
        std::vector< bool >  blockSequence(blocks.size() + 1, false, {});
        std::vector<panmanUtils::Node*> nodesFromTipToRoot;
        getNodesFromTipToRoot(node, nodesFromTipToRoot);
        getBlockSequence(nodesFromTipToRoot, blockSequence);

        // Blocks length
        std::unordered_map<int, int> blockLengths;

        // Expanding blocks only if exist in tip 
        std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence(blocks.size() + 1);
        std::vector< bool >  blockExists(blocks.size() + 1, false, {});
        std::vector< bool >  blockStrand(blocks.size() + 1, true, {});


        int32_t maxBlockId = 0;

        // Create consensus sequence of blocks
        for(size_t i = 0; i < blocks.size(); i++) {
            int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
            blockLengths[primaryBlockId] = 0;
            maxBlockId = std::max(maxBlockId, primaryBlockId);
            if (blockSequence[primaryBlockId]) {
                // int len = 0;
                for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
                    bool endFlag = false;
                    for(size_t k = 0; k < 8; k++) {
                        const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                        if(nucCode == 0) {
                            endFlag = true;
                            break;
                        }
                        // len++;
                        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
                        sequence[primaryBlockId].push_back({nucleotide, {}});
                    }

                    if(endFlag) {
                        break;
                    }
                }
                // End character to incorporate for gaps at the end
                sequence[primaryBlockId].push_back({'x', {}});
                // blockLengths[primaryBlockId] += len;
            } else {
                int len = 0;
                for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
                    bool endFlag = false;
                    for(size_t k = 0; k < 8; k++) {
                        const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                        if(nucCode == 0) {
                            endFlag = true;
                            break;
                        }
                        len++;
                        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
                    }

                    if(endFlag) {
                        break;
                    }
                }
                blockLengths[primaryBlockId] += len;
            }
        }

        sequence.resize(maxBlockId + 1);
        blockExists.resize(maxBlockId + 1);
        blockStrand.resize(maxBlockId + 1);

        // Assigning nucleotide gaps in blocks
        for(size_t i = 0; i < gaps.size(); i++) {
            int32_t primaryBId = (gaps[i].primaryBlockId);
            int32_t secondaryBId = (gaps[i].secondaryBlockId);
            if (blockSequence[primaryBId]){
                for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                    int len = gaps[i].nucGapLength[j];
                    int pos = gaps[i].nucPosition[j];
                    sequence[primaryBId][pos].second.resize(len, '-');
                    // blockLengths[primaryBId] += len;
                }
            } else {
                int len=0;
                for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                    len += gaps[i].nucGapLength[j];
                }
                blockLengths[primaryBId] += len;
            }
        }

        std::string line = printFASTAUltraFastHelper(blockSequence, blockLengths, nodesFromTipToRoot, sequence, blockExists, blockStrand, aligned, rootSeq, panMATStart, panMATEnd, allIndex);
        // nodeMutexes[node->identifier].lock();
        std::lock_guard<std::mutex> guard(printMutex);
        fout << line << "\n";
        // nodeMutexes[node->identifier].unlock();
        // break;
    });
    // }
}

std::pair<std::vector<std::string>, std::vector<int>> panmanUtils::Tree::extractSequenceHelper(
                            const std::vector<bool>& blockSequence,
                            std::unordered_map<int, int>& blockLengths,
                            const std::vector<panmanUtils::Node*>& nodesFromTipToRootIn, 
                            std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                            std::vector<bool>& blockExists, 
                            std::vector<bool>& blockStrand, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {

    // reverse traversal of nodes
    std::vector<panmanUtils::Node*> nodesFromTipToRoot(nodesFromTipToRootIn.size());
    for (int i = 0; i < nodesFromTipToRoot.size(); i++){
        nodesFromTipToRoot[i] = nodesFromTipToRootIn[nodesFromTipToRootIn.size() - 1- i];
    }

    for (auto node: nodesFromTipToRoot){
        // Block Mutations
        for(auto mutation: node->blockMutation) {
            int32_t primaryBlockId = mutation.primaryBlockId;
            bool type = mutation.blockMutInfo;
            bool inversion = mutation.inversion;
            if (blockSequence[primaryBlockId]) {
                if(type == 1) {
                    // insertion
                    bool oldStrand;
                    bool oldMut;
                    oldStrand = blockStrand[primaryBlockId];
                    oldMut = blockExists[primaryBlockId];
                    blockExists[primaryBlockId] = true;
                    // if insertion of inverted block takes place, the strand is backwards
                    blockStrand[primaryBlockId] = !inversion;
                } else {
                    bool oldMut;
                    bool oldStrand;
                    if(inversion) {
                        // This means that this is not a deletion, but instead an inversion
                        oldStrand = blockStrand[primaryBlockId];
                        oldMut = blockExists[primaryBlockId];
                        blockStrand[primaryBlockId] = !oldStrand;
                        if(oldMut != true) {
                            // std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
                        }
                    } else {
                        // Actually a deletion
                        oldStrand = blockStrand[primaryBlockId];
                        oldMut = blockExists[primaryBlockId];
                        blockExists[primaryBlockId] = false;
                        // resetting strand to true during deletion
                        blockStrand[primaryBlockId] = true;
                    }
                }
            }
        }

        // Nuc mutations
        for(size_t i = 0; i < node->nucMutation.size(); i++) {
            int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

            if (blockSequence[primaryBlockId]) {
            // if (rootSeq && (primaryBlockId>=std::get<0>(panMATStart) && primaryBlockId<=std::get<0>(panMATEnd)) && (secondaryBlockId<=std::get<1>(panMATStart) && secondaryBlockId<=std::get<1>(panMATEnd)) ) {
                int32_t nucPosition = node->nucMutation[i].nucPosition;
                int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
                uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
                char newVal = '-';

                if(type < 3) {
                    // Either S, I or D
                    int len = ((node->nucMutation[i].mutInfo) >> 4);

                    if(primaryBlockId >= sequence.size()) {
                        std::cout << primaryBlockId << " " << sequence.size() << std::endl;
                    }

                    if(type == panmanUtils::NucMutationType::NS) {
                        // Substitution
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition+j].first = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::NI) {
                        // Insertion
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                sequence[primaryBlockId][nucPosition+j].first = newVal;
                            }
                        }
                    } else if(type == panmanUtils::NucMutationType::ND) {
                        // Deletion
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                                sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                                sequence[primaryBlockId][nucPosition+j].first = '-';
                            }
                        }
                    }
                } else {
                    if(type == panmanUtils::NucMutationType::NSNPS) {
                        // SNP Substitution
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = newVal;
                        }
                    } else if(type == panmanUtils::NucMutationType::NSNPI) {
                        // SNP Insertion
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = newVal;
                        }
                    } else if(type == panmanUtils::NucMutationType::NSNPD) {
                        // SNP Deletion
                        if(nucGapPosition != -1) {
                            char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                            sequence[primaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            char oldVal = sequence[primaryBlockId][nucPosition].first;
                            sequence[primaryBlockId][nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    
    // Store sequence
    panmanUtils::Node* tipNode = nodesFromTipToRoot[nodesFromTipToRoot.size()-1];

    int offset = 0;
    if(!aligned && circularSequences.find(tipNode->identifier) != circularSequences.end()) {
        // If MSA is to be printed, offset doesn't matter
        offset = circularSequences[tipNode->identifier];
    }
    std::vector< std::vector< std::pair< char, std::vector< char > > > > sequencePrint = sequence;
    std::vector< bool > blockExistsPrint = blockExists;
    std::vector< bool > blockStrandPrint = blockStrand;

    if(rotationIndexes.find(tipNode->identifier) != rotationIndexes.end() && rotationIndexes[tipNode->identifier] != 0) {
        int ctr = -1, rotInd = 0;
        for(size_t i = 0; i < blockExistsPrint.size(); i++) {
            if(blockExistsPrint[i]) {
                ctr++;
            }
            if(ctr == rotationIndexes[tipNode->identifier]) {
                rotInd = i;
                break;
            }
        }
        // std::cout << "rotating" << std::endl;
        rotate(sequencePrint.begin(), sequencePrint.begin() + rotInd, sequencePrint.end());
        rotate(blockExistsPrint.begin(), blockExistsPrint.begin() + rotInd, blockExistsPrint.end());
        rotate(blockStrandPrint.begin(), blockStrandPrint.begin() + rotInd, blockStrandPrint.end());
    }

    if(sequenceInverted.find(tipNode->identifier) != sequenceInverted.end() && sequenceInverted[tipNode->identifier]) {
        // std::cout << "inverting" << std::endl;
        reverse(sequencePrint.begin(), sequencePrint.end());
        reverse(blockExistsPrint.begin(), blockExistsPrint.end());
        reverse(blockStrandPrint.begin(), blockStrandPrint.end());
    }
    
    return panmanUtils::printSequenceLinesNewer(sequencePrint, blockLengths, blockExistsPrint, blockStrandPrint, 70, aligned, offset, false);
}

std::pair<std::vector<std::string>, std::vector<int>> panmanUtils::Tree::extractSingleSequence(panmanUtils::Node* node, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex) {

    // Create a stringstream for each thread to avoid race conditions on fout
    if (node->children.size() != 0) {
        std::cerr << "Node is not a tip" << std::endl;
        exit(0);
    }
    // Get block sequnece of the Tip
    std::vector< bool >  blockSequence(blocks.size() + 1, false, {});
    std::vector<panmanUtils::Node*> nodesFromTipToRoot;
    getNodesFromTipToRoot(node, nodesFromTipToRoot);
    getBlockSequence(nodesFromTipToRoot, blockSequence);

    // Blocks length
    std::unordered_map<int, int> blockLengths;

    // Expanding blocks only if exist in tip 
    std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence(blocks.size() + 1);
    std::vector< bool >  blockExists(blocks.size() + 1, false, {});
    std::vector< bool >  blockStrand(blocks.size() + 1, true, {});


    int32_t maxBlockId = 0;

    // Create consensus sequence of blocks
    for(size_t i = 0; i < blocks.size(); i++) {
        int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
        blockLengths[primaryBlockId] = 0;
        maxBlockId = std::max(maxBlockId, primaryBlockId);
        if (blockSequence[i]) {
            int len = 0;
            for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
                bool endFlag = false;
                for(size_t k = 0; k < 8; k++) {
                    const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

                    if(nucCode == 0) {
                        endFlag = true;
                        break;
                    }
                    len++;
                    const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

                    sequence[primaryBlockId].push_back({nucleotide, {}});
                }

                if(endFlag) {
                    break;
                }
            }
            // End character to incorporate for gaps at the end
            sequence[primaryBlockId].push_back({'x', {}});
            blockLengths[primaryBlockId] += len;
        } else {
            int len = 0;
            for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
                bool endFlag = false;
                for(size_t k = 0; k < 8; k++) {
                    const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                    if(nucCode == 0) {
                        endFlag = true;
                        break;
                    }
                    len++;
                }

                if(endFlag) {
                    break;
                }
            }
            blockLengths[primaryBlockId] += len;
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    // Assigning nucleotide gaps in blocks
    for(size_t i = 0; i < gaps.size(); i++) {
        int32_t primaryBId = (gaps[i].primaryBlockId);
        int32_t secondaryBId = (gaps[i].secondaryBlockId);
        if (blockSequence[primaryBId]){
            for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                int len = gaps[i].nucGapLength[j];
                int pos = gaps[i].nucPosition[j];
                sequence[primaryBId][pos].second.resize(len, '-');
                blockLengths[primaryBId] += len;
            }
        } else {
            int len=0;
            for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                len += gaps[i].nucGapLength[j];
            }
            blockLengths[primaryBId] += len;
        }
    }

    return extractSequenceHelper(blockSequence, blockLengths, nodesFromTipToRoot, sequence, blockExists, blockStrand, aligned, rootSeq, panMATStart, panMATEnd, allIndex);
}
