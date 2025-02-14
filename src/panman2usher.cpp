#include "panmanUtils.hpp"

void getCoordMap(panmanUtils::Tree* panmanTree, std::vector<std::vector<std::pair<int, std::vector<int>>>> &globalCoords_t) {
    const std::vector<panmanUtils::Block> &blocks = panmanTree->blocks;
    const std::vector<panmanUtils::GapList> &gaps = panmanTree->gaps;

    // blocks
    for (size_t block_id = 0; block_id < blocks.size(); block_id++) {
        int32_t blockId = ((int32_t)blocks[block_id].primaryBlockId);
        for (size_t nuc_pos = 0; nuc_pos < blocks[block_id].consensusSeq.size(); nuc_pos++) {
            bool endFlag = false;
            for (size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[block_id].consensusSeq[nuc_pos]) >> (4 * (7 - k))) & 15);
                if (nucCode == 0) {
                    endFlag = true;
                    break;
                }
                globalCoords_t[blockId].push_back({0, {}});
            }
            if (endFlag){
                break;
            }
        }
        globalCoords_t[blockId].push_back({0, {}}); // do I need this?
    }

    // nuc gaps
    for (size_t i = 0; i < gaps.size(); i++) {
        int32_t blockId = (gaps[i].primaryBlockId);
        for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            globalCoords_t[blockId][pos].second.resize(len, 0);
        }
    }
    
    // Assigning coordinates
    int index = 1;
    for (size_t blockId = 0; blockId < globalCoords_t.size(); blockId++){
        for (size_t j = 0; j < globalCoords_t[blockId].size(); j++){
            for (size_t k = 0; k < globalCoords_t[blockId][j].second.size(); k++){
                globalCoords_t[blockId][j].second[k] = index;
                index++;
            }
            globalCoords_t[blockId][j].first = index;
            index++;
        }
    }
    return;

}

void getPseudoRoot(panmanUtils::Tree* panmanTree, std::vector<std::vector<std::pair<char, std::vector<char>>>> &pseudoRoot) {
    const std::vector<panmanUtils::Block> &blocks = panmanTree->blocks;
    const std::vector<panmanUtils::GapList> &gaps = panmanTree->gaps;

    // blocks
    for (size_t block_id = 0; block_id < blocks.size(); block_id++) {
        int32_t blockId = ((int32_t)blocks[block_id].primaryBlockId);
        for (size_t nuc_pos = 0; nuc_pos < blocks[block_id].consensusSeq.size(); nuc_pos++) {
            bool endFlag = false;
            for (size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[block_id].consensusSeq[nuc_pos]) >> (4 * (7 - k))) & 15);
                if (nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
                pseudoRoot[blockId].push_back({nucleotide, {}});
            }
            if (endFlag){
                break;
            }
        }
        pseudoRoot[blockId].push_back({'x', {}}); // do I need this?
    }

    // nuc gaps
    for (size_t i = 0; i < gaps.size(); i++) {
        int32_t blockId = (gaps[i].primaryBlockId);
        for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            pseudoRoot[blockId][pos].second.resize(len, '-');
        }
    }
    
    return;

}

int8_t get_nuc_id (char nuc) {
    int8_t ret = 0b1111;
    switch(nuc) {
    case 'a':
    case 'A':
        ret = 0b1;
        break;
    case 'c':
    case 'C':
        ret = 0b10;
        break;
    case 'g':
    case 'G':
        ret = 0b100;
        break;
    case 't':
    case 'T':
        ret = 0b1000;
        break;
    case 'R':
        ret = 0b101;
        break;
    case 'Y':
        ret = 0b1010;
        break;
    case 'S':
        ret = 0b110;
        break;
    case 'W':
        ret = 0b1001;
        break;
    case 'K':
        ret = 0b1100;
        break;
    case 'M':
        ret = 0b11;
        break;
    case 'B':
        ret = 0b1110;
        break;
    case 'D':
        ret = 0b1101;
        break;
    case 'H':
        ret = 0b1011;
        break;
    case 'V':
        ret = 0b111;
    case 'n':
    case 'N':
    default:
        ret = 0b1111;
        break;
    }
    return ret;
}

// Sets bits at positions specified by nuc_vec to 1 in int8
int8_t get_nuc_id (std::vector<int8_t> nuc_vec) {
    int8_t ret = 0;
    int8_t one = 1;
    for (auto nuc: nuc_vec) {
        assert((nuc >= 0) && (nuc <=3));
        ret += (one << nuc);
    }
    return ret;
}

// Convert nuc_id back to IUPAC base
char get_nuc (int8_t nuc_id) {
    char ret = 'N';
    //assert ((nuc_id >= 1) && (nuc_id <= 15));
    switch(nuc_id) {
    case 1:
        ret = 'A';
        break;
    case 2:
        ret = 'C';
        break;
    case 3:
        ret = 'M';
        break;
    case 4:
        ret = 'G';
        break;
    case 5:
        ret = 'R';
        break;
    case 6:
        ret = 'S';
        break;
    case 7:
        ret = 'V';
        break;
    case 8:
        ret = 'T';
        break;
    case 9:
        ret = 'W';
        break;
    case 10:
        ret = 'Y';
        break;
    case 11:
        ret = 'H';
        break;
    case 12:
        ret = 'K';
        break;
    case 13:
        ret = 'D';
        break;
    case 14:
        ret = 'B';
        break;
    default:
        ret = 'N';
        break;
    }
    return ret;
}

// A:0, C:1, G:2, T:3
int8_t get_nt (int8_t nuc_id) {
    int8_t ret = 0;
    switch(nuc_id) {
    case 1:
        ret = 0;
        break;
    case 2:
        ret = 1;
        break;
    case 4:
        ret = 2;
        break;
    case 8:
        ret = 3;
        break;
    default:
        ret = -1;
        break;
    }
    return ret;
}

std::vector<int8_t> get_nuc_vec (char c) {
    switch (c) {
    case 'a':
    case 'A':
        return std::vector<int8_t> {0};
    case 'c':
    case 'C':
        return std::vector<int8_t> {1};
    case 'g':
    case 'G':
        return std::vector<int8_t> {2};
    case 't':
    case 'T':
        return std::vector<int8_t> {3};
    case 'R':
        return std::vector<int8_t> {0,2};
    case 'Y':
        return std::vector<int8_t> {1,3};
    case 'S':
        return std::vector<int8_t> {1,2};
    case 'W':
        return std::vector<int8_t> {0,3};
    case 'K':
        return std::vector<int8_t> {2,3};
    case 'M':
        return std::vector<int8_t> {0,1};
    case 'B':
        return std::vector<int8_t> {1,2,3};
    case 'D':
        return std::vector<int8_t> {0,2,3};
    case 'H':
        return std::vector<int8_t> {0,1,3};
    case 'V':
        return std::vector<int8_t> {0,1,2};
    case 'n':
    case 'N':
        return std::vector<int8_t> {0,1,2,3};
    default:
        return std::vector<int8_t> {0,1,2,3};
    }
}
std::vector<int8_t> get_nuc_vec_from_id (int8_t nuc_id) {
    return get_nuc_vec(get_nuc(nuc_id));
}

void getNodeDFS(Parsimony::data &data, panmanUtils::Node* node, 
    std::vector<std::vector<std::pair<int, std::vector<int>>>> &globalCoords_t, 
    std::vector<std::vector<std::pair<char, std::vector<char>>>> &pseudoRoot, 
    std::vector<std::vector<std::pair<char, std::vector<char>>>> &sequence,
    std::vector< bool >  &blockExists,
    std::vector< bool >  &blockStrand){
    // write nuc mutations 
    auto mutation_list = data.add_node_mutations();

    std::vector< std::tuple< int32_t, bool, bool, bool, bool > > blockMutationInfo;

    // Block Mutations
    for(auto mutation: node->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;
        if(type == 1) {
            // insertion
            bool oldStrand;
            bool oldMut;
            oldStrand = blockStrand[primaryBlockId];
            oldMut = blockExists[primaryBlockId];
            blockExists[primaryBlockId] = true;
            // if insertion of inverted block takes place, the strand is backwards
            blockStrand[primaryBlockId] = !inversion;
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, true, !inversion) );
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
                blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, oldMut, !oldStrand) );
            } else {
                // Actually a deletion
                oldStrand = blockStrand[primaryBlockId];
                oldMut = blockExists[primaryBlockId];
                blockExists[primaryBlockId] = false;

                // resetting strand to true during deletion
                blockStrand[primaryBlockId] = true;
            }
            blockMutationInfo.push_back( std::make_tuple(mutation.primaryBlockId, oldMut, oldStrand, false, true) );

        }   
    }

    // For backtracking. primaryBlockId, secondaryBlockId, pos, gapPos, (oldVal, newVal) in substitution, ('-', newVal) in insertion, (oldVal, '-') in deletion
    std::vector< std::tuple< int32_t, int, int, char, char > > mutationInfo;

    for (int i=0; i<node->nucMutation.size(); i++) {
        int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
        int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;
        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = node->nucMutation[i].type();
        char newVal = '-';

        if(type < 3) {
            // Either S, I or D
            int len = node->nucMutation[i].length();

            if(primaryBlockId >= sequence.size()) {
                std::cout << primaryBlockId << " " << sequence.size() << std::endl;
            }

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                        newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getNucCode(j));
                        sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition+j]);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition+j]));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                        newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getNucCode(j));
                        sequence[primaryBlockId][nucPosition+j].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition+j].first);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition+j].first));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                        newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getNucCode(j));
                        sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, newVal));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition+j]);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition+j]));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                        newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getNucCode(j));
                        sequence[primaryBlockId][nucPosition+j].first = newVal;
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition+j].first);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition+j].first));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
                // Deletion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
                        sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j, oldVal, '-'));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition+j]);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition+j]));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        auto mut = mutation_list->add_mutation();
                        char oldVal = sequence[primaryBlockId][nucPosition+j].first;
                        sequence[primaryBlockId][nucPosition+j].first = '-';
                        mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
                        
                        mut->set_position(globalCoords_t[primaryBlockId][nucPosition+j].first);
                        mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                        mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition+j].first));
                        for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                            mut->add_mut_nuc(nuc);
                        }
                    }
                }
            }
        } else {
            if(type == panmanUtils::NucMutationType::NSNPS) {
                // SNP Substitution
                newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getFirstNucCode());
                if(nucGapPosition != -1) {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                    sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition]);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition]));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                } else {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].first;
                    sequence[primaryBlockId][nucPosition].first = newVal;
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].first);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].first));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NSNPI) {
                // SNP Insertion
                newVal = panmanUtils::getNucleotideFromCode(node->nucMutation[i].getFirstNucCode());
                if(nucGapPosition != -1) {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                    sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition]);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition]));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                } else {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].first;
                    sequence[primaryBlockId][nucPosition].first = newVal;
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].first);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].first));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NSNPD) {
                // SNP Deletion
                if(nucGapPosition != -1) {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
                    sequence[primaryBlockId][nucPosition].second[nucGapPosition] = '-';
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].second[nucGapPosition]);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].second[nucGapPosition]));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                } else {
                    auto mut = mutation_list->add_mutation();
                    char oldVal = sequence[primaryBlockId][nucPosition].first;
                    sequence[primaryBlockId][nucPosition].first = '-';
                    mutationInfo.push_back(std::make_tuple(primaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
                    
                    mut->set_position(globalCoords_t[primaryBlockId][nucPosition].first);
                    mut->set_par_nuc(panmanUtils::getCodeFromNucleotide(oldVal));
                    mut->set_ref_nuc(panmanUtils::getCodeFromNucleotide(pseudoRoot[primaryBlockId][nucPosition].first));
                    for (auto nuc: get_nuc_vec_from_id(panmanUtils::getCodeFromNucleotide(newVal))) {
                        mut->add_mut_nuc(nuc);
                    }
                }
            }
        }
    }

    for(panmanUtils::Node* child: node->children) {
        getNodeDFS(data, child, globalCoords_t, pseudoRoot, sequence, blockExists, blockStrand);
    }

    // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        blockExists[std::get<0>(mutation)] = std::get<1>(mutation);
        blockStrand[std::get<0>(mutation)] = std::get<2>(mutation);

    }

    // Undo nuc mutations when current node and its subtree have been processed
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<2>(mutation) != -1) {
            sequence[std::get<0>(mutation)][std::get<1>(mutation)].second[std::get<2>(mutation)] = std::get<3>(mutation);
        } else {
            sequence[std::get<0>(mutation)][std::get<1>(mutation)].first = std::get<3>(mutation);
        }

    }
}

void panmanUtils::panmanToUsher(panmanUtils::Tree* panmanTree, std::string refName, std::string filename,std::string refSeq) {
    std::vector<std::vector<std::pair<int, std::vector<int>>>> globalCoords_t;
    std::vector<std::vector<std::pair<char, std::vector<char>>>> pseudoRoot;
    
    getCoordMap(panmanTree, globalCoords_t);
    getPseudoRoot(panmanTree, pseudoRoot);

    std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence = pseudoRoot;

    std::vector< bool >  blockExists(panmanTree->blocks.size() + 1, false, {});
    std::vector< bool >  blockStrand(panmanTree->blocks.size() + 1, true, {});

    panmanUtils::Node* root = panmanTree->root;
    
    // Write Usher
    Parsimony::data data;
    data.set_newick(panmanTree->getNewickString(root));

    getNodeDFS(data, root, globalCoords_t, pseudoRoot, sequence, blockExists, blockStrand);

    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outbuf;

    if (filename.find(".gz\0") != std::string::npos) {
        try {
            outbuf.push(boost::iostreams::gzip_compressor());
            outbuf.push(outfile);
            std::ostream outstream(&outbuf);
            data.SerializeToOstream(&outstream);
            boost::iostreams::close(outbuf);
            outfile.close();
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        data.SerializeToOstream(&outfile);
        outfile.close();
    }

    return;
}