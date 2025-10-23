#include "panmanUtils.hpp"

#define epsilon 0.0001

struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        auto [a, b, c] = t;
        std::size_t h1 = std::hash<int>{}(a);
        std::size_t h2 = std::hash<int>{}(b);
        std::size_t h3 = std::hash<int>{}(c);
        // Combine the hashes
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct Hash {
    std::size_t operator()(const std::tuple<int, int, std::string>& t) const {
        auto [a, b, c] = t;
        std::size_t h1 = std::hash<int>{}(a);
        std::size_t h2 = std::hash<int>{}(b);
        std::size_t h3 = std::hash<std::string>{}(c);
        // Combine the hashes
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct Mut {
    int pos;
    int gapPos;
    int8_t type;
    int len;
    std::string chars;

    Mut(int pos, int gapPos, int8_t type, int len, char chars) {
        this->pos=pos;
        this->gapPos=gapPos;
        this->type=type;
        this->len=len;
        this->chars=chars;
    }

    Mut(int pos, int gapPos, int8_t type, int len, std::string chars) {
        this->pos=pos;
        this->gapPos=gapPos;
        this->type=type;
        this->len=len;
        this->chars=chars;
    }
};

void dfs(panmanUtils::Node* currNode, int& afrC, int& nonAfrC) {
    if (currNode->children.size() == 0) {
        if (currNode->annotations.size() > 0){
            if (currNode->annotations[0] == "African"){
                afrC++;
            } else if (currNode->annotations[0] == "Non-African"){
                nonAfrC++;
            }
        }
        return;
    }
    for (auto child: currNode->children){
        dfs(child, afrC, nonAfrC);
    }
}

void modifyCoordinate(  panmanUtils::Node* node, 
                        std::unordered_map<panmanUtils::Node*, std::vector<Mut>>& allMuts, 
                        std::unordered_map<int, std::pair<int, int>>& PanMATToRefCoordinateMap) {

    
    std::vector<Mut> nodeMuts;
    for(size_t i = 0; i < node->nucMutation.size(); i++) {
        int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;

        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            // Either S, I or D
            int len = ((node->nucMutation[i].mutInfo) >> 4);

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                for(int j = 0; j < len; j++) {
                    newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                    auto pos = PanMATToRefCoordinateMap[nucPosition+j];
                    nodeMuts.emplace_back(pos.first, pos.second, type, 1, newVal);
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                for(int j = 0; j < len; j++) {
                    newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                    auto pos = PanMATToRefCoordinateMap[nucPosition+j];
                    nodeMuts.emplace_back(pos.first, pos.second, type, 1, newVal);
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
                // Deletion
                for(int j = 0; j < len; j++) {
                    auto pos = PanMATToRefCoordinateMap[nucPosition+j];
                    nodeMuts.emplace_back(pos.first, pos.second, type, 1, '-');
                }
            }
        } else {
            std::cout << "unexpected\n";
        }
    }

    /* join mutations if possible */
    std::sort(nodeMuts.begin(), nodeMuts.end(), [](const Mut& a, const Mut& b) {
        if(a.pos == b.pos) {
            return a.gapPos < b.gapPos;
        }
        return a.pos < b.pos;
    });

    if (nodeMuts.size() == 0){
        return;
    }

    std::vector<Mut> joinedMuts;
    Mut& baseMut = nodeMuts[0];
    int len=baseMut.len;
    std::string chars=baseMut.chars;
    for(size_t i = 1; i < nodeMuts.size(); i++) {
        if( (nodeMuts[i].type != nodeMuts[i-1].type) || 
            ((nodeMuts[i].pos != nodeMuts[i-1].pos + nodeMuts[i-1].len) &&
            ((nodeMuts[i].pos == nodeMuts[i-1].pos) && (nodeMuts[i].gapPos != nodeMuts[i-1].gapPos + nodeMuts[i-1].len)))
         ){
            joinedMuts.emplace_back(baseMut.pos, baseMut.gapPos, baseMut.type, len, chars);
            baseMut = nodeMuts[i];
            len = baseMut.len;
            chars = baseMut.chars;
        } else {
            len++;
            chars += nodeMuts[i].chars;
        } 
    }
    joinedMuts.emplace_back(baseMut.pos, baseMut.gapPos, baseMut.type, len, chars);

    allMuts[node]=joinedMuts;

}

void panmanUtils::Tree::ratioTest(std::ostream& fout) {

    auto node = this->root;

    // Get block sequnece of the Tip
    std::vector< bool >  blockSequence(blocks.size() + 1, false, {});

    // Blocks length
    std::unordered_map<int, int> blockLengths;

    
    for(auto mutation: node->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;
        if(type == 1) {
            // insertion
            blockSequence[primaryBlockId] = true;
        } else {
            // deletion
            if(!inversion) {
                blockSequence[primaryBlockId] = false;
            }
        }
    }

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
            }
        } else {
            int len=0;
            for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                len += gaps[i].nucGapLength[j];
            }
            blockLengths[primaryBId] += len;
        }
    }

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

    /* PanMAT to Reference Coordinate */
    std::unordered_map<int, std::pair<int, int>> PanMATToRefCoordinateMap;

    bool aligned = true;

    int nucCo=0, nucGapCo=-1;
    for(size_t i = 0; i < blockExists.size(); i++) {
        // Non-gap block - the only type being used currently
        if(blockExists[i]) {
            if(blockStrand[i]) {
                for(size_t j = 0; j < sequence[i].size(); j++) {
                    if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        PanMATToRefCoordinateMap[j]=std::make_pair(nucCo++,-1);
                        nucGapCo=-1;
                    } else if(aligned) {
                        PanMATToRefCoordinateMap[j]=std::make_pair(nucCo,++nucGapCo);
                    }
                }
            } // Reverse strand not handled
        } 
    }

    
    /* Convert Mutations to New Coordinates */
    std::unordered_map<Node*, std::vector<Mut>> allMuts;

    int c=0;
    for (auto &node: allNodes){
        if (node.second != this->root) {
            modifyCoordinate(node.second, allMuts, PanMATToRefCoordinateMap);
        }
    }

    /* Print Annotations
    for (auto &node: allNodes){
        if (node.second != this->root) {
            std::cout << node.first << ": ";
            for (auto i: node.second->annotations){
                std::cout << i << " --- ";
            }
            std::cout << std::endl;
        }
    }
    */

    // /* List of all mutations (mutation inde) */
    // std::vector<std::tuple<int, int8_t, int, std::string, std::string>> mutations; /*pos, type, len, chars, node*/
    // for (auto &node: allNodes){
    //     if (node.second != this->root) {
    //         for (auto &mut: allMuts[node.second]){
    //             mutations.emplace_back(mut.pos+mut.gapPos, mut.type, mut.len, mut.chars, node.first);
    //         }
    //     }
    // }

    int totalAfr=0;
    int totalNonAfr=0;
    for (auto node:this->allNodes){
        if (node.second->children.size() == 0) {
            if (node.second->annotations.size() > 0){
                if (node.second->annotations[0] == "African"){
                    totalAfr++;
                } else if (node.second->annotations[0] == "Non-African"){
                    totalNonAfr++;
                }
            }
        }
    }

    std::cout << "Total African Tips: " << totalAfr << std::endl;
    std::cout << "Total Non-African Tips: " << totalNonAfr << std::endl;

    std::unordered_map<std::string, std::pair<int, int>> nodeAncMap;
    for (auto node:this->allNodes){
        if (node.second->identifier != this->root->identifier) {
            int afrC=0;
            int nonAfrC=0;
            panmanUtils::Node* currNode = node.second;
            dfs(currNode, afrC, nonAfrC);
            nodeAncMap[node.first]=std::make_pair(afrC, nonAfrC);
        }
    }

    /* List of all mutations (mutation inde) */
    std::vector<std::tuple<int, int8_t, int, std::string, double, double, double>> mutationsFinal; /*pos, type, len, chars, node*/
    for (auto &node: allNodes){
        if (node.second != this->root) {
            for (auto &mut: allMuts[node.second]){
                double fracAfr = (double)nodeAncMap[node.first].first / (double)totalAfr;
                double fracNonAfr = (double)nodeAncMap[node.first].second /(double) totalNonAfr;
                double ratioTest = fracAfr/(fracNonAfr+epsilon);
                mutationsFinal.emplace_back(mut.pos+mut.gapPos, mut.type, mut.len, mut.chars, 
                                        fracAfr, 
                                        fracNonAfr,
                                        ratioTest);
            }
        }
    }

    std::sort(mutationsFinal.begin(), mutationsFinal.end(), [](const std::tuple<int, int8_t, int, std::string, int, int, int>& a,
                                                     const std::tuple<int, int8_t, int, std::string, int, int, int>& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    for (auto &mut: mutationsFinal){
        fout << std::get<0>(mut) << "\t" 
             << (int)std::get<1>(mut) << "\t"
             << std::get<2>(mut) << "\t"
             << std::get<4>(mut) << "\t"
             << std::get<5>(mut) << "\t"
             << std::get<6>(mut) << std::endl;
    }

    std::cout << "Total Mutations: " << mutationsFinal.size() << std::endl;
    

    return;

}
