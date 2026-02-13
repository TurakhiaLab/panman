#include "panmanUtils.hpp"

#define epsilon 0.0001

#include <iomanip>
#include <cmath>

// using the error function (since chi2_1 ~ Z^2)
double chi_square_p_value(double chi2_stat) {
    double z = std::sqrt(chi2_stat);
    // CDF for chi-square(1) = 2 * Phi(z) - 1, where Phi is the normal CDF
    double p = std::erfc(z / std::sqrt(2));  // survival function 1 - Phi(z)
    return p; // two-tailed p-value
}

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
    std::size_t operator()(const std::pair<int, int8_t>& t) const {
        auto [a, b] = t;
        std::size_t h1 = std::hash<int>{}(a);
        std::size_t h2 = std::hash<int>{}(b);
        // Combine the hashes
        return h1 ^ (h2 << 1);
    }
};

struct HashTupleToInt {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        auto [a, b, c] = t;
        std::size_t h1 = std::hash<int>{}(a);
        std::size_t h2 = std::hash<int>{}(b);
        std::size_t h3 = std::hash<int>{}(c);
        // Combine the hashes
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct HashPairToBool {
    std::size_t operator()(const std::pair<int, char>& t) const {
        auto [a, b] = t;
        std::size_t h1 = std::hash<int>{}(a);
        std::size_t h2 = std::hash<int>{}(b);
        // Combine the hashes
        return h1 ^ (h2 << 1);
    }
};

struct Mut {
    int pos;
    int gapPos;
    int8_t type;
    int len;
    std::string chars;

    // Default constructor needed for containers (e.g. unordered_map::operator[])
    Mut() : pos(0), gapPos(0), type(0), len(0), chars() {}

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

struct Var {
    int pos;
    int8_t type;
    char chars;

    // Default constructor so maps/vectors can default-construct
    Var() : pos(0), type(0), chars('-') {}

    Var(int pos, int8_t type, char chars) {
        this->pos=pos;
        this->type=type;
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

void getVar( panmanUtils::Node* node,
             const std::vector<panmanUtils::Node*>& path,
             std::unordered_map<int, Var>& varMap) {

    for (auto currNode: path){
        for(size_t i = 0; i < currNode->nucMutation.size(); i++) {
            int32_t primaryBlockId = currNode->nucMutation[i].primaryBlockId;

            int32_t nucPosition = currNode->nucMutation[i].nucPosition;
            int32_t nucGapPosition = currNode->nucMutation[i].nucGapPosition;
            uint32_t type = (currNode->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3) {
                // Either S, I or D
                int len = ((currNode->nucMutation[i].mutInfo) >> 4);

                if(type == panmanUtils::NucMutationType::NS) {
                    // Substitution
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((currNode->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        varMap[nucPosition+j]=Var(nucPosition+j, type, newVal);
                    }
                } else if(type == panmanUtils::NucMutationType::NI) {
                    // Insertion
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((currNode->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        varMap[nucPosition+j]=Var(nucPosition+j, type, newVal);
                    }
                } else if(type == panmanUtils::NucMutationType::ND) {
                    // Deletion
                    for(int j = 0; j < len; j++) {
                        varMap[nucPosition+j]=Var(nucPosition+j, type, '-');
                    }
                }
            } else {
                std::cout << "unexpected\n";
            }
        }
    }
}

void getRefSeq(panmanUtils::Node* node, 
    std::vector<std::pair<int, std::vector<std::pair<char, std::vector<char>>>>>& refSeq,
    std::vector< bool >  &blockExists,
    std::vector< bool >  &blockStrand
){
    for(size_t i = 0; i < node->nucMutation.size(); i++) {
        int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
        if (!blockExists[primaryBlockId]) continue;

        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            int len = ((node->nucMutation[i].mutInfo) >> 4);

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        refSeq[primaryBlockId].second[nucPosition].second[j] = newVal;
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        refSeq[primaryBlockId].second[nucPosition+j].first = newVal;
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        refSeq[primaryBlockId].second[nucPosition].second[j] = newVal;
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        refSeq[primaryBlockId].second[nucPosition+j].first = newVal;
                    }
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
                // Deletion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        refSeq[primaryBlockId].second[nucPosition].second[j] = '-';
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        refSeq[primaryBlockId].second[nucPosition+j].first = '-';  
                    }
                }
            }
        } else {
            if(type == panmanUtils::NucMutationType::NSNPS || type == panmanUtils::NucMutationType::NSNPI) {
                // SNP Substitution
                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                refSeq[primaryBlockId].second[nucPosition].first = newVal;    
            } else if(type == panmanUtils::NucMutationType::NSNPD) {
                // SNP Deletion
                refSeq[primaryBlockId].second[nucPosition].first = '-';   
            }
        }
    }
}


void getVariations(panmanUtils::Node* node, 
    std::unordered_map<std::tuple<int,int, int>, int, HashTupleToInt> &mapFromPanmanToLinearCoordinate, 
    std::vector< bool >  &blockExists,
    std::vector< bool >  &blockStrand,
    std::map<int, char>& variations
){


    for(size_t i = 0; i < node->nucMutation.size(); i++) {
        int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
        if (!blockExists[primaryBlockId]) continue;

        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';

        if(type < 3) {
            int len = ((node->nucMutation[i].mutInfo) >> 4);

            if(type == panmanUtils::NucMutationType::NS) {
                // Substitution
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)]]=newVal;
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, nucGapPosition)]]=newVal;
                    }
                }
            } else if(type == panmanUtils::NucMutationType::NI) {
                // Insertion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)]]=newVal;
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, nucGapPosition)]]=newVal;

                    }
                }
            } else if(type == panmanUtils::NucMutationType::ND) {
                // Deletion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition+j)]]='-';    
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition+j, nucGapPosition)]]='-';    
                    }
                }
            }
        } else {
            if(type == panmanUtils::NucMutationType::NSNPS || type == panmanUtils::NucMutationType::NSNPI) {
                // SNP Substitution
                newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)]]=newVal;    
            } else if(type == panmanUtils::NucMutationType::NSNPD) {
                // SNP Deletion
                variations[mapFromPanmanToLinearCoordinate[std::make_tuple(primaryBlockId, nucPosition, nucGapPosition)]]='-';
            }
        }
    }
}

int presentInAfr(std::tuple<int, int8_t, char, std::string>& var, 
                std::unordered_map<panmanUtils::Node*, std::unordered_map<int, Var>>& allVars,
                std::string ancestry="") {
    int O = 0;
    for (auto &entry: allVars){
        if (ancestry != "" && entry.first->annotations.size() > 0 && entry.first->annotations[0] != ancestry){
            continue;
        }
        if (entry.second.find(std::get<0>(var)) != entry.second.end()){
            Var v = entry.second[std::get<0>(var)];
            if (v.type == std::get<1>(var) && v.chars == std::get<2>(var)){
                O++;
            }
        }
        
    }
    return O;
}

void applyBlockMutations(panmanUtils::Node* node, 
                        std::vector< bool >&  blockExists,
                        std::vector< bool >&  blockStrand
){
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
        } else {
            bool oldMut;
            bool oldStrand;
            if(inversion) {
                // This means that this is not a deletion, but instead an inversion
                oldStrand = blockStrand[primaryBlockId];
                oldMut = blockExists[primaryBlockId];
                blockStrand[primaryBlockId] = !oldStrand;
                
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

void panmanUtils::Tree::ratioTest(std::ostream& fout) {

    auto node = this->root;
    std::vector<std::pair<int, std::vector<std::pair<char, std::vector<char>>>>> refSeq(this->blocks.size() + 1);

    std::vector<std::pair<int, std::vector<std::pair<int, int>>>> sequence(this->blocks.size() + 1);
    for (size_t block_id = 0; block_id < blocks.size(); block_id++) {
        int32_t blockId = ((int32_t)blocks[block_id].primaryBlockId);
        sequence[block_id].first=blockId;
        refSeq[block_id].first=blockId;
        int nucPos=0;
        for (size_t nuc_pos = 0; nuc_pos < blocks[block_id].consensusSeq.size(); nuc_pos++) {
            bool endFlag = false;
            for (size_t k = 0; k < 8; k++) {
                const int nucCode = (((blocks[block_id].consensusSeq[nuc_pos]) >> (4 * (7 - k))) & 15);
                if (nucCode == 0) {
                    endFlag = true;
                    break;
                }
                sequence[block_id].second.push_back(std::make_pair(nucPos++, -1));
                refSeq[block_id].second.push_back(std::make_pair(getNucleotideFromCode(nucCode), std::vector<char>()));
            }
            if (endFlag){
                break;
            }
            
        }
    }

    for (const auto &g: this->gaps) {
        auto blockId=g.primaryBlockId;
        for (int i=0;i<g.nucPosition.size();i++){
            sequence[blockId].second[g.nucPosition[i]].second = g.nucGapLength[i];
            refSeq[blockId].second[g.nucPosition[i]].second.resize(g.nucGapLength[i]);
            assert(sequence[blockId].second[g.nucPosition[i]].first == g.nucPosition[i]);
        }
    }

    /* Convert PanMAN to linear (msa) coordinate */
    std::unordered_map<std::tuple<int,int, int>, int, HashTupleToInt> mapFromPanmanToLinearCoordinate;

    const auto &blocks = this->blocks;
    const auto &gaps = this->gaps;

    int index=1; // VCF index starts from 1

    for (size_t block_id = 0; block_id < sequence.size(); block_id++) {
        auto blockID=sequence[block_id].first;
        for (size_t nuc_pos = 0; nuc_pos < sequence[block_id].second.size(); nuc_pos++) {
            auto nucPos=sequence[block_id].second[nuc_pos].first;
            if (sequence[block_id].second[nuc_pos].second == -1)
                mapFromPanmanToLinearCoordinate[std::make_tuple(blockID, nucPos, -1)]=index++;
            else{
                for (size_t gap_pos=0; gap_pos < sequence[block_id].second[nuc_pos].second; gap_pos++){
                    mapFromPanmanToLinearCoordinate[std::make_tuple(blockID, nucPos, gap_pos)]=index++;
                }
            }
        }
    }
    
    // Maintain block Sequence to decide if mutation at a node needs to be considered
    std::vector< bool >  blockExists(blocks.size() + 1, false, {});
    std::vector< bool >  blockStrand(blocks.size() + 1, true, {});
    applyBlockMutations(this->root, blockExists, blockStrand);

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
 
    /* Working on varitations */
    std::unordered_map<Node*, std::map<int, char>> allVars;
    
    for (const auto &node: allNodes){
        if (node.second->children.size() != 0) continue;
        std::vector<Node*> pathToRoot;
        // copy of blockExists and blockStrand, since we are not reversing block mutation
        auto copyBlockExists=blockExists;
        auto copyBlockStrand=blockStrand;

        auto tmpNode=node.second;
        while(tmpNode!=this->root){
            pathToRoot.push_back(tmpNode);
            tmpNode=tmpNode->parent;
        }

        std::reverse(pathToRoot.begin(), pathToRoot.end());

        std::map<int, char> localVars;
        for (const auto& nodeInPath: pathToRoot) {
            applyBlockMutations(nodeInPath, copyBlockExists, copyBlockStrand);
            getVariations(nodeInPath, mapFromPanmanToLinearCoordinate, blockExists, blockStrand, localVars);
        }
        allVars[node.second]=localVars;
    }

    /* Extract reference MSA to find MSA to Ref Coordinate map */
    std::unordered_map<int, int> mapFromMSACoordinateToRefCoordinate;
    for (const auto &node: allNodes){
        if (node.second->identifier != "GRCh38#0#chr1") continue;
        std::cout << "FOund reference node\n";
        std::vector<Node*> pathToRoot;
        // copy of blockExists and blockStrand, since we are not reversing block mutation
        auto copyBlockExists=blockExists;
        auto copyBlockStrand=blockStrand;

        auto tmpNode=node.second;
        while(tmpNode!=this->root){
            pathToRoot.push_back(tmpNode);
            tmpNode=tmpNode->parent;
        }
        pathToRoot.push_back(tmpNode); // add root node too here as we are trying to extract tip node sequence
        std::reverse(pathToRoot.begin(), pathToRoot.end());

        for (const auto& nodeInPath: pathToRoot) {
            applyBlockMutations(nodeInPath, copyBlockExists, copyBlockStrand);
            getRefSeq(nodeInPath, refSeq, blockExists, blockStrand);
        }
    }

    int refPos=1;
    std::string refSequence="";
    for(int bId=0; bId<refSeq.size(); bId++) {
        for(int nPos=0; nPos<refSeq[bId].second.size(); nPos++) {
            // reference base
            char refBase = refSeq[bId].second[nPos].first;
            // map to linear coordinate
            int msaPos = mapFromPanmanToLinearCoordinate[std::make_tuple(refSeq[bId].first, nPos, -1)]; 
            mapFromMSACoordinateToRefCoordinate[msaPos] = refPos; // since ref is same as msa here
            if (refBase != '-') {
                refPos++;
            }
            refSequence.push_back(refBase);

            // inserted bases
            for(int gPos=0; gPos<refSeq[bId].second[nPos].second.size(); gPos++) {
                char insBase = refSeq[bId].second[nPos].second[gPos];
                int msaInsPos = mapFromPanmanToLinearCoordinate[std::make_tuple(refSeq[bId].first, nPos, gPos)];
                mapFromMSACoordinateToRefCoordinate[msaInsPos] = refPos; // since ref is same as msa here
                if (insBase != '-') {
                    refPos++;
                }
                refSequence.push_back(insBase);
            }
        }
    }

    



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
    
    std::unordered_map<std::pair<int, char>, bool, HashPairToBool> visitedVar;
    for (auto &nodeVar: allVars){
        for (auto &var: nodeVar.second) {
            int pos=var.first;
            char c=var.second;
            if (visitedVar.find(std::make_pair(var.first, var.second))!=visitedVar.end()) continue;
            visitedVar[std::make_pair(var.first, var.second)]=1;
            double O=0.,O_dash=0.,E=0.;

            for (auto &itrNodeVar: allVars){
                if (itrNodeVar.second.find(pos) != itrNodeVar.second.end() && itrNodeVar.second[pos]==c) {
                    E++;
                    if (itrNodeVar.first->annotations[0] == "African") O++;
                }
            }
            O_dash=E-O;
            if (O == 0 || E == 0){
                continue;
            }

            E=E*(double)totalAfr/((double)(totalAfr + totalNonAfr));

            double chi2 = std::pow(O - E, 2) / E;
            double p_value = chi_square_p_value(chi2);
            double frac_A = (double)O / (double)totalAfr;
            double frac_nonA = (double)O_dash / (double)totalNonAfr;
            double ratio = frac_nonA > 0 ? frac_A / frac_nonA : INFINITY;

            fout << pos << "\t" 
                << mapFromMSACoordinateToRefCoordinate[pos] << "\t"
                << refSequence[mapFromMSACoordinateToRefCoordinate[pos]-1] << "\t"
                << c << "\t"
                << O << "\t"
                << O_dash << "\t"
                << E << "\t"
                << std::fixed << std::setprecision(4) << chi2 << "\t"
                << std::scientific << std::setprecision(4) << p_value << "\t"
                << std::fixed <<  std::setprecision(4) << frac_A << "\t"
                << std::fixed <<  std::setprecision(4) << frac_nonA << "\t"
                << std::fixed <<  std::setprecision(4) << ratio << std::endl;
            

        }
        
    }

    return;

}

// void panmanUtils::Tree::ratioTest(std::ostream& fout) {

//     auto node = this->root;

//     // Get block sequnece of the Tip
//     std::vector< bool >  blockSequence(blocks.size() + 1, false, {});

//     // Blocks length
//     std::unordered_map<int, int> blockLengths;

    
//     for(auto mutation: node->blockMutation) {
//         int32_t primaryBlockId = mutation.primaryBlockId;
//         bool type = mutation.blockMutInfo;
//         bool inversion = mutation.inversion;
//         if(type == 1) {
//             // insertion
//             blockSequence[primaryBlockId] = true;
//         } else {
//             // deletion
//             if(!inversion) {
//                 blockSequence[primaryBlockId] = false;
//             }
//         }
//     }

//     // Expanding blocks only if exist in tip 
//     std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence(blocks.size() + 1);
//     std::vector< bool >  blockExists(blocks.size() + 1, false, {});
//     std::vector< bool >  blockStrand(blocks.size() + 1, true, {});


//     int32_t maxBlockId = 0;

//     // Create consensus sequence of blocks
//     for(size_t i = 0; i < blocks.size(); i++) {
//         int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
//         blockLengths[primaryBlockId] = 0;
//         maxBlockId = std::max(maxBlockId, primaryBlockId);
//         if (blockSequence[primaryBlockId]) {
//             for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
//                 bool endFlag = false;
//                 for(size_t k = 0; k < 8; k++) {
//                     const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

//                     if(nucCode == 0) {
//                         endFlag = true;
//                         break;
//                     }
//                     // len++;
//                     const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
//                     sequence[primaryBlockId].push_back({nucleotide, {}});
//                 }

//                 if(endFlag) {
//                     break;
//                 }
//             }
//             // End character to incorporate for gaps at the end
//             sequence[primaryBlockId].push_back({'x', {}});
//             // blockLengths[primaryBlockId] += len;
//         } else {
//             int len = 0;
//             for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
//                 bool endFlag = false;
//                 for(size_t k = 0; k < 8; k++) {
//                     const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
//                     if(nucCode == 0) {
//                         endFlag = true;
//                         break;
//                     }
//                     len++;
//                     const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
//                 }

//                 if(endFlag) {
//                     break;
//                 }
//             }
//             blockLengths[primaryBlockId] += len;
//         }
//     }

//     sequence.resize(maxBlockId + 1);
//     blockExists.resize(maxBlockId + 1);
//     blockStrand.resize(maxBlockId + 1);

//     // Assigning nucleotide gaps in blocks
//     for(size_t i = 0; i < gaps.size(); i++) {
//         int32_t primaryBId = (gaps[i].primaryBlockId);
//         int32_t secondaryBId = (gaps[i].secondaryBlockId);
//         if (blockSequence[primaryBId]){
//             for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
//                 int len = gaps[i].nucGapLength[j];
//                 int pos = gaps[i].nucPosition[j];
//                 sequence[primaryBId][pos].second.resize(len, '-');
//             }
//         } else {
//             int len=0;
//             for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
//                 len += gaps[i].nucGapLength[j];
//             }
//             blockLengths[primaryBId] += len;
//         }
//     }

//     for(auto mutation: node->blockMutation) {
//         int32_t primaryBlockId = mutation.primaryBlockId;
//         bool type = mutation.blockMutInfo;
//         bool inversion = mutation.inversion;
//         if (blockSequence[primaryBlockId]) {
//             if(type == 1) {
//                 // insertion
//                 bool oldStrand;
//                 bool oldMut;
//                 oldStrand = blockStrand[primaryBlockId];
//                 oldMut = blockExists[primaryBlockId];
//                 blockExists[primaryBlockId] = true;
//                 // if insertion of inverted block takes place, the strand is backwards
//                 blockStrand[primaryBlockId] = !inversion;
//             } else {
//                 bool oldMut;
//                 bool oldStrand;
//                 if(inversion) {
//                     // This means that this is not a deletion, but instead an inversion
//                     oldStrand = blockStrand[primaryBlockId];
//                     oldMut = blockExists[primaryBlockId];
//                     blockStrand[primaryBlockId] = !oldStrand;
//                     if(oldMut != true) {
//                         // std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
//                     }
//                 } else {
//                     // Actually a deletion
//                     oldStrand = blockStrand[primaryBlockId];
//                     oldMut = blockExists[primaryBlockId];
//                     blockExists[primaryBlockId] = false;
//                     // resetting strand to true during deletion
//                     blockStrand[primaryBlockId] = true;
//                 }
//             }
//         }
//     }

//     // Nuc mutations
//     for(size_t i = 0; i < node->nucMutation.size(); i++) {
//         int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
//         int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

//         if (blockSequence[primaryBlockId]) {
//             int32_t nucPosition = node->nucMutation[i].nucPosition;
//             int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
//             uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
//             char newVal = '-';

//             if(type < 3) {
//                 // Either S, I or D
//                 int len = ((node->nucMutation[i].mutInfo) >> 4);

//                 if(primaryBlockId >= sequence.size()) {
//                     std::cout << primaryBlockId << " " << sequence.size() << std::endl;
//                 }

//                 if(type == panmanUtils::NucMutationType::NS) {
//                     // Substitution
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                             newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                             newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId][nucPosition+j].first = newVal;
//                         }
//                     }
//                 } else if(type == panmanUtils::NucMutationType::NI) {
//                     // Insertion
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                             newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                             newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
//                             sequence[primaryBlockId][nucPosition+j].first = newVal;
//                         }
//                     }
//                 } else if(type == panmanUtils::NucMutationType::ND) {
//                     // Deletion
//                     if(nucGapPosition != -1) {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition+j];
//                             sequence[primaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
//                         }
//                     } else {
//                         for(int j = 0; j < len; j++) {
//                             char oldVal = sequence[primaryBlockId][nucPosition+j].first;
//                             sequence[primaryBlockId][nucPosition+j].first = '-';
//                         }
//                     }
//                 }
//             } else {
//                 if(type == panmanUtils::NucMutationType::NSNPS) {
//                     // SNP Substitution
//                     newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                     } else {
//                         char oldVal = sequence[primaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId][nucPosition].first = newVal;
//                     }
//                 } else if(type == panmanUtils::NucMutationType::NSNPI) {
//                     // SNP Insertion
//                     newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition] = newVal;
//                     } else {
//                         char oldVal = sequence[primaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId][nucPosition].first = newVal;
//                     }
//                 } else if(type == panmanUtils::NucMutationType::NSNPD) {
//                     // SNP Deletion
//                     if(nucGapPosition != -1) {
//                         char oldVal = sequence[primaryBlockId][nucPosition].second[nucGapPosition];
//                         sequence[primaryBlockId][nucPosition].second[nucGapPosition] = '-';
//                     } else {
//                         char oldVal = sequence[primaryBlockId][nucPosition].first;
//                         sequence[primaryBlockId][nucPosition].first = '-';
//                     }
//                 }
//             }
//         }
//     }

//     /* PanMAT to Reference Coordinate */
//     std::unordered_map<int, std::pair<int, int>> PanMATToRefCoordinateMap;

//     bool aligned = true;

//     int nucCo=0, nucGapCo=-1;
//     for(size_t i = 0; i < blockExists.size(); i++) {
//         // Non-gap block - the only type being used currently
//         if(blockExists[i]) {
//             if(blockStrand[i]) {
//                 for(size_t j = 0; j < sequence[i].size(); j++) {
//                     if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
//                         PanMATToRefCoordinateMap[j]=std::make_pair(nucCo++,-1);
//                         nucGapCo=-1;
//                     } else if(aligned) {
//                         PanMATToRefCoordinateMap[j]=std::make_pair(nucCo,++nucGapCo);
//                     }
//                 }
//             } // Reverse strand not handled
//         } 
//     }

    
//     /* Working on varitations */
//     std::unordered_map<Node*, std::unordered_map<int, Var>> allVars;
//     std::unordered_map<Node*, std::vector<Node*>> paths;

//     for (auto &node: allNodes){
//         if (node.second->children.size() == 0) {
//             std::vector<Node*> path;
//             panmanUtils::Node* currNode = node.second;
//             while (currNode != nullptr){
//                 path.push_back(currNode);
//                 currNode = currNode->parent;
//             }
//             path.pop_back(); // remove root
//             std::reverse(path.begin(), path.end());
//             paths[node.second]=path;
//         }
//     }

//     for (auto &node: paths){
//         std::cout << node.first->identifier << "\t";
//         getVar(node.first, node.second, allVars[node.first]);
//         if (allVars[node.first].find(308098) != allVars[node.first].end() && node.first->annotations.size() > 0 && node.first->annotations[0] == "African"){
//             std::cout << "1" << std::endl;
//         } else if (node.first->annotations.size() > 0 && node.first->annotations[0] != "African" && node.first->annotations[0] != "Non-African"){
//             std::cout << "2" << std::endl;
//         } else {
//             std::cout << "0" << std::endl;
//         }
//     }


//     /* To compute mutations for all nodes
//     std::unordered_map<Node*, std::vector<Mut>> allMuts;
//     int c=0;
//     for (auto &node: allNodes){
//         if (node.second != this->root) {
//             modifyCoordinate(node.second, allMuts, PanMATToRefCoordinateMap);
//         }
//     }
//     */


//     /* Print Annotations
//     for (auto &node: allNodes){
//         if (node.second != this->root) {
//             std::cout << node.first << ": ";
//             for (auto i: node.second->annotations){
//                 std::cout << i << " --- ";
//             }
//             std::cout << std::endl;
//         }
//     }
//     */

//     /* List of all variations (mutations) */
//     std::vector<std::tuple<int, int8_t, char, std::string>> variations; /*pos, type, chars, node*/
//     std::unordered_map<int, Var> tempVarMap;
//     for (auto &node: paths){
//         for (auto &var: allVars[node.first]){
//             if (tempVarMap.find(var.first) != tempVarMap.end()){
//                 Var existingVar = tempVarMap[var.first];
//                 if (existingVar.type == var.second.type && existingVar.chars == var.second.chars){
//                     continue;
//                 }
//             }
//             tempVarMap[var.first]=var.second;
//             variations.emplace_back(var.second.pos, var.second.type, var.second.chars, node.first->identifier);
//         }
//     }

//     int totalAfr=0;
//     int totalNonAfr=0;
//     for (auto node:this->allNodes){
//         if (node.second->children.size() == 0) {
//             if (node.second->annotations.size() > 0){
//                 if (node.second->annotations[0] == "African"){
//                     totalAfr++;
//                 } else if (node.second->annotations[0] == "Non-African"){
//                     totalNonAfr++;
//                 }
//             }
//         }
//     }

//     std::cout << "Total African Tips: " << totalAfr << std::endl;
//     std::cout << "Total Non-African Tips: " << totalNonAfr << std::endl;

//     std::unordered_map<std::string, std::pair<int, int>> nodeAncMap;
//     for (auto node:this->allNodes){
//         if (node.second->identifier != this->root->identifier) {
//             int afrC=0;
//             int nonAfrC=0;
//             panmanUtils::Node* currNode = node.second;
//             dfs(currNode, afrC, nonAfrC);
//             nodeAncMap[node.first]=std::make_pair(afrC, nonAfrC);
//         }
//     }

//     /* List of all variations (mutation inde) */
//     // std::vector<std::tuple<int, int8_t, char, double, double, double>> mutationsFinal; /*pos, type, chars, node*/
//     // for (auto &node: allNodes){
//     //     if (node.second != this->root) {
//     //         for (auto &mut: allMuts[node.second]){
//     //             double fracAfr = (double)nodeAncMap[node.first].first / (double)totalAfr;
//     //             double fracNonAfr = (double)nodeAncMap[node.first].second /(double) totalNonAfr;
//     //             double ratioTest = fracAfr/(fracNonAfr+epsilon);
//     //             mutationsFinal.emplace_back(mut.pos+mut.gapPos, mut.type, mut.len, mut.chars, 
//     //                                     fracAfr, 
//     //                                     fracNonAfr,
//     //                                     ratioTest);
//     //         }
//     //     }
//     // }

//     // std::sort(mutationsFinal.begin(), mutationsFinal.end(), [](const std::tuple<int, int8_t, int, std::string, int, int, int>& a,
//     //                                                  const std::tuple<int, int8_t, int, std::string, int, int, int>& b) {
//     //     return std::get<0>(a) < std::get<0>(b);
//     // });

//     // // for (auto &mut: mutationsFinal){
//     // //     fout << std::get<0>(mut) << "\t" 
//     // //          << (int)std::get<1>(mut) << "\t"
//     // //          << std::get<2>(mut) << "\t"
//     // //          << std::get<4>(mut) << "\t"
//     // //          << std::get<5>(mut) << "\t"
//     // //          << std::get<6>(mut) << std::endl;
//     // // }

//     // std::cout << "Total Mutations: " << mutationsFinal.size() << std::endl;
    
//     // std::cout << "\nChi-square statistic\t" 
//     //           << "p-value (approx, df=1)\t"
//     //           << "African fraction\t"
//     //           << "Non-African fraction\t"
//     //           << "Ratio (frac_A / frac_nonA)\n";

//     return;
//     std::sort(variations.begin(), variations.end(), [](const std::tuple<int, int8_t, char, std::string>& a,
//                                                      const std::tuple<int, int8_t, char, std::string>& b) {
//         return std::get<0>(a) < std::get<0>(b);
//     });
    
//     std::unordered_map<std::pair<int, int8_t>, std::tuple<int, int8_t, std::string, int, int, int, double, double, double, double, double>, Hash> finalVars;
//     std::cout << variations.size() << " variations found.\n";
//     for (auto &var: variations){
//         int O = presentInAfr(var, allVars, "African");
//         int O_dash = presentInAfr(var, allVars, "Non-African");
//         int E = (double)(presentInAfr(var, allVars)/(double)(totalAfr + totalNonAfr)) * (double)totalAfr;

//         if (O == 0 || O_dash == 0 or E == 0){
//             continue;
//         }

//         double chi2 = std::pow(O - E, 2) / E;
//         double p_value = chi_square_p_value(chi2);
//         double frac_A = (double)O / (double)totalAfr;
//         double frac_nonA = (double)O_dash / (double)totalNonAfr;
//         double ratio = frac_nonA > 0 ? frac_A / frac_nonA : INFINITY;

//         // fout << std::get<0>(var) << "\t" 
//         //  << panmanUtils::getNucleotideFromCode(std::get<1>(var)) << "\t"
//         //  << std::get<2>(var) << "\t"
//         //  << O << "\t"
//         //  << O_dash << "\t"
//         //  << E << "\t"
//         //  << std::fixed << std::setprecision(4) << chi2 << "\t"
//         //  << std::scientific << std::setprecision(4) << p_value << "\t"
//         //  << std::fixed <<  std::setprecision(4) << frac_A << "\t"
//         //  << std::fixed <<  std::setprecision(4) << frac_nonA << "\t"
//         //  << std::fixed <<  std::setprecision(4) << ratio << std::endl;

//         if (finalVars.find(std::make_pair(std::get<0>(var)-1, std::get<1>(var))) == finalVars.end()){
//             finalVars[std::make_pair(std::get<0>(var), std::get<1>(var))]=std::make_tuple(std::get<0>(var), std::get<1>(var), std::get<2>(var), O, O_dash, E, chi2, p_value, frac_A, frac_nonA, ratio);
//         } else {
//             auto existingVar = finalVars[std::make_pair(std::get<0>(var)-1, std::get<1>(var))];
//             finalVars[std::make_pair(std::get<0>(var)-1, std::get<1>(var))]=std::make_tuple(std::get<0>(var)-1, std::get<1>(var), std::get<2>(existingVar)+std::get<2>(var), O, O_dash, E, chi2, p_value, frac_A, frac_nonA, ratio);
//         }
        
//     }

//     return;
//     // final round of compression
//     for (auto &entry: finalVars){
//         auto var = entry.second;
//         if (finalVars.find(std::make_pair(std::get<0>(var)+std::get<2>(var).size(), std::get<1>(var))) != finalVars.end()){
//             auto nextVar = finalVars[std::make_pair(std::get<0>(var)+std::get<2>(var).size(), std::get<1>(var))];
//             finalVars[std::make_pair(std::get<0>(var), std::get<1>(var))]=std::make_tuple(std::get<0>(var), std::get<1>(var), std::get<2>(var)+std::get<2>(nextVar), 
//                                                                                             std::get<3>(var), 
//                                                                                             std::get<4>(var), 
//                                                                                             std::get<5>(var), 
//                                                                                             std::get<6>(var), 
//                                                                                             std::get<7>(var), 
//                                                                                             std::get<8>(var), 
//                                                                                             std::get<9>(var), 
//                                                                                             std::get<10>(var));
//             finalVars.erase(std::make_pair(std::get<0>(var)+std::get<2>(var).size(), std::get<1>(var)));
//         }
//     }

//     // for (auto &entry: finalVars){
//     //     auto var = entry.second;
//     //     fout << std::get<0>(var) << "\t" 
//     //          << panmanUtils::getNucleotideFromCode(std::get<1>(var)) << "\t"
//     //          << std::get<2>(var) << "\t"
//     //          << std::get<3>(var) << "\t"
//     //          << std::get<4>(var) << "\t"
//     //          << std::get<5>(var) << "\t"
//     //          << std::fixed << std::setprecision(4) << std::get<6>(var) << "\t"
//     //          << std::scientific << std::setprecision(4) << std::get<7>(var) << "\t"
//     //          << std::fixed <<  std::setprecision(4) << std::get<8>(var) << "\t"
//     //          << std::fixed <<  std::setprecision(4) << std::get<9>(var) << "\t"
//     //          << std::fixed <<  std::setprecision(4) << std::get<10>(var) << std::endl;
//     // }
    
//     return;

// }


// Chi Square Test

/*
O = Number of Afr Sequeunce carryiung the mutation
E = Fraction of population carrying the allele * Total number of Afr sequences
0-E/E
Discard: 50% bases conserve use it as single alelle 
Look for SNP first and then marege.
*/