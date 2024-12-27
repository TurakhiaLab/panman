#include "panmanUtils.hpp"

int getTotalParsimonyParallelHelper(panmanUtils::Node* root, panmanUtils::NucMutationType nucMutType, panmanUtils::BlockMutationType blockMutType) {
    int totalMutations = 0;

    if(nucMutType != panmanUtils::NucMutationType::NNONE) {
        totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->nucMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
            for(int i = r.begin(); i != r.end(); i++) {
                if(((root->nucMutation[i].mutInfo) & 0x7) == nucMutType) {
                    if(nucMutType == panmanUtils::NucMutationType::NS) {
                        init += ((root->nucMutation[i].mutInfo) >> 4); // Length of contiguous mutation in case of substitution
                    } else {
                        init++;
                    }
                }
            }
            return init;
        }, [&](int x, int y) {
            return x + y;
        });
    }

    if(blockMutType == panmanUtils::BlockMutationType::BIn) {
        totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
            for(int i = r.begin(); i != r.end(); i++) {
                if(root->blockMutation[i].inversion == true) {
                    init++;
                }
            }
            return init;
        }, [&](int x, int y) {
            return x + y;
        });
    } else if(blockMutType != panmanUtils::BlockMutationType::NONE) {
        totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
            for(int i = r.begin(); i != r.end(); i++) {
                // If not an inversion and mut type matches. Inversion is marked by blockMutInfo = deletion and inversion = true
                if((blockMutType == panmanUtils::BlockMutationType::BI || root->blockMutation[i].inversion == false) && root->blockMutation[i].blockMutInfo == blockMutType) {
                    init++;
                }
            }
            return init;
        }, [&](int x, int y) {
            return x + y;
        });
    }

    totalMutations += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->children.size()), 0, [&](tbb::blocked_range<int>& r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++) {
            init += getTotalParsimonyParallelHelper(root->children[i], nucMutType, blockMutType);
        }
        return init;
    },
    [](int x, int y) -> int {
        return x+y;
    });

    return totalMutations;
}

int panmanUtils::Tree::getTotalParsimonyParallel(NucMutationType nucMutType, BlockMutationType blockMutType) {

    return getTotalParsimonyParallelHelper(root, nucMutType, blockMutType);

}

std::tuple<int, int, int> getBlockMutationsParallelHelper(panmanUtils::Node* root) {
    std::tuple<int, int, int> muts(0,0,0);

    std::get<0>(muts) += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++) {
            if(root->blockMutation[i].blockMutInfo == panmanUtils::BlockMutationType::BI) {
                init++;
            }
        }
        return init;
    }, [&](int x, int y) {
        return x + y;
    });

    std::get<1>(muts) += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++) {
            if(root->blockMutation[i].inversion == false && root->blockMutation[i].blockMutInfo == panmanUtils::BlockMutationType::BD) {
                init++;
            }
        }
        return init;
    }, [&](int x, int y) {
        return x + y;
    });
    
    std::get<2>(muts) += tbb::parallel_reduce(tbb::blocked_range<int>(0, root->blockMutation.size()), 0, [&](tbb::blocked_range<int> r, int init) -> int{
        for(int i = r.begin(); i != r.end(); i++) {
            if(root->blockMutation[i].inversion == true && root->blockMutation[i].blockMutInfo == panmanUtils::BlockMutationType::BD) {
                init++;
            }
        }
        return init;
    }, [&](int x, int y) {
        return x + y;
    });
    

    for (auto i=0; i< root->children.size(); i++) {
        std::tuple<int, int, int> child_muts = getBlockMutationsParallelHelper(root->children[i]);
        std::get<0>(muts) += std::get<0>(child_muts);
        std::get<1>(muts) += std::get<1>(child_muts);
        std::get<2>(muts) += std::get<2>(child_muts);
    }

    return muts;
}

std::tuple<int, int> getOtherBlockMutationsParallelHelper(
    panmanUtils::Node* root, 
    std::vector< bool >& blockExists,
    std::vector< bool >& blockStrand,
    std::vector<std::vector<uint32_t>>& dups,
    std::vector<uint32_t>& dupsPos) {
    
    std::tuple<int, int> muts(0,0);
    std::vector< bool > blockExistsParent = blockExists;
    // For reversing block mutations - primary block id, secondary block id, old mutation, old strand, new mutation, new strand
    std::vector< std::tuple< int32_t, bool, bool, bool, bool > > blockMutationInfo;

    // Block Mutations
    for(auto mutation: root->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        // std::vector<uint32_t> copies = map_[blocks[]]

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
                    std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
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

    // std::cout << (blockExists ==  blockExistsParent) << std::endl;

    for(auto mutation: root->blockMutation) {
        int32_t primaryBlockId = mutation.primaryBlockId;
        bool type = mutation.blockMutInfo;

        if (type==1) {
            std::vector<uint32_t> localDups = dups[dupsPos[primaryBlockId]];
            for (auto d: localDups){
                if (d!=primaryBlockId & blockExists[d] & blockExistsParent[d]) {
                    std::get<0>(muts) += 1;
                    break;
                }

                if (d!=primaryBlockId & !blockExists[d] & blockExistsParent[d]) {
                    std::get<1>(muts) += 1;
                    break;
                }
            }            
        }
    }

    for(panmanUtils::Node* child: root->children) {
        std::tuple<int, int> mutsChild = getOtherBlockMutationsParallelHelper(child, blockExists, blockStrand, dups, dupsPos);
        std::get<0>(muts) += std::get<0>(mutsChild);
        std::get<1>(muts) += std::get<1>(mutsChild);
    }

     // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        blockExists[std::get<0>(mutation)] = std::get<1>(mutation);
        blockStrand[std::get<0>(mutation)] = std::get<2>(mutation);
    }

    return muts;
}

struct VectorHash {
    std::size_t operator()(const std::vector<uint32_t>& vec) const {
        std::size_t hash = 0;
        for (uint32_t num : vec) {
            hash ^= std::hash<int>()(num) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

void panmanUtils::Tree::getBlockMutationsParallel() {
    //insertions, deletions, inversions
    std::tuple<int, int, int> muts = getBlockMutationsParallelHelper(root);
    std::cout << "Total Block Insertoins: " <<  std::get<0>(muts) << std::endl;
    std::cout << "Total Block Deletions: " <<  std::get<1>(muts) << std::endl;
    std::cout << "Total Block Inversion: " <<  std::get<2>(muts) << std::endl;

    // get duplicate blocks mapping (consensus to blockIDs)
    std::unordered_map<std::vector<uint32_t>, std::vector<uint32_t>, VectorHash> map_;
    for (auto i=0; i<blocks.size(); i++) {
        std::vector<uint32_t> consensus = blocks[i].consensusSeq;
        map_[consensus].push_back(blocks[i].primaryBlockId);
    }

    std::unordered_map<std::vector<uint32_t>, int, VectorHash> mapIndex;
    std::vector<std::vector<uint32_t>> dups(map_.size());
    std::vector<uint32_t> dupsPos(blocks.size());

    int index = 0;
    for(auto &a: map_){
        mapIndex[a.first] = index;
        for (auto &b: a.second){
            dupsPos[b] = index;
        }
        dups[index] = a.second;
        index++;
    }

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< bool >  blockExists(blocks.size(), false);
    std::vector< bool >  blockStrand(blocks.size(), true);

    std::tuple<int, int> otherMuts = getOtherBlockMutationsParallelHelper(root, blockExists, blockStrand, dups, dupsPos);
    std::cout << "Total Block Duplications: " <<  std::get<0>(otherMuts) << std::endl;
    std::cout << "Total Block Translocation: " <<  std::get<1>(otherMuts) << std::endl;
}

void panmanUtils::Tree::printSummary(std::ostream &out) {

    out << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    out << "Total Samples in Tree: " << m_numLeaves << std::endl;
    out << "Total Substitutions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NS) << std::endl;
    out << "Total Insertions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NI, panmanUtils::BlockMutationType::BI) << std::endl;
    out << "Total Deletions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::ND, panmanUtils::BlockMutationType::BD) << std::endl;
    out << "Total Inversions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NNONE, panmanUtils::BlockMutationType::BIn) << std::endl;
    // out << "Total SNP Substitutions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPS) << std::endl;
    // out << "Total SNP Insertions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPI) << std::endl;
    // out << "Total SNP Deletions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPD) << std::endl;
    out << "Max Tree Depth: " << m_maxDepth << std::endl;
    out << "Mean Tree Depth: " << m_meanDepth << std::endl;

    getBlockMutationsParallel();

}
