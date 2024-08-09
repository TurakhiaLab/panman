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

void panmanUtils::Tree::printSummary(std::ostream &out) {

    out << "Total Nodes in Tree: " << m_currInternalNode + m_numLeaves << std::endl;
    out << "Total Samples in Tree: " << m_numLeaves << std::endl;
    out << "Total Substitutions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NS) << std::endl;
    out << "Total Insertions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NI, panmanUtils::BlockMutationType::BI) << std::endl;
    out << "Total Deletions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::ND, panmanUtils::BlockMutationType::BD) << std::endl;
    out << "Total Inversions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NNONE, panmanUtils::BlockMutationType::BIn) << std::endl;
    out << "Total SNP Substitutions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPS) << std::endl;
    out << "Total SNP Insertions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPI) << std::endl;
    out << "Total SNP Deletions: " << getTotalParsimonyParallel(panmanUtils::NucMutationType::NSNPD) << std::endl;
    out << "Max Tree Depth: " << m_maxDepth << std::endl;
    out << "Mean Tree Depth: " << m_meanDepth << std::endl;

}
