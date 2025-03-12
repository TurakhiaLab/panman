
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs() {
    int imputedNs = 0;
    for (const auto& child: root->children) {
        imputeSubtree(root, imputedNs);
    }

    std::cout << "Imputed " << imputedNs << " SNPs to N" << std::endl;
}

void panmanUtils::Tree::imputeSubtree(panmanUtils::Node* node, int& imputedNs) {
    imputedNs += imputeAllSubstitutionsWithNs(node->nucMutation);
    for (const auto& child: node->children) {
        imputeSubtree(child, imputedNs);
    }
}

const int panmanUtils::Tree::imputeSubstitution(std::vector<panmanUtils::NucMut>& nucMutation, const NucMut& mutToN) {
    // Get rid of the old mutation in the node's list
    std::vector<NucMut>::iterator oldIndex = std::find(nucMutation.begin(), nucMutation.end(), mutToN);
    oldIndex = nucMutation.erase(oldIndex);
    int subNs = mutToN.length();

    // Possible MNP
    if (mutToN.type() == panmanUtils::NucMutationType::NS) {
        std::vector<NucMut> snps;
        // Add non-N mutations back in (for MNPs which are partially N)
        for(int i = 0; i < mutToN.length(); i++) {
            if (mutToN.getNucCode(i) != panmanUtils::NucCode::N) {
                snps.push_back(NucMut(mutToN, i));
            }
        }
        nucMutation.insert(oldIndex, snps.begin(), snps.end());
        // These SNPs were not erased
        subNs -= snps.size();
    }
    return subNs;
}

const int panmanUtils::Tree::imputeAllSubstitutionsWithNs(std::vector<panmanUtils::NucMut>& nucMutation) {
    int totalImputedNs = 0;
    // Loop over vector backwards as elements will be erased
    for (int i = nucMutation.size() - 1; i >= 0; i--) {
        panmanUtils::NucMut curMut = nucMutation[i];

        if (curMut.isSubstitution()) {
            for (int i = 0; i < curMut.length(); i++) {
                if (curMut.getNucCode(i) == panmanUtils::NucCode::N) {
                    totalImputedNs += imputeSubstitution(nucMutation, curMut);
                    break;
                }
            }
        }
    }
    return totalImputedNs;
}