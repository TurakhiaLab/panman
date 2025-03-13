
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs() {
    int imputedNs = 0;
    for (const auto& child: root->children) {
        imputeSubtree(root, imputedNs);
    }

    std::string pluralS = (imputedNs == 1) ? "" : "s";
    std::cout << "Imputed " << imputedNs << " SNP" << pluralS << " to N" << std::endl;
}

void panmanUtils::Tree::imputeSubtree(panmanUtils::Node* node, int& imputedNs) {
    imputedNs += imputeSubstitutions(node->nucMutation);
    for (const auto& child: node->children) {
        imputeSubtree(child, imputedNs);
    }
}

int panmanUtils::imputeSubstitutions(std::vector<panmanUtils::NucMut>& nucMutation) {
    int totalImputedNs = 0;
    // Will copy mutations back from oldMuts to nucMutation after processing
    std::vector<panmanUtils::NucMut> oldMuts = nucMutation;
    nucMutation.clear();

    for (const auto& curMut: oldMuts) {
        if (curMut.isSubstitution()) {
            // Store non-N substitutions seen so far, which must be retained
            panmanUtils::NucMut newMut = panmanUtils::NucMut(curMut, 0);

            for (int i = 0; i < curMut.length(); i++) {
                int8_t curNucCode = curMut.getNucCode(i);
                if (curNucCode == panmanUtils::NucCode::N) {
                    // Add previous non-N MNP if present
                    if (newMut.length() != 0) nucMutation.push_back(newMut);

                    // Reset non-N MNP
                    newMut = panmanUtils::NucMut(curMut, i);
                    totalImputedNs++;
                } else {
                    // Must save this non-N substitution
                    newMut.appendNuc(curNucCode);
                }
            }

            // Add non-N MNP if present
            if (newMut.length() != 0) nucMutation.push_back(newMut);
        } else {
            // Non-substitutions are not imputed
            nucMutation.push_back(curMut);
        }
    }
    
    return totalImputedNs;
}