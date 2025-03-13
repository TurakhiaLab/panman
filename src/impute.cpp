
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

bool panmanUtils::canImpute(int8_t oldNuc, int8_t newNuc) {
    // Using table from https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    switch(newNuc) {
    case panmanUtils::NucCode::Y: // C or T
        return oldNuc == panmanUtils::NucCode::C || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::R: // A or G
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::G;
    case panmanUtils::NucCode::W: // A or T
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::S: // C or G
        return oldNuc == panmanUtils::NucCode::C || oldNuc == panmanUtils::NucCode::G;
    case panmanUtils::NucCode::K: // G or T
        return oldNuc == panmanUtils::NucCode::G || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::M: // A or C
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::C;
    case panmanUtils::NucCode::B: // C, G, or T
        return oldNuc == panmanUtils::NucCode::C || oldNuc == panmanUtils::NucCode::G || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::D: // A, G, or T
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::G || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::H: // A, C, or T
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::C || oldNuc == panmanUtils::NucCode::T;
    case panmanUtils::NucCode::V: // A, C, or G
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::C || oldNuc == panmanUtils::NucCode::G;
    case panmanUtils::NucCode::N: // Any nucleotide (A, C, G, T)
        return oldNuc == panmanUtils::NucCode::A || oldNuc == panmanUtils::NucCode::C || 
               oldNuc == panmanUtils::NucCode::G || oldNuc == panmanUtils::NucCode::T;
    default:
        return false;
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