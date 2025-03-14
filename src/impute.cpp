
#include "panmanUtils.hpp"

void panmanUtils::Tree::impute() {
    // Prepare current-state trackers
    sequence_t sequence;
    blockExists_t blockExists;
    blockStrand_t blockStrand;
    getSequenceFromReference(sequence, blockExists, blockStrand, root->identifier);

    int imputedBases = 0;
    for (const auto& child: root->children) {
        imputeSubtree(root, sequence, imputedBases);
    }

    std::string pluralS = (imputedBases == 1) ? "" : "s";
    std::cout << "Imputed " << imputedBases << " SNP" << pluralS << std::endl;
}

void panmanUtils::Tree::imputeSubtree(panmanUtils::Node* node, sequence_t& sequence, int& imputedBases) {
    imputedBases += imputeSubstitutions(node->nucMutation, sequence);

    // Only bother applying mutations etc. if there's more subtree to search
    if (!node->children.empty()) {
        auto mutationInfo = panmanUtils::applyNucMut(node->nucMutation, sequence);

        for (const auto& child: node->children) {
            imputeSubtree(child, sequence, imputedBases);
        }

        panmanUtils::undoNucMut(sequence, mutationInfo);
    }
}

bool panmanUtils::canImpute(char oldNuc, int8_t newNuc) {
    // Using table from https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    switch(newNuc) {
    case panmanUtils::NucCode::Y: // C or T
        return oldNuc == 'C' || oldNuc == 'T';
    case panmanUtils::NucCode::R: // A or G
        return oldNuc == 'A' || oldNuc == 'G';
    case panmanUtils::NucCode::W: // A or T
        return oldNuc == 'A' || oldNuc == 'T';
    case panmanUtils::NucCode::S: // C or G
        return oldNuc == 'C' || oldNuc == 'G';
    case panmanUtils::NucCode::K: // G or T
        return oldNuc == 'G' || oldNuc == 'T';
    case panmanUtils::NucCode::M: // A or C
        return oldNuc == 'A' || oldNuc == 'C';
    case panmanUtils::NucCode::B: // C, G, or T
        return (oldNuc == 'C' || oldNuc == 'G' || oldNuc == 'T' 
             || oldNuc == 'Y' || oldNuc == 'S' || oldNuc == 'K');
    case panmanUtils::NucCode::D: // A, G, or T
        return (oldNuc == 'A' || oldNuc == 'G' || oldNuc == 'T' 
             || oldNuc == 'R' || oldNuc == 'W' || oldNuc == 'K');
    case panmanUtils::NucCode::H: // A, C, or T
        return (oldNuc == 'A' || oldNuc == 'C' || oldNuc == 'T' 
             || oldNuc == 'Y' || oldNuc == 'W' || oldNuc == 'M');
    case panmanUtils::NucCode::V: // A, C, or G
        return (oldNuc == 'A' || oldNuc == 'C' || oldNuc == 'G' 
             || oldNuc == 'R' || oldNuc == 'S' || oldNuc == 'M');
    case panmanUtils::NucCode::N: // Any nucleotide
        return oldNuc != '-';
    default:
        return false;
    }
}

int panmanUtils::imputeSubstitutions(std::vector<panmanUtils::NucMut>& nucMutation, const sequence_t& sequence) {
    int totalImputedBases = 0;
    // Will copy mutations back from oldMuts to nucMutation after processing
    std::vector<panmanUtils::NucMut> oldMuts = nucMutation;
    nucMutation.clear();

    for (const auto& curMut: oldMuts) {
        if (curMut.isSubstitution()) {
            // Store non-impute substitutions seen so far, which must be retained
            panmanUtils::NucMut newMut = panmanUtils::NucMut(curMut, 0);

            for (int i = 0; i < curMut.length(); i++) {
                int8_t curNucCode = curMut.getNucCode(i);
                panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

                if (panmanUtils::canImpute(curPos.getSequenceBase(sequence), curNucCode)) {
                    // Add previous non-imputed MNP if present
                    if (newMut.length() != 0) nucMutation.push_back(newMut);

                    // Reset non-imputed MNP
                    newMut = panmanUtils::NucMut(curMut, i + 1);
                    totalImputedBases++;
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
    
    return totalImputedBases;
}