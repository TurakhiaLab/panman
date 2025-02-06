
#include "panmanUtils.hpp"

// Convenience function to check if a the ith base in mut.nucs is N
bool isNucN(panmanUtils::NucMut mut, int offset) {
    // Peel away layers to extract a single nucleotide
    int curNucCode = (mut.nucs >> (4*(5-offset))) & 0xF;
    return (panmanUtils::getNucleotideFromCode(curNucCode) == 'N');
}

void panmanUtils::Tree::imputeNs() {
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > > substitutions;
    std::vector< std::pair< panmanUtils::Node*, std::vector<panmanUtils::NucMut> > > insertions;
    findMutationsToN(root, substitutions, insertions);

    std::cout << substitutions.size() << " substitutions to impute" << std::endl;
    for (const auto& toImpute: substitutions) {
        imputeSNV(toImpute.first, toImpute.second);
    }
    
    std::cout << insertions.size() << " insertions to impute" << std::endl;
    for (const auto& toImpute: insertions) {
        imputeInsertion(toImpute.first, toImpute.second);
    }
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > >& substitutions,
        std::vector< std::pair< panmanUtils::Node*, std::vector<panmanUtils::NucMut> > >& insertions) {
    if (node == nullptr) return;

    std::vector<panmanUtils::NucMut> curInsertion;
    bool curInsertionHasNs = false;

    for (const auto& curMut: node->nucMutation) {
        // Based on printSingleNodeHelper() in fasta.cpp
        uint32_t type = curMut.type();

        // Does this mutation have Ns?
        bool hasNs = false;
        for(int i = 0; i < curMut.length(); i++) {
            hasNs |= isNucN(curMut, i);
        }

        // Save mutation if relevant
        if (type == panmanUtils::NucMutationType::NSNPS
            || type == panmanUtils::NucMutationType::NS) {
            if (hasNs) {
                substitutions.push_back(std::make_pair(node, curMut));
            }
        } else if (type == panmanUtils::NucMutationType::NSNPI
                   || type == panmanUtils::NucMutationType::NI) {
            if (curInsertion.empty() || !curInsertion.back().samePosExceptGap(curMut)) {
                // This insertion can't be merged with the previous one
                if (curInsertionHasNs) {
                    insertions.push_back(std::make_pair(node, curInsertion));
                }

                curInsertion = std::vector<NucMut>();
                curInsertionHasNs = false;
            }
            
            // In all cases, add this insertion to the currently growing insertion
            curInsertion.push_back(curMut);
            curInsertionHasNs |= hasNs;
        }
    }

    // Handle last insertion just in case
    if (curInsertionHasNs) {
        insertions.push_back(std::make_pair(node, curInsertion));
    }

    for(auto child: node->children) {
        findMutationsToN(child, substitutions, insertions);
    }
}

void panmanUtils::Tree::imputeSNV(panmanUtils::Node* node, NucMut mutToN) {
    if (node == nullptr) return;

    // Get rid of the old mutation in the node's list
    std::vector<NucMut>::iterator oldIndex = std::find(node->nucMutation.begin(), node->nucMutation.end(), mutToN);
    node->nucMutation.erase(oldIndex);

    // Possible MNP
    if (mutToN.type() == panmanUtils::NucMutationType::NS) {
        // Add non-N mutations back in (for MNPs which are partially N)
        for(int i = 0; i < mutToN.length(); i++) {
            if (isNucN(mutToN, i)) {
                node->nucMutation.push_back(NucMut(mutToN, i));
            }
        }
    }
}

void panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, std::vector<panmanUtils::NucMut> mutToN) {
    // Algorithm only works if niblings are available
    if (node == nullptr || node->parent == nullptr) return;

    // Determine total insertion length (NI may be multibase, NSNPI is one base)
    int imputeLength = 0;
    for (const auto& curMut: mutToN) imputeLength += curMut.length();

    std::cout << "Imputing indel (length " << imputeLength << ") for " << node->identifier << " pos (" << mutToN[0].primaryBlockId;
    std::cout << ", " << mutToN[0].nucPosition << ", " << mutToN[0].nucGapPosition << ")" << std::endl;

    Node* sourceNibling = nullptr;
    std::vector<panmanUtils::NucMut> niblingMut;

    // Find nibling with insertion of identical length/position
    for (const auto& sibling: node->parent->children) {
        // Don't look at the current node's children
        if (sibling != node) {
            for (const auto& nibling: sibling->children) {
                // Temporary storage for the correct mutation
                std::vector<panmanUtils::NucMut> curInsertion;
                int curInsertionLength = 0;

                for (const auto& curMut: nibling->nucMutation) {
                    uint32_t type = curMut.type();
                    bool isInsertion = (type == panmanUtils::NucMutationType::NSNPI 
                                        || type == panmanUtils::NucMutationType::NI);
                    if (isInsertion && curMut.samePosExceptGap(mutToN[0])) {
                        curInsertion.push_back(curMut);
                        curInsertionLength += curMut.length();
                    }
                }

                if (curInsertionLength == imputeLength) {
                    sourceNibling = nibling;
                    niblingMut = curInsertion; 
                }
            }
        }
    }

    // Failed to find an appropriate nibling
    if (sourceNibling == nullptr) {
        return;
    }
    std::cout << "Nibling found: " << sourceNibling->identifier << std::endl;
    // ahh have to also handle the aunt/uncle case
    // Erase insertion from current node
    // Erase insertion from nibling
    // Add insertion to parent (likely the exact one from nibling)
    // Add deletion to nibling's sibling (same as one from nibling, but changed to deletion)
}