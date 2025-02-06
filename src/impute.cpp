
#include "panmanUtils.hpp"

std::vector<bool> posWithN(const panmanUtils::NucMut mut) {
    int mutLength = mut.mutInfo >> 4;
    std::vector<bool> isN = std::vector<bool>(mutLength);
    for(int i = 0; i < mutLength; i++) {
        // Peel away layers to extract a single nucleotide
        isN[i] = (panmanUtils::getNucleotideFromCode((mut.nucs >> (4*(5-i))) & 0xF) == 'N');
    }
    return isN;
}

void panmanUtils::Tree::imputeNs() {
    std::cout << "Imputing a tree" << std::endl;

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

    // We only care about N/missing nucleotides
    int codeForN = panmanUtils::getCodeFromNucleotide('N');

    for (const auto& curMut: node->nucMutation) {
        // Based on printSingleNodeHelper() in fasta.cpp
        uint32_t type = (curMut.mutInfo & 0x7);
        // <3 NX types may be multibase, >= 3 NSNPX types are all one base
        int len = type < 3 ? (curMut.mutInfo >> 4) : 1;

        if (type == panmanUtils::NucMutationType::NSNPS
            || type == panmanUtils::NucMutationType::NS) {
            std::vector<bool> isN = posWithN(curMut);
            // If there's an N somewhere, store this SNP
            if (std::find(begin(isN), end(isN), true) != end(isN)) {
                substitutions.push_back(std::make_pair(node, curMut));
            }
        } else if (type == panmanUtils::NucMutationType::NSNPI
                   || type == panmanUtils::NucMutationType::NI) {
            if (insertions.empty() || insertions.back().first != node
                || !insertions.back().second.back().samePosExceptGap(curMut)) {
                // This insertion can't be merged with the previous one, so make a new entry
                insertions.push_back(std::make_pair(node, std::vector<panmanUtils::NucMut>()));
            }
            // In all cases, add this insertion to the currently growing insertion
            insertions.back().second.push_back(curMut);
        }
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
    if ((mutToN.mutInfo & 0x7) == panmanUtils::NucMutationType::NS) {
        // Add non-N mutations back in (for MNPs which are partially N)
        std::vector<bool> isN = posWithN(mutToN);
        for(int i = 0; i < isN.size(); i++) {
            if (!isN[i]) {
                node->nucMutation.push_back(NucMut(mutToN, i));
            }
        }
    }
}

void panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, std::vector<panmanUtils::NucMut> mutToN) {
    std::cout << "Imputing indel (length " << mutToN.size() << ") for " << node->identifier << " pos (" << mutToN[0].primaryBlockId;
    std::cout << ", " << mutToN[0].nucPosition << ", " << mutToN[0].nucGapPosition << ")" << std::endl;

    // Algorithm only works if niblings are available
    if (node == nullptr || node->parent == nullptr) return;

    Node* nibling = nullptr;
    //std::vector<panmanUtils::NucMut> niblingMut = std::vector<panmanUtils::NucMut>(mutToN.size());
    // Find nibling with insertion of identical length/position
    // Erase insertion from current node
    // Erase insertion from nibling
    // Add insertion to parent (likely the exact one from nibling)
    // Add deletion to nibling's sibling (same as one from nibling, but changed to deletion)
}