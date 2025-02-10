
#include "panmanUtils.hpp"

// Convenience function to check if a the ith base in mut.nucs is N
bool isNucN(panmanUtils::NucMut mut, int offset) {
    // Peel away layers to extract a single nucleotide
    int curNucCode = (mut.nucs >> (4*(5-offset))) & 0xF;
    return (panmanUtils::getNucleotideFromCode(curNucCode) == 'N');
}

void panmanUtils::Tree::imputeNs() {
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > > substitutions;
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > > insertions;
    std::unordered_map< panmanUtils::IndelPosition, std::unordered_set<std::string> > allInsertions;
    findMutationsToN(root, substitutions, insertions, allInsertions);

    std::cout << substitutions.size() << " substitutions to impute" << std::endl;
    for (const auto& toImpute: substitutions) {
        imputeSNV(toImpute.first, toImpute.second);
    }
    
    std::cout << insertions.size() << " insertions to impute" << std::endl;
    for (const auto& toImpute: insertions) {
        imputeInsertion(toImpute.first, toImpute.second, allInsertions);
    }
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > >& substitutions,
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > >& insertions,
        std::unordered_map< panmanUtils::IndelPosition, std::unordered_set<std::string> >& allInsertions) {
    if (node == nullptr) return;

    // Default value
    panmanUtils::IndelPosition curInsertion = panmanUtils::IndelPosition();

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
                substitutions.emplace_back(node, curMut);
            }
        } else if (type == panmanUtils::NucMutationType::NSNPI
                   || type == panmanUtils::NucMutationType::NI) {
            if (curInsertion.indelLength == -1) {
                curInsertion = panmanUtils::IndelPosition(curMut, hasNs);
            } else if (!curInsertion.mergeIndels(curMut, hasNs)) {
                // This insertion can't be merged with the previous one
                if (curInsertion.hasSpecialNucs) {
                    insertions.emplace_back(node, curInsertion);
                }
                if (allInsertions.find(curInsertion) == allInsertions.end()) {
                    allInsertions.emplace(curInsertion, std::unordered_set<std::string>());
                }
                allInsertions[curInsertion].emplace(node->identifier);

                curInsertion = panmanUtils::IndelPosition(curMut, hasNs);
            }
        }
    }

    // Handle last insertion just in case
    if (curInsertion.indelLength != -1) {
        if (allInsertions.find(curInsertion) == allInsertions.end()) {
            allInsertions.emplace(curInsertion, std::unordered_set<std::string>());
        }
        allInsertions[curInsertion].emplace(node->identifier);
        if (curInsertion.hasSpecialNucs) {
            insertions.emplace_back(node, curInsertion);
        }
    }

    for(auto child: node->children) {
        findMutationsToN(child, substitutions, insertions, allInsertions);
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

void panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, panmanUtils::IndelPosition mutToN,
    std::unordered_map< panmanUtils::IndelPosition, std::unordered_set<std::string> >& allInsertions) {
    if (node == nullptr) return;

    if (findNearbyInsertion(node->parent, mutToN, 10 - node->branchLength, node, allInsertions)) {
        std::cout << "Found insertion" << std::endl;
    } else {
        std::cout << "No insertion" << std::endl;
    }
}

panmanUtils::Node* panmanUtils::Tree::findNearbyInsertion(
    panmanUtils::Node* node, panmanUtils::IndelPosition mutToN, int allowedDistance, panmanUtils::Node* ignore,
    std::unordered_map< panmanUtils::IndelPosition, std::unordered_set<std::string> >& allInsertions) {
    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nullptr;

    panmanUtils::IndelPosition nonNMut = mutToN;
    if (allInsertions.find(nonNMut) != allInsertions.end() && 
        allInsertions[nonNMut].find(node->identifier) != allInsertions[nonNMut].end()) {
        return node;
    }

    // Try children
    panmanUtils::Node* found = nullptr;
    for (const auto& child: node->children) {
        if (child != ignore) {
            found = findNearbyInsertion(child, mutToN, allowedDistance - child->branchLength, node, allInsertions);
            if (found != nullptr) return found;
        }
    }
    // Try parent
    if (node->parent != ignore) {
        found = findNearbyInsertion(node->parent, mutToN, allowedDistance - node->branchLength, node, allInsertions);
    }
    return found;
}