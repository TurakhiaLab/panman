
#include "panmanUtils.hpp"

// Convenience function to check if a the ith base in mut.nucs is N
bool isNucN(panmanUtils::NucMut mut, int offset) {
    // Peel away layers to extract a single nucleotide
    int curNucCode = (mut.nucs >> (4*(5-offset))) & 0xF;
    return (panmanUtils::getNucleotideFromCode(curNucCode) == 'N');
}

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > > substitutions;
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > > insertions;
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> > allInsertions;
    findMutationsToN(root, substitutions, insertions, allInsertions);

    std::cout << substitutions.size() << " substitutions to impute" << std::endl;
    for (const auto& toImpute: substitutions) {
        imputeSNV(toImpute.first, toImpute.second);
    }
    
    std::cout << insertions.size() << " insertions to impute" << std::endl;
    for (const auto& toImpute: insertions) {
        imputeInsertion(toImpute.first, toImpute.second, allowedIndelDistance, allInsertions);
    }
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > >& substitutions,
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > >& insertions,
        std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions) {
    if (node == nullptr) return;

    std::vector<panmanUtils::IndelPosition> nodeInsertions;
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
                substitutions.emplace_back(node, curMut);
            }
        } else if (type == panmanUtils::NucMutationType::NSNPI
                   || type == panmanUtils::NucMutationType::NI) {
            if (nodeInsertions.empty()) {
                nodeInsertions.emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionHasNs = hasNs;
            } else if (nodeInsertions.back().mergeIndels(curMut)) {
                // Current insertion was merged with the previous one.
                curInsertionHasNs |= hasNs;
            } else {
                // Start new insertion
                if (curInsertionHasNs) {
                    insertions.emplace_back(node, nodeInsertions.back());
                }

                nodeInsertions.emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionHasNs = hasNs;
            }
        }
    }

    // Handle last insertion just in case
    if (curInsertionHasNs) {
        insertions.emplace_back(node, nodeInsertions.back());
    }
    // Convert vector to set for easier lookup
    allInsertions.emplace(node->identifier, 
        std::unordered_set<panmanUtils::IndelPosition>(nodeInsertions.begin(), nodeInsertions.end()));

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

// Simplify the changing mutations, e.g. cancel out insertions and deletions
// TODO: implement
panmanUtils::MutationList simplifyMutations(panmanUtils::MutationList pathMutations) {
    return pathMutations;
}

void panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, panmanUtils::IndelPosition mutToN, int allowedDistance,
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions) {
    if (node == nullptr) return;

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    MutationList bestMutationList = MutationList();

    for (const auto& nearby: findNearbyInsertions(node, mutToN, allowedDistance, nullptr, allInsertions)) {
        std::cout << nearby.first->identifier << std::endl;
        if (nearby.first != node) {
            panmanUtils::MutationList simpleMutations = simplifyMutations(nearby.second);

            // Parsimony improvement score is the decrease in mutation count
            int blockImprovement = node->blockMutation.size() - simpleMutations.blockMutation.size();
            int nucImprovement = node->nucMutation.size() - simpleMutations.nucMutation.size();
            int totalImprovement = blockImprovement + nucImprovement;
            std::cout << blockImprovement << ", " << nucImprovement << std::endl;

            if (blockImprovement >= 0 && nucImprovement >= 0 && totalImprovement > bestParsimonyImprovement) {
                bestParsimonyImprovement = blockImprovement + nucImprovement;
                bestNewParent = nearby.first;
                bestMutationList = simpleMutations;
            }
        }
    }

    if (bestNewParent != nullptr) {
        moveNode(node, bestNewParent, bestMutationList);
    } else {
        std::cout << "Failed to find insertion" << std::endl;
    }
}

const std::vector< std::pair< panmanUtils::Node*, panmanUtils::MutationList > > panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, panmanUtils::IndelPosition mutToN, int allowedDistance, panmanUtils::Node* ignore,
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions) {

    std::vector< std::pair< panmanUtils::Node*, panmanUtils::MutationList > > nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    panmanUtils::IndelPosition nonNMut = mutToN;
    std::string curID = node->identifier;
    if (allInsertions[curID].find(mutToN) != allInsertions[curID].end()) {
        nearbyInsertions.emplace_back(node, panmanUtils::MutationList());
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore) {
            for (const auto& nearby: findNearbyInsertions(child, mutToN, allowedDistance - child->branchLength, node, allInsertions)) {
                // Add mutations to get to child
                nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(MutationList(child, false)));
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutToN, allowedDistance - node->branchLength, node, allInsertions)) {
            // Add mutations to get to parent (which must be reversed)
            nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(MutationList(node, true)));
        }
    }
    std::cout << nearbyInsertions.size() << " possibilities from " << node->identifier << std::endl;
    return nearbyInsertions;
}

void panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList pathMutations) {
    std::cout << "Moving " << toMove->identifier << " to be a child of " << newParent->identifier << std::endl;

    // TODO: implement
}