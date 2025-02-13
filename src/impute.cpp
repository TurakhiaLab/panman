
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > > substitutions;
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > > insertions;
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> > allInsertions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > > originalNucs;
    findMutationsToN(root, substitutions, insertions, allInsertions, originalNucs);

    std::cout << substitutions.size() << " substitutions to impute" << std::endl;
    for (const auto& toImpute: substitutions) {
        imputeSNV(toImpute.first, toImpute.second);
    }
    
    std::cout << insertions.size() << " insertions to impute" << std::endl;
    for (const auto& toImpute: insertions) {
        imputeInsertion(toImpute.first, toImpute.second, allowedIndelDistance, allInsertions, originalNucs);
    }
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > >& substitutions,
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > >& insertions,
        std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return;

    std::vector<panmanUtils::IndelPosition> nodeInsertions;
    bool curInsertionHasNs = false;

    for (const auto& curMut: node->nucMutation) {
        // Does this mutation have Ns?
        bool hasNs = false;
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNuc = curMut.getNucCode(i);
            hasNs |= (curNuc == panmanUtils::NucCode::N);
            if (curMut.isDeletion() || curMut.isSubstitution()) {
                // set original nucleotide
            }
        }

        // Save mutation if relevant
        if (curMut.isSubstitution()) {
            if (hasNs) {
                substitutions.emplace_back(node, curMut);
            }
        } else if (curMut.isInsertion()) {
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
        findMutationsToN(child, substitutions, insertions, allInsertions, originalNucs);
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
            if (mutToN.getNucCode(i) != panmanUtils::NucCode::N) {
                node->nucMutation.push_back(NucMut(mutToN, i));
            }
        }
    }
}

void panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, panmanUtils::IndelPosition mutToN, int allowedDistance,
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return;

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    MutationList bestMutationList = MutationList();

    for (const auto& nearby: findNearbyInsertions(node, mutToN, allowedDistance, nullptr, allInsertions, originalNucs)) {
        std::cout << nearby.first->identifier << std::endl;
        if (nearby.first != node) {
            panmanUtils::MutationList simpleMutations = nearby.second.copy();
            simpleMutations.nucMutation = consolidateNucMutations(simpleMutations.nucMutation);
            simpleMutations.blockMutation = consolidateBlockMutations(simpleMutations.blockMutation);
            // Ignore block mutations for now
            if (simpleMutations.blockMutation.size() > 0) continue;

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
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {

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
            for (const auto& nearby: findNearbyInsertions(child, mutToN, allowedDistance - child->branchLength, 
                                                          node, allInsertions, originalNucs)) {
                // Add mutations to get to child (which must be reversed)
                panmanUtils::MutationList toAdd = panmanUtils::MutationList(child, true, originalNucs[child->identifier]);
                nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(toAdd));
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutToN, allowedDistance - node->branchLength,
                                                      node, allInsertions, originalNucs)) {
            // Add mutations to get to parent
            panmanUtils::MutationList toAdd = panmanUtils::MutationList(node, false, originalNucs[node->identifier]);
            nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(toAdd));
        }
    }
    std::cout << nearbyInsertions.size() << " possibilities from " << node->identifier << std::endl;
    return nearbyInsertions;
}

void panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList pathMutations) {
    std::cout << "Moving " << toMove->identifier << " to be a child of " << newParent->identifier << std::endl;
    
    if (toMove->parent != nullptr) {
        std::vector<Node*>::iterator position = std::find(toMove->parent->children.begin(), toMove->parent->children.end(), toMove);
        toMove->parent->children.erase(position);
    }
    newParent->children.emplace_back(toMove);

    toMove->branchLength = 1;
    toMove->level = newParent->level + 1;
    toMove->parent = newParent;
    toMove->nucMutation = pathMutations.nucMutation;
    toMove->blockMutation = pathMutations.blockMutation;
    toMove->annotations = std::vector<std::string>();
    toMove->isComMutHead = false;
}