
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > > substitutions;
    std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > > insertions;
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> > allInsertions;
    std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher >  curBlockSeqs;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > > originalNucs;
    
    for (const auto& curBlock: blocks) {
        for (int i = 0; i < curBlock.consensusSeq.size(); i++) {
            bool endFlag = false;
            for (int j = 0; j < 8; j++) {
                int curNucCode = ((curBlock.consensusSeq[i] >> (4*(7 - j))) & 15);
                if(curNucCode == 0) {
                    endFlag = true;
                    break;
                } else {
                    curBlockSeqs[Coordinate(i*8 + j, 0, curBlock.primaryBlockId, curBlock.secondaryBlockId)] = (int8_t) curNucCode;
                }
            }

            if(endFlag) break;
        }
    }
    
    findMutationsToN(root, substitutions, insertions, allInsertions, curBlockSeqs, originalNucs);

    for (const auto& toImpute: substitutions) {
        imputeSNV(toImpute.first, toImpute.second);
    }
    std::cout << "Imputed " << substitutions.size() << "/" << substitutions.size() << " SNPs/MNPs to N" << std::endl;
    
    size_t insertionImputeSuccesses = 0;
    for (const auto& toImpute: insertions) {
        insertionImputeSuccesses += imputeInsertion(toImpute.first, toImpute.second, allowedIndelDistance, allInsertions, originalNucs);
    }
    std::cout << "Imputed " << insertionImputeSuccesses << "/" << insertions.size() << " insertions to N" << std::endl;
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::NucMut > >& substitutions,
        std::vector< std::pair< panmanUtils::Node*, panmanUtils::IndelPosition > >& insertions,
        std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
        std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher >& curBlockSeqs,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return;

    std::vector<panmanUtils::IndelPosition> nodeInsertions;
    bool curInsertionHasNs = false;

    for (const auto& curMut: node->nucMutation) {
        // Does this mutation have Ns?
        bool hasNs = false;
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            hasNs |= (curNucCode == panmanUtils::NucCode::N);

            if (curBlockSeqs.find(curPos) == curBlockSeqs.end()) {
                curBlockSeqs[curPos] = panmanUtils::NucCode::MISSING;
            }
            originalNucs[node->identifier][curPos] = curBlockSeqs[curPos];
            curBlockSeqs[curPos] = curNucCode;
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
        findMutationsToN(child, substitutions, insertions, allInsertions, curBlockSeqs, originalNucs);
    }

    // Undo mutations before passing back up the tree
    for (const auto& curMut: node->nucMutation) {
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            curBlockSeqs[curPos] = originalNucs[node->identifier][curPos];
        }
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

bool panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, panmanUtils::IndelPosition mutToN, int allowedDistance,
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return false;

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    MutationList bestMutationList = MutationList();

    for (const auto& nearby: findNearbyInsertions(node, mutToN, allowedDistance, nullptr, allInsertions, originalNucs)) {
        if (nearby.first != node) {
            panmanUtils::MutationList simpleMutations = nearby.second.copy();
            std::reverse(simpleMutations.nucMutation.begin(), simpleMutations.nucMutation.end());
            simpleMutations.nucMutation = consolidateNucMutations(simpleMutations.nucMutation);
            simpleMutations.blockMutation = consolidateBlockMutations(simpleMutations.blockMutation);
            // Ignore block mutations for now
            if (simpleMutations.blockMutation.size() > 0) continue;

            // Parsimony improvement score is the decrease in mutation count
            int nucImprovement = 0;
            for (const auto& curMut: simpleMutations.nucMutation) nucImprovement -= curMut.length();
            for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();

            if (nucImprovement > bestParsimonyImprovement) {
                bestParsimonyImprovement = nucImprovement;
                bestNewParent = nearby.first;
                bestMutationList = simpleMutations;
            }
        }
    }

    if (bestNewParent != nullptr) {
        moveNode(node, bestNewParent, bestMutationList);
        return true;
    } else {
        return false;
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
    return nearbyInsertions;
}

void panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList pathMutations) {
    // Make dummy parent from grandparent -> dummy -> newParent
    panmanUtils::Node* dummyParent = new Node(newParent, newInternalNodeId());
    allNodes[dummyParent->identifier] = dummyParent;

    newParent->parent->removeChild(newParent);
    newParent->parent = dummyParent;
    dummyParent->children = {newParent};
    adjustLevels(newParent);

    // newParent now has a 0-length branch from the dummy
    newParent->nucMutation.clear();
    newParent->blockMutation.clear();
    newParent->branchLength = 0;

    // Move node to be a child of the dummy
    toMove->parent->removeChild(toMove);

    toMove->parent = dummyParent;
    dummyParent->children.emplace_back(toMove);
    adjustLevels(toMove);

    // TODO: figure out how branch length works
    toMove->branchLength = 1;
    toMove->nucMutation = pathMutations.nucMutation;
    toMove->blockMutation = pathMutations.blockMutation;
    // TODO: what is that?
    toMove->isComMutHead = false;
}