
#include "panmanUtils.hpp"

// Build map of all coordinates to their reference/consensus nucleotide
std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > getNucs(std::vector<panmanUtils::Block> blocks) {
    std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > nucs;

    for (const auto& curBlock: blocks) {
        for (int i = 0; i < curBlock.consensusSeq.size(); i++) {
            bool endFlag = false;
            for (int j = 0; j < 8; j++) {
                int curNucCode = ((curBlock.consensusSeq[i] >> (4*(7 - j))) & 15);
                if(curNucCode == 0) {
                    endFlag = true;
                    break;
                } else {
                    nucs[panmanUtils::Coordinate(i*8 + j, 0, curBlock.primaryBlockId, 
                                                 curBlock.secondaryBlockId)] = (int8_t) curNucCode;
                }
            }

            if(endFlag) break;
        }
    }

    return nucs;
}

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< std::string, panmanUtils::NucMut > > substitutions;
    std::unordered_map< std::string, std::vector<panmanUtils::IndelPosition > > insertions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > > allInsertions;
    std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher >  curNucs = getNucs(blocks);
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > > originalNucs;
    
    findMutationsToN(root, substitutions, insertions, allInsertions, curNucs, originalNucs);

    for (const auto& toImpute: substitutions) {
        imputeSNV(allNodes[toImpute.first], toImpute.second);
    }
    std::cout << "Imputed " << substitutions.size() << "/" << substitutions.size() << " SNPs/MNPs to N" << std::endl;
    
    size_t insertionImputeSuccesses = 0;
    for (const auto& toImpute: insertions) {
        insertionImputeSuccesses += imputeInsertion(
            allNodes[toImpute.first], toImpute.second, allowedIndelDistance, allInsertions, originalNucs);
    }
    std::cout << "Imputed " << insertionImputeSuccesses << "/" << insertions.size() << " insertions to N" << std::endl;

    size_t numLeaves;
    size_t totalLeafDepth;
    fixLevels(root, numLeaves, totalLeafDepth);
    m_meanDepth = totalLeafDepth / numLeaves;
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
        std::unordered_map< std::string, std::vector<panmanUtils::IndelPosition> >& insertions,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
        std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher >& curNucs,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return;

    std::string curID = node->identifier;
    insertions[curID] = std::vector<panmanUtils::IndelPosition>();
    allInsertions[curID] = std::unordered_map< panmanUtils::IndelPosition, int32_t >();
    int curInsertionNumNs = 0;

    for (const auto& curMut: node->nucMutation) {
        int numNs = 0;
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            numNs += (curNucCode == panmanUtils::NucCode::N);

            if (curNucs.find(curPos) == curNucs.end()) {
                curNucs[curPos] = panmanUtils::NucCode::MISSING;
            }
            originalNucs[node->identifier][curPos] = curNucs[curPos];
            curNucs[curPos] = curNucCode;
        }

        // Save mutation if relevant
        if (curMut.isSubstitution()) {
            if (numNs > 0) {
                substitutions.emplace_back(node->identifier, curMut);
            }
        } else if (curMut.isInsertion()) {
            if (insertions[curID].empty()) {
                insertions[curID].emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionNumNs = numNs;
            } else if (insertions[curID].back().mergeIndels(curMut)) {
                // Current insertion was merged with the previous one.
                curInsertionNumNs += numNs;
            } else {
                allInsertions[curID].emplace(insertions[curID].back(), curInsertionNumNs);
                if (curInsertionNumNs == 0) {
                    insertions[curID].erase(insertions[curID].end() - 1);
                }
                // Start new insertion
                insertions[curID].emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionNumNs = numNs;
            }
        }
    }
    
    if (insertions[curID].empty()) {
        // No insertions with Ns
        insertions.erase(curID);
    } else {
        // Handle last insertion
        allInsertions[curID].emplace(insertions[curID].back(), curInsertionNumNs);
        if (curInsertionNumNs == 0) {
            insertions[curID].erase(insertions[curID].end() - 1);
            if (insertions[curID].empty()) insertions.erase(curID);
        }
    }

    if (allInsertions[curID].empty()) allInsertions.erase(curID);

    for(auto child: node->children) {
        findMutationsToN(child, substitutions, insertions, allInsertions, curNucs, originalNucs);
    }

    // Undo mutations before passing back up the tree
    for (const auto& curMut: node->nucMutation) {
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            curNucs[curPos] = originalNucs[node->identifier][curPos];
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

bool panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node,
    const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return false;

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    MutationList bestMutationList = MutationList();

    for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance, node, allInsertions, originalNucs)) {
        if (nearby.first != node->identifier) {
            panmanUtils::MutationList simpleMutations = nearby.second.concat(panmanUtils::MutationList(node, false, originalNucs[node->identifier]));
            std::reverse(simpleMutations.nucMutation.begin(), simpleMutations.nucMutation.end());
            simpleMutations.nucMutation = consolidateNucMutations(simpleMutations.nucMutation);
            simpleMutations.blockMutation = consolidateBlockMutations(simpleMutations.blockMutation);
            // Ignore block mutations for now
            if (simpleMutations.blockMutation.size() > 0) continue;

            // Parsimony improvement score is the decrease in mutation count
            int nucImprovement = 0;
            for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();
            for (const auto& curMut: simpleMutations.nucMutation) nucImprovement -= curMut.length();

            if (nucImprovement > bestParsimonyImprovement) {
                bestParsimonyImprovement = nucImprovement;
                bestNewParent = allNodes[nearby.first];
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

const std::unordered_map< std::string, panmanUtils::MutationList > panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {

    std::unordered_map< std::string, panmanUtils::MutationList > nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (allInsertions[curID].find(curMut) != allInsertions[curID].end()) {
            // Only use if this insertion has non-N nucleotides to contribute
            if (allInsertions[curID][curMut] < curMut.length) {
                nearbyInsertions.emplace(node->identifier, panmanUtils::MutationList());
            }
            break;
        }
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore) {
            for (const auto& nearby: findNearbyInsertions(child, mutsToN, allowedDistance - child->branchLength, 
                                                          node, allInsertions, originalNucs)) {
                // Add mutations to get to child (which must be reversed)
                panmanUtils::MutationList toAdd = panmanUtils::MutationList(child, true, originalNucs[child->identifier]);
                nearbyInsertions.emplace(nearby.first, nearby.second.concat(toAdd));
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance - node->branchLength,
                                                      node, allInsertions, originalNucs)) {
            // Add mutations to get to parent
            panmanUtils::MutationList toAdd = panmanUtils::MutationList(node, false, originalNucs[node->identifier]);
            nearbyInsertions.emplace(nearby.first, nearby.second.concat(toAdd));
        }
    }
    return nearbyInsertions;
}

void panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList mutList) {
    // Make dummy parent from grandparent -> dummy -> newParent
    panmanUtils::Node* dummyParent = new Node(newParent, newInternalNodeId());
    allNodes[dummyParent->identifier] = dummyParent;

    newParent->changeParent(dummyParent);
    Node* oldParent = toMove->changeParent(dummyParent);

    // Merge old parent with single remaining child
    if (oldParent->children.size() == 1) {
        mergeNodes(oldParent, oldParent->children[0]);
    }

    // newParent now has a 0-length branch from the dummy
    newParent->nucMutation.clear();
    newParent->branchLength = 0;

    // TODO: figure out how branch length works
    toMove->branchLength = 1;
    toMove->nucMutation = mutList.nucMutation;
}