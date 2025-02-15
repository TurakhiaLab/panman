
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
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> > allInsertions;
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
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
        std::unordered_map< std::string, std::vector<panmanUtils::IndelPosition> >& insertions,
        std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
        std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher >& curNucs,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return;

    std::string curID = node->identifier;
    insertions[curID] = std::vector<panmanUtils::IndelPosition>();
    allInsertions[curID] = std::unordered_set<panmanUtils::IndelPosition>();
    bool curInsertionHasNs = false;

    for (const auto& curMut: node->nucMutation) {
        // Does this mutation have Ns?
        bool hasNs = false;
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            hasNs |= (curNucCode == panmanUtils::NucCode::N);

            if (curNucs.find(curPos) == curNucs.end()) {
                curNucs[curPos] = panmanUtils::NucCode::MISSING;
            }
            originalNucs[node->identifier][curPos] = curNucs[curPos];
            curNucs[curPos] = curNucCode;
        }

        // Save mutation if relevant
        if (curMut.isSubstitution()) {
            if (hasNs) {
                substitutions.emplace_back(node->identifier, curMut);
            }
        } else if (curMut.isInsertion()) {
            if (insertions[curID].empty()) {
                insertions[curID].emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionHasNs = hasNs;
            } else if (insertions[curID].back().mergeIndels(curMut)) {
                // Current insertion was merged with the previous one.
                curInsertionHasNs |= hasNs;
            } else {
                allInsertions[curID].emplace(insertions[curID].back());
                if (!curInsertionHasNs) {
                    insertions[curID].erase(insertions[curID].end() - 1);
                }
                // Start new insertion
                insertions[curID].emplace_back(panmanUtils::IndelPosition(curMut));
                curInsertionHasNs = hasNs;
            }
        }
    }
    
    if (insertions[curID].empty()) {
        // No insertions with Ns
        insertions.erase(curID);
    } else {
        // Handle last insertion
        allInsertions[curID].emplace(insertions[curID].back());
        if (!curInsertionHasNs) {
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
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {
    if (node == nullptr) return false;
    std::cout << "imputing " << node->identifier << " insertions: ";
    for (const auto& curMut: mutsToN) {
        std::cout << "length " << curMut.length << " (" << curMut.pos.primaryBlockId << "," << curMut.pos.secondaryBlockId;
        std::cout << "," << curMut.pos.nucPosition << "," << curMut.pos.nucGapPosition << ")";
    }
    std::cout << std::endl;

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    MutationList bestMutationList = MutationList();

    for (const auto& nearby: findNearbyInsertions(node, mutsToN, allowedDistance, nullptr, allInsertions, originalNucs)) {
        if (nearby.first != node->identifier) {
            panmanUtils::MutationList simpleMutations = nearby.second.copy();
            std::reverse(simpleMutations.nucMutation.begin(), simpleMutations.nucMutation.end());
            simpleMutations.nucMutation = consolidateNucMutations(simpleMutations.nucMutation);
            simpleMutations.blockMutation = consolidateBlockMutations(simpleMutations.blockMutation);
            // Ignore block mutations for now
            if (simpleMutations.blockMutation.size() > 0) continue;

            // Parsimony improvement score is the decrease in mutation count
            int nucImprovement = 0;
            int numNs = 0;
            for (const auto& curMut: node->nucMutation) {
                nucImprovement += curMut.length();
                for (int i = 0; i < curMut.length(); i++) {
                    if (curMut.getNucCode(i) == panmanUtils::NucCode::N) numNs++;
                }
            }
            for (const auto& curMut: simpleMutations.nucMutation) {
                nucImprovement -= curMut.length();
                for (int i = 0; i < curMut.length(); i++) {
                    if (curMut.getNucCode(i) == panmanUtils::NucCode::N) numNs--;
                }
            }

            if (nucImprovement > bestParsimonyImprovement && numNs > 0) {
                bestParsimonyImprovement = nucImprovement;
                bestNewParent = allNodes[nearby.first];
                bestMutationList = simpleMutations;
            } else if (nucImprovement < 0) {
                std::cout << nearby.first << " did not improve parsimony" << std::endl;
            }
        }
    }

    if (bestNewParent != nullptr) {
        std::cout << node->identifier << ": "; 
        for (const auto& curMut: node->nucMutation) {
            std::cout << curMut.type() << ",";
            for (int i = 0; i < curMut.length(); i++) {
                std::cout << panmanUtils::getNucleotideFromCode(curMut.getNucCode(i));
            }
            std::cout << " ";
        }
        std::cout << std::endl;
        
        std::cout << bestNewParent->identifier << ": "; 
        for (const auto& curMut: bestNewParent->nucMutation) {
            std::cout << curMut.type() << ",";
            for (int i = 0; i < curMut.length(); i++) {
                std::cout << panmanUtils::getNucleotideFromCode(curMut.getNucCode(i));
            }
            std::cout << " ";
        }
        std::cout << std::endl;
        moveNode(node, bestNewParent, bestMutationList);
        return true;
    } else {
        return false;
    }
}

const std::unordered_map< std::string, panmanUtils::MutationList > panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    std::unordered_map< std::string, std::unordered_set<panmanUtils::IndelPosition> >& allInsertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t, panmanUtils::CoordinateHasher > >& originalNucs) {

    std::unordered_map< std::string, panmanUtils::MutationList > nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (allInsertions[curID].find(curMut) != allInsertions[curID].end()) {
            nearbyInsertions.emplace(node->identifier, panmanUtils::MutationList());
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
    std::cout << "moving " << toMove->identifier << " under " << newParent->identifier << std::endl;

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
    toMove->nucMutation = mutList.nucMutation;
    toMove->blockMutation = mutList.blockMutation;
    // TODO: what is that?
    toMove->isComMutHead = false;
}