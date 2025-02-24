
#include "panmanUtils.hpp"

// Concatenate two vectors of NucMuts
std::vector<panmanUtils::NucMut> concat(const std::vector<panmanUtils::NucMut>& first, const std::vector<panmanUtils::NucMut>& second) {
    std::vector<panmanUtils::NucMut> temp(first.size() + second.size());
    temp.insert(temp.end(), first.begin(), first.end());
    temp.insert(temp.end(), second.begin(), second.end());
    return temp;
}

// Build map of all coordinates to their reference/consensus nucleotide
std::unordered_map< panmanUtils::Coordinate, int8_t > getNucs(std::vector<panmanUtils::Block> blocks) {
    std::unordered_map< panmanUtils::Coordinate, int8_t > nucs;

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
    std::unordered_map< panmanUtils::Coordinate, int8_t >  curNucs = getNucs(blocks);
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > originalNucs;
    
    // Make one pre-order pass over the tree, building lookup tables
    fillImputationLookupTables(root, substitutions, insertions, allInsertions, curNucs, originalNucs);

    // Impute all substitutions (100% success rate)
    for (const auto& toImpute: substitutions) {
        imputeSubstitution(allNodes[toImpute.first], toImpute.second);
    }
    std::cout << "Imputed " << substitutions.size() << "/" << substitutions.size() << " SNPs/MNPs to N" << std::endl;
    
    // Attempt to impute insertions
    // Find possible places to move nodes to for insertion imputation
    std::unordered_map< std::string, std::pair< panmanUtils::Node*, std::vector<panmanUtils::NucMut> > > toMove;
    for (const auto& toImpute: insertions) {
        toMove[toImpute.first] = findInsertionImputationMove(
            allNodes[toImpute.first], toImpute.second, allowedIndelDistance, allInsertions, originalNucs);
    }
    // Make all moves
    std::vector<panmanUtils::Node*> oldParents;
    for (const auto& curMove: toMove) {
        if (curMove.second.first != nullptr) {
            Node* curNode = allNodes[curMove.first];
            if (curNode->parent != nullptr) oldParents.push_back(curNode->parent);

            std::string dummyParentID = moveNode(curNode, curMove.second.first, curMove.second.second);

            // Transfer insertions from curMove.second.first to new dummy parent
            std::string newSibling = curMove.second.first->identifier;
            insertions[dummyParentID] = insertions[newSibling];
            allInsertions[dummyParentID] = allInsertions[newSibling];
            insertions[newSibling].clear();
            allInsertions[newSibling].clear();

            // Re-calculate insertions for the moved node
            insertions.erase(curMove.first);
            allInsertions.erase(curMove.first);
            fillImputationLookupTablesHelper(curNode, substitutions, insertions, allInsertions, curNucs, originalNucs);
        }
    }
    std::cout << "Moved " << oldParents.size() << "/" << insertions.size() << " nodes with insertions to N" << std::endl;
    int totalInsertions = 0;
    int descImputations = 0;
    // Attempt to impute from descendants
    for (const auto& toImpute: insertions) {
        for (const auto& curMut: toImpute.second) {
            totalInsertions++;
            descImputations += imputeFromDescendant(allNodes[toImpute.first], curMut, allowedIndelDistance, allInsertions);
        }
    }
    std::cout << "Imputed " << descImputations << " /" << totalInsertions << " insertions from descendants" << std::endl;
    // Compress parents with single children left over from moves
    for (const auto& curParent: oldParents) {
        if (curParent->children.size() == 1) {
            std::string parentID = curParent->identifier;
            std::string childID = curParent->children[0]->identifier;
            mergeNodes(curParent, curParent->children[0]);
            // Combine child & parent insertions
            insertions[childID].insert(insertions[childID].end(), insertions[parentID].begin(), insertions[parentID].end());
            allInsertions[childID].insert(allInsertions[parentID].begin(), allInsertions[parentID].end());
        }
    }

    // Fix depth/level attributes, post-move
    size_t numLeaves;
    size_t totalLeafDepth;
    fixLevels(root, numLeaves, totalLeafDepth);
    m_meanDepth = totalLeafDepth / numLeaves;
}

const void panmanUtils::Tree::fillImputationLookupTables(panmanUtils::Node* node, 
        std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
        std::unordered_map< std::string, std::vector<panmanUtils::IndelPosition> >& insertions,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
        std::unordered_map< panmanUtils::Coordinate, int8_t >& curNucs,
        std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {
    if (node == nullptr) return;

    fillImputationLookupTablesHelper(node, substitutions, insertions, allInsertions, curNucs, originalNucs);

    for(auto child: node->children) {
        fillImputationLookupTables(child, substitutions, insertions, allInsertions, curNucs, originalNucs);
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

const void panmanUtils::Tree::fillImputationLookupTablesHelper(panmanUtils::Node* node, 
    std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::vector<panmanUtils::IndelPosition> >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    std::unordered_map< panmanUtils::Coordinate, int8_t >& curNucs,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {
    
    std::string curID = node->identifier;
    std::vector< std::pair< panmanUtils::IndelPosition, int32_t > > curNodeInsertions;
    // Will store the parent's nucleotide at all positions with an insertion, to allow reversability
    originalNucs[node->identifier] = std::unordered_map< panmanUtils::Coordinate, int8_t >();

    for (const auto& curMut: node->nucMutation) {
        int numNs = 0;
        // Handle each base in the current mutation
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            numNs += (curNucCode == panmanUtils::NucCode::N);

            // If this position was previous unknown, set its prior state to missing
            if (curNucs.find(curPos) == curNucs.end()) {
                curNucs[curPos] = panmanUtils::NucCode::MISSING;
            }
            // Mutate this position, but store the original value
            originalNucs[node->identifier][curPos] = curNucs[curPos];
            curNucs[curPos] = curNucCode;
        }

        // Save mutation if relevant
        if (curMut.isSubstitution()) {
            if (numNs > 0) {
                substitutions.emplace_back(node->identifier, curMut);
            }
        } else if (curMut.isInsertion()) {
            if (!curNodeInsertions.empty() && curNodeInsertions.back().first.mergeIndels(curMut)) {
                // Current insertion was merged with the previous one.
                curNodeInsertions.back().second += numNs;
            } else {
                // Start new insertion
                curNodeInsertions.emplace_back(panmanUtils::IndelPosition(curMut), numNs);
            }
        }
    }

    // Transfer curNodeInsertions into insertions and allInsertions
    // TODO: now that allInsertions stores numNs, could get away with storing only it
    insertions[curID] = std::vector<panmanUtils::IndelPosition>();
    allInsertions[curID] = std::unordered_map< panmanUtils::IndelPosition, int32_t >();

    for (const auto& curInsert: curNodeInsertions) {
        allInsertions[curID].emplace(curInsert.first, curInsert.second);
        if (curInsert.second > 0) insertions[curID].push_back(curInsert.first);
    }

    // If no insertions have Ns, don't store this node as needing insertions
    if (insertions[curID].empty()) insertions.erase(curID);
}

const void panmanUtils::Tree::imputeSubstitution(panmanUtils::Node* node, NucMut mutToN) {
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

const std::pair< panmanUtils::Node*, std::vector<panmanUtils::NucMut> > panmanUtils::Tree::findInsertionImputationMove(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {
    // Certain cases are simply impossible
    if (node == nullptr || !node->blockMutation.empty()) {
        return std::make_pair(nullptr, std::vector<panmanUtils::NucMut>());
    }

    // Tracking best new position so far
    int bestParsimonyImprovement = -1;
    Node* bestNewParent = nullptr;
    std::vector<panmanUtils::NucMut> bestNewMuts;

    for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance, node, allInsertions, originalNucs)) {
        if (nearby.first != node->identifier) {
            std::vector<panmanUtils::NucMut> curNewMuts = concat(nearby.second, node->nucMutation);
            std::reverse(curNewMuts.begin(), curNewMuts.end());
            curNewMuts = consolidateNucMutations(curNewMuts);

            // Parsimony improvement score is the decrease in mutation count
            int nucImprovement = 0;
            for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();
            for (const auto& curMut: curNewMuts) nucImprovement -= curMut.length();

            if (nucImprovement > bestParsimonyImprovement) {
                bestParsimonyImprovement = nucImprovement;
                bestNewParent = allNodes[nearby.first];
                bestNewMuts = curNewMuts;
            }
        }
    }

    return std::make_pair(bestNewParent, bestNewMuts);
}

const std::unordered_map< std::string, std::vector<panmanUtils::NucMut> > panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {

    std::unordered_map< std::string, std::vector<panmanUtils::NucMut> > nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (allInsertions.at(curID).find(curMut) != allInsertions.at(curID).end()) {
            // Only use if this insertion has non-N nucleotides to contribute
            if (allInsertions.at(curID).at(curMut) < curMut.length) {
                nearbyInsertions.emplace(node->identifier, std::vector<panmanUtils::NucMut>());
            }
            break;
        }
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore && child->blockMutation.empty()) {
            // Add mutations to get to child (which must be reversed)
            std::vector<panmanUtils::NucMut> toAdd = child->nucMutation;
            reverseNucMutations(toAdd, originalNucs.at(child->identifier));
            for (const auto& nearby: findNearbyInsertions(child, mutsToN, allowedDistance - child->branchLength, 
                                                          node, allInsertions, originalNucs)) {
                nearbyInsertions[nearby.first] = concat(nearby.second, toAdd);
            }
        }
    }
    // Try parent
    if (node->parent != ignore && node->blockMutation.empty()) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance - node->branchLength,
                                                      node, allInsertions, originalNucs)) {
            // Add mutations to get to parent
            nearbyInsertions[nearby.first] = concat(nearby.second, node->nucMutation);
        }
    }
    return nearbyInsertions;
}

std::string panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, std::vector<panmanUtils::NucMut> newMuts) {
    panmanUtils::Node* oldParent = toMove->parent;
    // Make dummy parent from grandparent -> dummy -> newParent
    panmanUtils::Node* dummyParent = new Node(newParent, newInternalNodeId());
    allNodes[dummyParent->identifier] = dummyParent;

    newParent->changeParent(dummyParent);
    toMove->changeParent(dummyParent);

    // newParent now has a 0-length branch from the dummy
    newParent->nucMutation.clear();
    newParent->branchLength = 0;

    // TODO: figure out how branch length works
    toMove->branchLength = 1;
    toMove->nucMutation = newMuts;

    return dummyParent->identifier;
}



const bool panmanUtils::Tree::imputeFromDescendant(panmanUtils::Node* node, 
    const panmanUtils::IndelPosition& mutToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions) {
    if (node == nullptr) return false;
    
    std::vector<panmanUtils::Node*> bestChildPath;
    int32_t bestChildLength;
    for (const auto& curPath: findChildInsertions(node, mutToN, allowedDistance, allInsertions)) {
        if (bestChildPath.size() > 1 || curPath.second < bestChildLength) {
            bestChildPath = curPath.first;
            bestChildLength = curPath.second;
        }
    }

    // Did imputing from this child succeed?
    // TODO: allow imputing from further children if the closests fails?
    bool success = false;
    if (!bestChildPath.empty()) {
        if (bestChildPath.size() >= 5) {
            panmanUtils::Node* child = bestChildPath[0];
            // Find where the parent/child insertions are
            int parentI, childI;
            for (int i = 0; i < node->nucMutation.size(); i++) {
                if (panmanUtils::Coordinate(node->nucMutation[i]) == mutToN.pos) {
                    parentI = i;
                    break;
                }
            }
            for (int i = 0; i < child->nucMutation.size(); i++) {
                if (panmanUtils::Coordinate(child->nucMutation[i]) == mutToN.pos) {
                    childI = i;
                    break;
                }
            }
            // Copy child values up to parent
            int parentJ = 0;
            int childJ = 0;
            int8_t parentNuc, childNuc;
            for (int i = 0; i < mutToN.length; i++) {
                parentNuc = node->nucMutation[parentI].getNucCode(parentJ);
                childNuc = child->nucMutation[childI].getNucCode(childJ);
                if (parentNuc == panmanUtils::NucCode::N && childNuc != panmanUtils::NucCode::N) {
                    node->nucMutation[parentI].changeNucCode(childNuc, parentJ);
                    success = true;
                }
                // Update pointers to move to the next inserted nucleotide
                parentJ++;
                if (parentJ == node->nucMutation[parentI].length()) {
                    parentI++;
                    parentJ = 0;
                }
                childJ++;
                if (childJ == child->nucMutation[childI].length()) {
                    childI++;
                    childJ = 0;
                }
            }
        } else {
            std::cout << "path from " << bestChildPath[0]->identifier << " to " << node->identifier << ": ";
        }
    }
    return success;
}

const std::vector< std::pair< std::vector<panmanUtils::Node*>, int32_t > > panmanUtils::Tree::findChildInsertions(
    panmanUtils::Node* node, const panmanUtils::IndelPosition& mutToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions) {
    
    std::vector< std::pair< std::vector<panmanUtils::Node*>, int32_t > > childInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return childInsertions;

    std::string curID = node->identifier;
    if (allInsertions.at(curID).find(mutToN) != allInsertions.at(curID).end()) {
        // Only use if this insertion has non-N nucleotides to contribute
        if (allInsertions.at(curID).at(mutToN) < mutToN.length) {
            std::vector<panmanUtils::Node*> path = {node};
            childInsertions.emplace_back(path, node->branchLength);
        }
    }

    for (const auto& child: node->children) {
        for (const auto& nearby: findChildInsertions(child, mutToN, allowedDistance - child->branchLength, allInsertions)) {
            std::vector<panmanUtils::Node*> path = nearby.first;
            path.push_back(node);
            childInsertions.emplace_back(path, nearby.second + node->branchLength);
        }
    }

    return childInsertions;
}