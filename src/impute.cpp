
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
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > > insertions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > originalNucs;
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > > wasBlockInv;
    
    // Make pre-order pass over the tree, building lookup tables
    fillImputationLookupTables(substitutions, insertions, originalNucs, wasBlockInv);

    // Impute all substitutions (100% success rate)
    for (const auto& toImpute: substitutions) {
        imputeSubstitution(allNodes[toImpute.first], toImpute.second);
    }
    std::cout << "Imputed " << substitutions.size() << "/" << substitutions.size() << " SNPs/MNPs to N" << std::endl;
    
    // Attempt to impute insertions

    // Store {node ID to move : {node to move under, new mutations}} for all imputation attempts
    std::unordered_map< std::string, std::pair< panmanUtils::Node*, panmanUtils::MutationList > > toMove;
    // Must filter "insertions" to just those with Ns
    std::vector<panmanUtils::IndelPosition> insertionsWithNs;
    int insertionImputationAttempts = 0;

    // Find possible places to move nodes to for insertion imputation
    for (const auto& toImpute: insertions) {
        insertionsWithNs.clear();
        for (const auto& curInsertion: toImpute.second) {
            if (curInsertion.second > 0) {
                insertionsWithNs.push_back(curInsertion.first);
            }
        }

        // Only attempt an imputation if necessary
        if (!insertionsWithNs.empty()) {
            insertionImputationAttempts++;
            toMove[toImpute.first] = findInsertionImputationMove(
                allNodes[toImpute.first], insertionsWithNs, allowedIndelDistance, insertions, originalNucs, wasBlockInv);
        }
    }

    // Make all moves
    std::vector<panmanUtils::Node*> oldParents;
    for (const auto& curMove: toMove) {
        if (curMove.second.first != nullptr) {
            Node* curNode = allNodes[curMove.first];
            if (curNode->parent != nullptr) oldParents.push_back(curNode->parent);

            moveNode(curNode, curMove.second.first, curMove.second.second);
        }
    }

    // Compress parents with single children left over from moves
    for (const auto& curParent: oldParents) {
        if (curParent->children.size() == 1) {
            mergeNodes(curParent, curParent->children[0]);
        }
    }

    std::cout << "Moved " << oldParents.size() << "/" << insertionImputationAttempts << " nodes with insertions to N" << std::endl;

    // Fix depth/level attributes, post-move
    size_t numLeaves;
    size_t totalLeafDepth;
    fixLevels(root, numLeaves, totalLeafDepth);
    m_meanDepth = totalLeafDepth / numLeaves;
}

const void panmanUtils::Tree::fillImputationLookupTables( 
    std::vector< std::pair < std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    
    // Prepare current-state trackers
    std::unordered_map< panmanUtils::Coordinate, int8_t > curNucs = getNucs(blocks);
     std::unordered_map< uint64_t, bool > isInv;
    for (const auto& curBlock: blocks) {
        isInv[curBlock.singleBlockID()] = false;
    }

    fillImputationLookupTablesHelper(root, substitutions, insertions, originalNucs, curNucs, wasBlockInv, isInv);
}

const void panmanUtils::Tree::fillImputationLookupTablesHelper(panmanUtils::Node* node, 
    std::vector< std::pair < std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map<panmanUtils::Coordinate, int8_t >& curNucs,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv,
    std::unordered_map< uint64_t, bool >& isInv) {

    if (node == nullptr) return;

    fillNucleotideLookupTables(node, substitutions, insertions, originalNucs, curNucs);
    fillBlockLookupTables(node, wasBlockInv, isInv);

    for(auto child: node->children) {
        fillImputationLookupTablesHelper(child, substitutions, insertions, originalNucs, curNucs, wasBlockInv, isInv);
    }

    // Undo mutations before passing back up the tree
    for (const auto& curMut: node->nucMutation) {
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            curNucs[curPos] = originalNucs[node->identifier][curPos];
        }
    }

    for (const auto& curMut: node->blockMutation) {
        if (curMut.inversion) { 
            uint64_t blockID = curMut.singleBlockID();
            isInv[blockID] = !isInv[blockID];
        }
    }
}

const void panmanUtils::Tree::fillNucleotideLookupTables(panmanUtils::Node* node, 
    std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map< panmanUtils::Coordinate, int8_t >& curNucs) {
    
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
    std::copy(curNodeInsertions.begin(), curNodeInsertions.end(), std::inserter(insertions[curID], insertions[curID].begin()));
}

const void panmanUtils::Tree::fillBlockLookupTables(panmanUtils::Node* node,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv,
    std::unordered_map< uint64_t, bool >& isInv) {
    
    wasBlockInv[node->identifier] = std::unordered_map< uint64_t, bool >();
    for (const auto& curMut: node->blockMutation) {
        // Store original state for all deletions
        if (curMut.isDeletion()) {
            uint64_t curID = curMut.singleBlockID();
            wasBlockInv[node->identifier][curID] = isInv[curID];
        }
    }
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

const std::pair< panmanUtils::Node*, panmanUtils::MutationList > panmanUtils::Tree::findInsertionImputationMove(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    // Certain cases are simply impossible
    if (node == nullptr) {
        return std::make_pair(nullptr, panmanUtils::MutationList());
    }

    // Tracking best new position so far
    int bestNucImprovement = -1;
    int bestBlockImprovement = 0;
    Node* bestNewParent = nullptr;
    panmanUtils::MutationList bestNewMuts;

    for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance, node, allInsertions, originalNucs, wasBlockInv)) {
        if (nearby.first != node->identifier) {
            panmanUtils::MutationList curNewMuts = nearby.second.concat(MutationList(node));
            curNewMuts.reverseMutations();
            curNewMuts.nucMutation = consolidateNucMutations(curNewMuts.nucMutation);
            curNewMuts.blockMutation = consolidateBlockMutations(curNewMuts.blockMutation);

            // Parsimony improvement score is the decrease in mutation count
            int nucImprovement = 0;
            int blockImprovement = node->blockMutation.size() - curNewMuts.blockMutation.size();
            for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();
            for (const auto& curMut: curNewMuts.nucMutation) nucImprovement -= curMut.length();

            if (nucImprovement > bestNucImprovement & blockImprovement >= bestBlockImprovement) {
                bestNucImprovement = nucImprovement;
                bestBlockImprovement = blockImprovement;
                bestNewParent = allNodes[nearby.first];
                bestNewMuts = curNewMuts;
            }
        }
    }

    return std::make_pair(bestNewParent, bestNewMuts);
}

const std::unordered_map< std::string, panmanUtils::MutationList > panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& allInsertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {

    std::unordered_map< std::string, panmanUtils::MutationList > nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (allInsertions.at(curID).find(curMut) != allInsertions.at(curID).end()) {
            // Only use if this insertion has non-N nucleotides to contribute
            if (allInsertions.at(curID).at(curMut) < curMut.length) {
                nearbyInsertions.emplace(node->identifier, MutationList());
            }
            break;
        }
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore) {
            // Add mutations to get to child (which must be reversed)
            panmanUtils::MutationList toAdd = MutationList(child);
            toAdd.invertMutations(originalNucs.at(child->identifier), wasBlockInv.at(child->identifier));
            for (const auto& nearby: findNearbyInsertions(child, mutsToN, allowedDistance - child->branchLength, 
                                                          node, allInsertions, originalNucs, wasBlockInv)) {
                nearbyInsertions[nearby.first] = nearby.second.concat(toAdd);
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance - node->branchLength,
                                                      node, allInsertions, originalNucs, wasBlockInv)) {
            // Add mutations to get to parent
            nearbyInsertions[nearby.first] = nearby.second.concat(MutationList(node));
        }
    }
    return nearbyInsertions;
}

void panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList newMuts) {
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
    toMove->nucMutation = newMuts.nucMutation;
    toMove->blockMutation = newMuts.blockMutation;
}