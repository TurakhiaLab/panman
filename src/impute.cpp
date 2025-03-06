
#include "panmanUtils.hpp"
#include <random>

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< std::string, panmanUtils::NucMut > > substitutions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > > insertions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > originalNucs;
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > > wasBlockInv;

    // Make pre-order pass over the tree, building lookup tables
    fillImputationLookupTables(substitutions, insertions, originalNucs, wasBlockInv);

    // Impute all substitutions (100% success rate)
    int totalSubNs = 0;
    for (const auto& toImpute: substitutions) {
        totalSubNs += imputeSubstitution(allNodes[toImpute.first]->nucMutation, toImpute.second);
    }
    std::cout << "Imputed " << totalSubNs << "/" << totalSubNs << " SNPs/MNPs to N" << std::endl;
    
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
            oldParents.push_back(curNode->parent);

            if (!moveNode(curNode, curMove.second.first, curMove.second.second)) {
                // The move failed
                oldParents.erase(oldParents.end() - 1);
            }
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
    sequence_t curSequence;
    blockExists_t blockExists;
    blockStrand_t blockStrand;
    getSequenceFromReference(curSequence, blockExists, blockStrand, root->identifier);

    // Never used, but needed so that the key is in the map
    insertions[root->identifier] = std::unordered_map< panmanUtils::IndelPosition, int32_t >();
    originalNucs[root->identifier] = std::unordered_map< panmanUtils::Coordinate, int8_t >();
    wasBlockInv[root->identifier] = std::unordered_map< uint64_t, bool >();

    for (const auto& child: root->children) {
        fillImputationLookupTablesHelper(child, substitutions, insertions, originalNucs, wasBlockInv, curSequence, blockStrand);
    }
}

const void panmanUtils::Tree::fillImputationLookupTablesHelper(panmanUtils::Node* node, 
    std::vector< std::pair < std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv,
    sequence_t& curSequence, blockExists_t& blockStrand) {

    if (node == nullptr) return;

    fillNucleotideLookupTables(node, curSequence, substitutions, insertions, originalNucs);
    fillBlockLookupTables(node, blockStrand, wasBlockInv);

    for(auto child: node->children) {
        fillImputationLookupTablesHelper(child, substitutions, insertions, originalNucs, wasBlockInv, curSequence, blockStrand);
    }

    // Undo mutations before passing back up the tree
    for (const auto& curMut: node->nucMutation) {
        for(int i = 0; i < curMut.length(); i++) {
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);
            char originalNuc = panmanUtils::getNucleotideFromCode(originalNucs[node->identifier][curPos]);
            curPos.setSequenceBase(curSequence, originalNuc);
        }
    }

    for (const auto& curMut: node->blockMutation) {
        if (curMut.inversion) { 
            if(curMut.secondaryBlockId != -1) {
                blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId] = 
                    !blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId];
            } else {
                blockStrand[curMut.primaryBlockId].first = !blockStrand[curMut.primaryBlockId].first;
            }
        }
    }
}

const void panmanUtils::Tree::fillNucleotideLookupTables(panmanUtils::Node* node, sequence_t& curSequence,
    std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {
    
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

            // Mutate this position, but store the original value
            originalNucs[node->identifier][curPos] = panmanUtils::getCodeFromNucleotide(curPos.getSequenceBase(curSequence));
            curPos.setSequenceBase(curSequence, panmanUtils::getNucleotideFromCode(curNucCode));
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

    // Store curNodeInsertions as this node's insertions
    std::copy(curNodeInsertions.begin(), curNodeInsertions.end(), 
              std::inserter(insertions[node->identifier], insertions[node->identifier].begin()));
}

const void panmanUtils::Tree::fillBlockLookupTables(panmanUtils::Node* node, blockExists_t& blockStrand,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    
    wasBlockInv[node->identifier] = std::unordered_map< uint64_t, bool >();
    for (const auto& curMut: node->blockMutation) {
        // Store original state for all deletions
        if (curMut.isDeletion()) {
            bool originalState;
            if(curMut.secondaryBlockId != -1) {
                originalState = !blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId];
            } else {
                originalState= !blockStrand[curMut.primaryBlockId].first;
            }
            uint64_t curID = curMut.singleBlockID();
            wasBlockInv[node->identifier][curMut.singleBlockID()] = originalState;
        }
    }
}

const int panmanUtils::Tree::imputeSubstitution(std::vector<panmanUtils::NucMut>& nucMutation, const NucMut& mutToN) {
    // Get rid of the old mutation in the node's list
    std::vector<NucMut>::iterator oldIndex = std::find(nucMutation.begin(), nucMutation.end(), mutToN);
    oldIndex = nucMutation.erase(oldIndex);
    int subNs = mutToN.length();

    // Possible MNP
    if (mutToN.type() == panmanUtils::NucMutationType::NS) {
        std::vector<NucMut> snps;
        // Add non-N mutations back in (for MNPs which are partially N)
        for(int i = 0; i < mutToN.length(); i++) {
            if (mutToN.getNucCode(i) != panmanUtils::NucCode::N) {
                snps.push_back(NucMut(mutToN, i));
            }
        }
        nucMutation.insert(oldIndex, snps.begin(), snps.end());
        // These SNPs were not erased
        subNs -= snps.size();
    }
    return subNs;
}

const void panmanUtils::Tree::imputeAllSubstitutionsWithNs(std::vector<panmanUtils::NucMut>& nucMutation) {
    // Loop over vector backwards as elements will be erased
    for (int i = nucMutation.size() - 1; i >= 0; i--) {
        panmanUtils::NucMut curMut = nucMutation[i];

        if (curMut.isSubstitution()) {
            for (int i = 0; i < curMut.length(); i++) {
                if (curMut.getNucCode(i) == panmanUtils::NucCode::N) {
                    imputeSubstitution(nucMutation, curMut);
                    break;
                }
            }
        }
    }
}

const std::pair< panmanUtils::Node*, panmanUtils::MutationList > panmanUtils::Tree::findInsertionImputationMove(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
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

    std::vector<panmanUtils::IndelPosition> toImpute;
    for (const auto& curInsertion: insertions.at(node->identifier)) {
        if (curInsertion.second > 0) toImpute.push_back(curInsertion.first);
    }
    
    for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance, node,
                                                  insertions, originalNucs, wasBlockInv)) {
        panmanUtils::MutationList curNewMuts = nearby.second.concat(MutationList(node));
        curNewMuts.nucMutation = consolidateNucMutations(curNewMuts.nucMutation);
        imputeAllSubstitutionsWithNs(curNewMuts.nucMutation);
        curNewMuts.blockMutation = consolidateBlockMutations(curNewMuts.blockMutation);

        // Parsimony improvement score is the decrease in mutation count
        int nucImprovement = 0;
        int blockImprovement = node->blockMutation.size() - curNewMuts.blockMutation.size();
        for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();
        for (const auto& curMut: curNewMuts.nucMutation) nucImprovement -= curMut.length();

        if (nucImprovement > bestNucImprovement & blockImprovement >= bestBlockImprovement) {
            bestNucImprovement = nucImprovement;
            bestBlockImprovement = blockImprovement;
            bestNewParent = nearby.first;
            bestNewMuts = curNewMuts;
        }
    }

    return std::make_pair(bestNewParent, bestNewMuts);
}

const std::vector<std::pair< panmanUtils::Node*, panmanUtils::MutationList >> panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {

    std::vector<std::pair< panmanUtils::Node*, panmanUtils::MutationList >> nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (insertions.at(curID).find(curMut) != insertions.at(curID).end()) {
            // Only use if this insertion has non-N nucleotides to contribute
            if (insertions.at(curID).at(curMut) < curMut.length) {
                nearbyInsertions.emplace_back(node, MutationList());
            }
            break;
        }
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore) {
            auto childPossibilities = findNearbyInsertions(
                child, mutsToN, allowedDistance - child->branchLength, 
                node, insertions, originalNucs, wasBlockInv);
            
            // Only bother with getting/inverting mutations if necessary
            if (!childPossibilities.empty()) {
                // Add mutations to get to child (which must be reversed)
                panmanUtils::MutationList toAdd = MutationList(child);
                toAdd.invertMutations(originalNucs.at(child->identifier), wasBlockInv.at(child->identifier));
                for (const auto& nearby: childPossibilities) {
                    nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(toAdd));
                }
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance - node->branchLength,
                                                      node, insertions, originalNucs, wasBlockInv)) {
            // Add mutations to get to parent
            nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(MutationList(node)));
        }
    }
    return nearbyInsertions;
}

bool panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList newMuts) {
    // Prevent looping
    if (newParent->isDescendant(toMove)) return false;

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
    return true;
}

void panmanUtils::Tree::testImputation(int p, int allowedIndelDistance) {
    auto masked = maskNs(p);
    imputeNs(allowedIndelDistance);
    std::tuple< int, int, int, int > subMaskedImputeCounts = checkImputeOfMasked(true, masked["SUB"]);
    std::tuple< int, int, int, int > insMaskedImputeCounts = checkImputeOfMasked(false, masked["INS"]);

    std::cout << "Substitution imputations:" << std::endl;
    std::cout << "\t" << std::get<0>(subMaskedImputeCounts) << " ignored (could not be checked)" << std::endl;
    std::cout << "\t" << std::get<1>(subMaskedImputeCounts) << " unimputed (remained N)" << std::endl;
    std::cout << "\t" << std::get<2>(subMaskedImputeCounts) << " wrong (imputed to different non-N base)" << std::endl;
    std::cout << "\t" << std::get<3>(subMaskedImputeCounts) << " correct (imputed to original base)" << std::endl;

    std::cout << "Insertion imputations:" << std::endl;
    std::cout << "\t" << std::get<0>(insMaskedImputeCounts) << " ignored (could not be checked)" << std::endl;
    std::cout << "\t" << std::get<1>(insMaskedImputeCounts) << " unimputed (remained N)" << std::endl;
    std::cout << "\t" << std::get<2>(insMaskedImputeCounts) << " wrong (imputed to different non-N base)" << std::endl;
    std::cout << "\t" << std::get<3>(insMaskedImputeCounts) << " correct (imputed to original base)" << std::endl;
}

std::unordered_map< std::string, std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > > panmanUtils::Tree::maskNs(int p) {
    // Set up random distribution
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,99);

    // Set up masking tracker
    std::unordered_map< std::string, std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > > masked;
    masked["SUB"] = std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >();
    masked["INS"] = std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >();

    for (const auto& node: allNodes) {
        // don't impute root
        if (node.second != root) {
            masked["SUB"][node.first] = std::unordered_map< panmanUtils::Coordinate, int8_t >();
            masked["INS"][node.first] = std::unordered_map< panmanUtils::Coordinate, int8_t >();

            for (auto& curMut: node.second->nucMutation) {
                if (!curMut.isDeletion()) {
                    for (int i = 0; i < curMut.length(); i++) {
                        if (dist(rng) < p) {
                            // Save original nucleotide
                            if (curMut.isSubstitution()) {
                                masked["SUB"][node.first][panmanUtils::Coordinate(curMut, i)] = curMut.getNucCode(i);
                            } else {
                                masked["INS"][node.first][panmanUtils::Coordinate(curMut, i)] = curMut.getNucCode(i);
                            }
                            // Mask
                            curMut.changeNucCode(panmanUtils::NucCode::N, i);
                        }
                    }
                }
            }
        }
    }

    return masked;
}

std::tuple< int, int, int, int > panmanUtils::Tree::checkImputeOfMasked(bool erasedIsCorrect,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > masked) {
    
    int ignore = 0;
    int unimputed = 0;
    int wrong = 0;
    int correct = 0; 
    for (const auto& curMasked: masked) {
        std::unordered_map< panmanUtils::Coordinate, int8_t > changeList = curMasked.second;
        if (allNodes.find(curMasked.first) != allNodes.end()) {
            for (const auto& curMut: allNodes[curMasked.first]->nucMutation) {
                for (int i = 0; i < curMut.length(); i++) {
                    panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

                    if (changeList.find(curPos) != changeList.end()) {
                        int8_t curCode = curMut.getNucCode(i);

                        // This mutation has an N-masked coordinate
                        if (curCode == panmanUtils::NucCode::N) {
                            unimputed++;
                        } else if (curCode == changeList[curPos]) {
                            correct++;
                        } else {
                            wrong++;
                        }

                        changeList.erase(curPos);
                    }
                }
            }

            if (erasedIsCorrect) {
                correct += changeList.size();
                changeList.clear();
            }
        }

        // Imputation that couldn't be checked for any reason
        ignore += changeList.size();
    }

    return std::make_tuple(ignore, unimputed, wrong, correct);
}