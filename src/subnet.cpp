#include "panmanUtils.hpp"

void panmanUtils::Tree::compressTreeParallel(panmanUtils::Node* node, size_t level, const std::set< std::string >& nodeIdsToDefinitelyInclude) {
    node->level = level;

    while(node->children.size() == 1) {
        if(nodeIdsToDefinitelyInclude.find(node->identifier) != nodeIdsToDefinitelyInclude.end() ||
                nodeIdsToDefinitelyInclude.find(node->children[0]->identifier)
                != nodeIdsToDefinitelyInclude.end()) {
            break;
        }
        mergeNodes(node, node->children[0]);
        auto oldVector = node->nucMutation;
        node->nucMutation = consolidateNucMutations(oldVector);
        if(!debugSimilarity(oldVector, node->nucMutation)) {
            printError("Inaccuracy observed in subtree extract. Please report to creators");
            return;
        }
    }

    if(node->children.size() == 0) {
        return;
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, (int)node->children.size()), [&](tbb::blocked_range<int> r) {
        for(int i = r.begin(); i < r.end(); i++) {
            while(node->children[i]->children.size() == 1) {
                if(nodeIdsToDefinitelyInclude.find(node->children[i]->identifier) != nodeIdsToDefinitelyInclude.end() ||
                        nodeIdsToDefinitelyInclude.find(node->children[i]->children[0]->identifier) != nodeIdsToDefinitelyInclude.end()) {
                    break;
                }
                mergeNodes(node->children[i], node->children[i]->children[0]);
            }
            // consolidate mutations of parent
            auto oldVector = node->children[i]->nucMutation;

            if(nodeIdsToDefinitelyInclude.find(node->children[i]->identifier) != nodeIdsToDefinitelyInclude.end() ||
                    (node->children[i]->children.size() == 1 &&
                     nodeIdsToDefinitelyInclude.find(node->children[i]->children[0]->identifier) != nodeIdsToDefinitelyInclude.end())) {
                node->children[i]->nucMutation = consolidateNucMutations(node->children[i]->nucMutation);
            }

            if(!debugSimilarity(oldVector, node->children[i]->nucMutation)) {
                printError("Inaccuracy observed in subtree extract. Please report to the"
                           "creators.");
                return;
            }

            compressTreeParallel(node->children[i], level + 1, nodeIdsToDefinitelyInclude);
        }
    });

}

panmanUtils::Node* subtreeExtractParallelHelper(panmanUtils::Node* node, const tbb::concurrent_unordered_map< panmanUtils::Node*, size_t >& ticks) {
    if(ticks.find(node) == ticks.end()) {
        return nullptr;
    }

    panmanUtils::Node* newNode = new panmanUtils::Node(node->identifier, node->branchLength);

    for(auto mutation: node->nucMutation) {
        newNode->nucMutation.push_back(mutation);
    }

    for(auto mutation: node->blockMutation) {
        newNode->blockMutation.push_back(mutation);
    }

    newNode->children.resize(node->children.size(), nullptr);

    // Extract children that have been ticked
    tbb::parallel_for(tbb::blocked_range(0, (int)node->children.size()), [&](tbb::blocked_range<int> r) {
        for(int i = r.begin(); i < r.end(); i++) {
            panmanUtils::Node* child = node->children[i];
            if(ticks.find(child) != ticks.end()) {

                panmanUtils::Node* newChild = subtreeExtractParallelHelper(child, ticks);

                newChild->parent = newNode;
                newNode->children[i] = newChild;
            }
        }
    });

    // Bring all children to front of array
    size_t i = 0, j = 0;
    while(j < newNode->children.size()) {
        if(newNode->children[j] != nullptr) {
            std::swap(newNode->children[i], newNode->children[j]);
            i++;
        }
        j++;
    }
    newNode->children.resize(i);

    return newNode;

}

panmanUtils::Node* panmanUtils::Tree::subtreeExtractParallel(std::vector< std::string > nodeIds, const std::set< std::string >& nodeIdsToDefinitelyInclude) {
    tbb::concurrent_vector< panmanUtils::Node* > requiredNodes;

    std::atomic<bool> idDoesntExist = false;

    tbb::parallel_for_each(nodeIds.begin(), nodeIds.end(), [&](std::string& id) {
        if(allNodes.find(id) != allNodes.end()) {
            requiredNodes.push_back(allNodes[id]);
        } else {
            idDoesntExist = true;
        }
    });

    if(idDoesntExist) {
        printError("Some of the specified node identifiers don't exist!!!");
        return nullptr;
    }

    tbb::concurrent_unordered_map< panmanUtils::Node*, size_t > ticks;

    tbb::parallel_for_each(requiredNodes.begin(), requiredNodes.end(), [&](panmanUtils::Node*& node) {
        Node* current = node;

        while(current != nullptr) {
            ticks[current]++;
            current = current->parent;
        }
    });

    panmanUtils::Node* newTreeRoot = subtreeExtractParallelHelper(root, ticks);

    compressTreeParallel(newTreeRoot, 1, nodeIdsToDefinitelyInclude);

    return newTreeRoot;
}


panmanUtils::TreeGroup* panmanUtils::TreeGroup::subnetworkExtract(std::unordered_map< int, std::vector< std::string > >& nodeIds) {
    // Temporary directory to store subtrees
    std::filesystem::create_directory("./temp_panman");
    std::vector< std::string > fileNames;
    for (size_t i = 0; i < trees.size(); i++) {
        std::vector< std::string > subtreeNodeIds = nodeIds[i];
        std::set< std::string > cplxMutationNodeIds;
        for (auto mutation: complexMutations) {
            if (mutation.treeIndex1 == i) {
                cplxMutationNodeIds.insert(mutation.sequenceId1);
            } else if(mutation.treeIndex2 == i) {
                cplxMutationNodeIds.insert(mutation.sequenceId2);
            } else if(mutation.treeIndex3 == i) {
                cplxMutationNodeIds.insert(mutation.sequenceId3);
            }
        }

        std::set< std::string > subtreeNodeIdSet(subtreeNodeIds.begin(), subtreeNodeIds.end());
        for (auto id: cplxMutationNodeIds) {
            subtreeNodeIdSet.insert(id);
        }

        subtreeNodeIds = std::vector< std::string >(subtreeNodeIdSet.begin(), subtreeNodeIdSet.end());

        std::string fileName = "./temp_panman/panmat_" + std::to_string(i) + ".panmat";

        std::ofstream outputFile(fileName);
        fileNames.push_back(fileName);

        boost::iostreams::filtering_streambuf< boost::iostreams::output>
        outPMATBuffer;
        outPMATBuffer.push(boost::iostreams::gzip_compressor());
        outPMATBuffer.push(outputFile);
        std::ostream outstream(&outPMATBuffer);

        trees[i].writeToFile(outstream, trees[i].subtreeExtractParallel(subtreeNodeIds, cplxMutationNodeIds));

        boost::iostreams::close(outPMATBuffer);
        outputFile.close();
    }

    std::ofstream outputFile;
    std::streambuf * buf;
    outputFile.open("./temp_panman/cplx");
    buf = outputFile.rdbuf();

    std::ostream fout (buf);

    printComplexMutations(fout);

    outputFile.close();

    std::ifstream cplxMutationFile("./temp_panman/cplx");

    std::vector< std::ifstream > files;
    for(auto f: fileNames) {
        files.emplace_back(f);
    }

    TreeGroup* subnetwork = new TreeGroup(files, cplxMutationFile);

    cplxMutationFile.close();

    // Delete temporary directory
    std::filesystem::remove_all("./temp_panman");

    return subnetwork;
}
