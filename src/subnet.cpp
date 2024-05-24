#include "panmanUtils.hpp"

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

panmanUtils::Node* panmanUtils::Tree::subtreeExtractParallel(std::vector< std::string > nodeIds) {
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

    compressTreeParallel(newTreeRoot, 1);

    return newTreeRoot;
}
