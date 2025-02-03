
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs() {
    std::cout << "Imputing a tree" << std::endl;

    std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > > snvs;
    std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > > insertions;
    findMutationsToN(root, snvs, insertions);

    for (const auto toImpute: snvs) {
        for (const auto& curMut: toImpute.second) {
            imputeSNV(toImpute.first, curMut, nullptr, 0, snvs);
        }
    }
    
    for (const auto& toImpute: insertions) {
        for (const auto& curMut: toImpute.second) {
            imputeInsertion(toImpute.first, curMut, nullptr, 0, insertions);
        }
    }
}

const void panmanUtils::Tree::findMutationsToN(panmanUtils::Node* node, 
        std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > >& snvs,
        std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > >& insertions) {
}

std::string panmanUtils::Tree::imputeSNV(panmanUtils::Node* node, panmanUtils::NucMutPosition mutToN, panmanUtils::Node* childToIgnore, int distanceSoFar,
        const std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > >& mutsToIgnore) {
    return "";
}

std::string panmanUtils::Tree::imputeInsertion(panmanUtils::Node* node, panmanUtils::NucMutPosition mutToN, panmanUtils::Node* childToIgnore, int distanceSoFar,
        const std::unordered_map< panmanUtils::Node*, std::unordered_set< panmanUtils::NucMutPosition > >& mutsToIgnore) {
    return "";
}