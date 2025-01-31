
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeViaMaskedFitch() {
    std::cout << "Imputing a tree" << std::endl;

    // Identify all coordinates which have a nucleotide substitution to N
    // Also get a list of all leaves (is that already somewhere?)
    // For each coordinate, loop over all leaves
        // If the leaf has an N, mask it, otherwise save its state
        // Mask internal nodes which have no unmasked children
        // Run Fitch forward pass
        // Unmask all nodes (update "states" as necessary)
        // Run Fitch backward pass
        // Assign mutations as makes sense

    // Identify all nodes which have an insertion of Ns
    // For each node/insertion pair, copy nearby insertions of same length/pos
}