#include <iostream>
#include <fstream>
#include <string>
#include <google/protobuf/message.h>

#include "PangenomeMAT.hpp"

#include "mutation_annotation_test_proto3_optional_new.pb.h"

// Function to calculate memory usage of a protobuf message
size_t CalculateMessageSize(const MATNew::tree& tree) {
    std::cout << tree.newick() << std::endl;
    return (size_t)0;
    // return message.ByteSizeLong();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_proto_file.pb" << std::endl;
        return 1;
    }

    // Read the protobuf file
    std::fstream input(argv[1], std::ios::in | std::ios::binary);
    if (!input) {
        std::cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    // Create a message object of your Protobuf type
    MATNew::tree tree;
    tree.ParseFromIstream(&input);
    // std::cout << "Newick String Size:" << tree.newick().size() << std::endl;
    // std::cout << "Nodes Size: " << tree.nodes().SpaceUsedExcludingSelf() << std::endl;
    // std::cout << "Consensus Map Size: " << tree.consensusseqmap().SpaceUsedExcludingSelf() << std::endl;
    // std::cout << "Gaps Size: " << tree.gaps().SpaceUsedExcludingSelf() << std::endl;
    // std::cout << "Block Gap Size: " << tree.blockgaps().ByteSizeLong() << std::endl;
    // std::cout << "Circular Offset Size: " << tree.circularsequences().SpaceUsedExcludingSelf() << std::endl;
    

    // Node Size
    int node_size = 0;
    int nucmut_size = 0;
    int blkmut_size = 0;
    for (auto &n: tree.nodes())
    {
        node_size += n.ByteSizeLong();
        for (auto &m: n.mutations())
        { 
            blkmut_size += m.ByteSizeLong();
            for (auto &c: m.nucmutation())
            {
                nucmut_size+=c.ByteSizeLong();
            } 
        }
    }
    std::cout << "Nodes Size (Accurate): " << node_size/(1024*1024) << " MB" << std::endl;
    std::cout << "NucMut Size (Accurate): " << nucmut_size/(1024*1024) << " MB" << std::endl;
    std::cout << "BlockMut Size (Accurate): " << (blkmut_size - nucmut_size)/(1024*1024) << " MB" << std::endl;


    // Consensus Size
    int consensus_size = 0;
    for (auto &n: tree.consensusseqmap())
    {
        consensus_size += n.ByteSizeLong();
    }
    std::cout << "Consensus Map Size (Accurate): " << consensus_size/(1024*1024) << " MB" << std::endl;

    // Gaps Size
    int gap_size = 0;
    for (auto &n: tree.gaps())
    {
        gap_size += n.ByteSizeLong();
    }
    std::cout << "Gap Size (Accurate): " << gap_size/(1024*1024) << " MB" << std::endl;

    // Block Gap Size
    // std::cout << "Block Gap Size: " << tree.blockgaps().ByteSizeLong() << std::endl;



    return 0;
}
