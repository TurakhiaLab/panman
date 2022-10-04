#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mutation_annotation.pb.h"

int main(){

	GOOGLE_PROTOBUF_VERIFY_VERSION;

	std::ifstream input("sample_mat.pb");

	MAT::tree mainTree;

	if(!mainTree.ParseFromIstream(&input)){
		std::cerr << "Could not read tree from input file." << std::endl;
		return -1;
	}

	std::cout << "Newick: " << mainTree.newick() << std::endl;
	std::cout << std::endl;
	std::cout << "nodes:\n";
	for(int i = 0; i < mainTree.nodes_size(); i++){
		std::cout << "node " << i << ":\n";
		std::cout << "nuc_mutation:\n";
		for(int j = 0; j < mainTree.nodes(i).nuc_mutation_size(); j++){
			std::cout << "(" << mainTree.nodes(i).nuc_mutation(j).position() << ", " << mainTree.nodes(i).nuc_mutation(j).gap_position() << ", "\
			<< mainTree.nodes(i).nuc_mutation(j).condensed() << ", " << mainTree.nodes(i).nuc_mutation(j).nucs() << ") ";
		}
		std::cout << "\n";
		for(int j = 0; j < mainTree.nodes(i).block_mutation().condensed_block_mut_size(); j++){
			std::cout << mainTree.nodes(i).block_mutation().condensed_block_mut(j) << " ";
		}
		std::cout << "\n";
	}

	std::cout << "\nGaps:\n";
	std::cout << "position:\n";
	for(int i = 0; i < mainTree.gaps().position_size(); i++){
		std::cout << mainTree.gaps().position(i) << " ";
	}
	std::cout << '\n';

	std::cout << "\nBlocks:\n";
	std::cout << "\nblockid:\n";
	for(int i = 0; i < mainTree.blocks().blockid_size(); i++){
		std::cout << mainTree.blocks().blockid(i) << " ";
	}
	std::cout << '\n';
	std::cout << "\nchromosome_name:\n";
	for(int i = 0; i < mainTree.blocks().chromosome_name_size(); i++){
		std::cout << mainTree.blocks().chromosome_name(i) << " ";
	}
	std::cout << '\n';

}