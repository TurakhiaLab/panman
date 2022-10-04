#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mutation_annotation.pb.h"

int main(){

	GOOGLE_PROTOBUF_VERIFY_VERSION;

	std::ofstream output("sample_mat.pb");

	MAT::tree mainTree;

	// First Node	
	MAT::nuc_mut firstNodeNucMutationOne;
	firstNodeNucMutationOne.set_position(2);
	firstNodeNucMutationOne.set_gap_position(18);
	firstNodeNucMutationOne.set_condensed(56);
	firstNodeNucMutationOne.set_nucs(7);

	MAT::nuc_mut firstNodeNucMutationTwo;
	firstNodeNucMutationTwo.set_position(7);
	firstNodeNucMutationTwo.set_gap_position(19);
	firstNodeNucMutationTwo.set_condensed(53);
	firstNodeNucMutationTwo.set_nucs(67);

	MAT::block_mut firstNodeBlockMutation;
	firstNodeBlockMutation.add_condensed_block_mut(34);
	firstNodeBlockMutation.add_condensed_block_mut(75);
	firstNodeBlockMutation.add_condensed_block_mut(33);

	MAT::node firstNode;
	firstNode.add_nuc_mutation();
	*firstNode.mutable_nuc_mutation(0) = firstNodeNucMutationOne;
	firstNode.add_nuc_mutation();
	*firstNode.mutable_nuc_mutation(1) = firstNodeNucMutationTwo;
	*firstNode.mutable_block_mutation() = firstNodeBlockMutation;

	// Second Node	
	MAT::nuc_mut secondNodeNucMutationOne;
	secondNodeNucMutationOne.set_position(34);
	secondNodeNucMutationOne.set_gap_position(43);
	secondNodeNucMutationOne.set_condensed(67);
	secondNodeNucMutationOne.set_nucs(84);

	MAT::nuc_mut secondNodeNucMutationTwo;
	secondNodeNucMutationTwo.set_position(32);
	secondNodeNucMutationTwo.set_gap_position(64);
	secondNodeNucMutationTwo.set_condensed(92);
	secondNodeNucMutationTwo.set_nucs(46);

	MAT::block_mut secondNodeBlockMutation;
	secondNodeBlockMutation.add_condensed_block_mut(56);
	secondNodeBlockMutation.add_condensed_block_mut(72);
	secondNodeBlockMutation.add_condensed_block_mut(58);

	MAT::node secondNode;
	secondNode.add_nuc_mutation();
	*secondNode.mutable_nuc_mutation(0) = secondNodeNucMutationOne;
	secondNode.add_nuc_mutation();
	*secondNode.mutable_nuc_mutation(1) = secondNodeNucMutationTwo;
	*secondNode.mutable_block_mutation() = secondNodeBlockMutation;

	MAT::gap_list treeGapList;
	treeGapList.add_position(23);
	treeGapList.add_position(19);
	treeGapList.add_position(1);	

	treeGapList.add_condensed(24);
	treeGapList.add_condensed(73);
	treeGapList.add_condensed(65);
	treeGapList.add_condensed(38);

	MAT::block_list treeBlockList;
	treeBlockList.add_blockid(1);
	treeBlockList.add_blockid(645);
	treeBlockList.add_blockid(6);

	treeBlockList.add_chromosome_name(2);
	treeBlockList.add_chromosome_name(0);
	treeBlockList.add_chromosome_name(1);

	mainTree.set_newick("abc");
	mainTree.add_nodes();
	*mainTree.mutable_nodes(0) = firstNode;
	mainTree.add_nodes();
	*mainTree.mutable_nodes(1) = secondNode;
	*mainTree.mutable_gaps() = treeGapList;
	*mainTree.mutable_blocks() = treeBlockList;

    if (!mainTree.SerializeToOstream(&output)) {
		std::cerr << "Failed to write output file." << std::endl;
		return -1;
    }

}