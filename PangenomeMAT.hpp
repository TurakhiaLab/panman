#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include "mutation_annotation.pb.h"


namespace PangenomeMAT {

	enum MutationType {
		S = 0,
		D = 1,
		I = 2,
		SNP = 3
	};

	void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);

	struct NucMut {
		NucMut(MAT::nuc_mut mutation){
			position = mutation.position();
			gapPosition = mutation.gap_position();
			condensed = mutation.condensed();
			nucs = mutation.nucs();
		}

		int32_t position;
		int32_t gapPosition;
		int32_t condensed; // The first two bits here should indicate the type of mutation and the other parameters would be read accordingly;
		int32_t nucs;
	};

	struct BlockMut {

		void loadFromProtobuf( MAT::block_mut mutation ){
			for(int i = 0; i < mutation.condensed_block_mut_size(); i++){
				condensedBlockMut.push_back(mutation.condensed_block_mut(i));
			}
		}

		std::vector< int32_t > condensedBlockMut;
	};

	struct Block {
		uint32_t blockId;
		std::vector< uint32_t > consensusSeq;
		uint32_t chromosomeName;
	};

	struct GapList {
		std::vector< int32_t > position;
		std::vector< int32_t > condensed;
	};

	class Node {
	public:
		Node(std::string id, float len);
		Node(std::string id, Node* par, float len);

		// Check this with professor
    	float branchLength;
    	size_t level;

    	std::string identifier;
		Node* parent;
		std::vector< Node* > children;
		std::vector< NucMut > nucMutation;
		BlockMut blockMutation;
	};

	class Tree {
	private:
		size_t currInternalNode;
		size_t numLeaves;
		size_t maxDepth;
		float meanDepth;

		Node* createTreeFromNewickString(std::string newick);
		void assignMutationsToNodes(Node* root, size_t currentIndex, std::vector< MAT::node >& nodes);
		int getTotalParsimony(MutationType type);
		int getTotalParsimonyParallel(MutationType type);

		std::unordered_map<std::string, Node*> allNodes;

		std::string newInternalNodeId() {
        	return "node_" + std::to_string(++currInternalNode);
    	}

	public:
		Tree(std::ifstream& fin);
		void printSummary();
		void printBfs(); // Temporary function. To be removed later;

		Node* root;
		GapList gap;
		std::vector< Block > block;
	};

};