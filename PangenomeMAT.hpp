#include <vector>
#include <string>
#include <fstream>

namespace PangenomeMAT {

	struct NucMut {
		int32_t position;
		int32_t gapPosition;
		int32_t condensed; // The first two bits here should indicate the type of mutation and the other parameters would be read accordingly;
		int32_t nucs;
	};

	struct BlockMut {
		std::vector< int32_t > condensedBlockMuts;
	};

	struct BlockList {
		std::vector< int32_t > blockIdsAndChromosomeNames; // Since blockId takes 3 bytes and chromosomeName takes one byte, we can store it efficiently
	};

	struct GapList {
		std::vector< int32_t > position;
		std::vector< int32_t > condensed;
	};

	class Node {
	public:
		std::string id;
		Node* parent;
		std::vector< Node* > children;
		std::vector< NucMut > nucMutations;
		BlockMut blockMutations;
	};

	class Tree {
	private:
		Node* getTreeFromNewickString(std::string newick);
		std::unordered_map<std::string, Node*> nodes;
		
	public:
		Tree(std::ifstream& fin);

		Node* root;
		GapList gaps;
		BlockList blocks;
	};

};