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
		std::vector< int32_t > condensedBlockMut;
	};

	struct BlockList {
		std::vector< int32_t > blockId;
		std::vector< int32_t > chromosome_name;
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
		std::vector< NucMut > nucMutation;
		BlockMut blockMutation;
	};

	class Tree {
	private:
		Node* getTreeFromNewickString(std::string newick);
		std::unordered_map<std::string, Node*> node;
		
	public:
		Tree(std::ifstream& fin);

		Node* root;
		GapList gap;
		BlockList block;
	};

};