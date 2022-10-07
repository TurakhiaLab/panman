#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>

namespace PangenomeMAT {

	void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);

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
		Node(std::string id, float len);
		Node(std::string id, Node* par, float len);

		// Check this with professor
		size_t level;
    	float branchLength;

    	std::string identifier;
		Node* parent;
		std::vector< Node* > children;
		std::vector< NucMut > nucMutation;
		BlockMut blockMutation;
	};

	class Tree {
	private:
		Node* createTreeFromNewickString(std::string newick);
		std::unordered_map<std::string, Node*> node;
		size_t currInternalNode;

		std::string newInternalNodeId() {
        	return "node_" + std::to_string(++currInternalNode);
    	}

	public:
		Tree(std::ifstream& fin);

		Node* root;
		GapList gap;
		BlockList block;
	};

};