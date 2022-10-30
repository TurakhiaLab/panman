#ifndef PANGENOME_MAT_HPP
#define PANGENOME_MAT_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include "mutation_annotation.pb.h"


namespace PangenomeMAT {

    enum NucMutationType {
        NS = 0,
        ND = 1,
        NI = 2,
        NSNPS = 3,
        NSNPI = 4,
        NSNPD = 5,
    };

    enum BlockMutationType {
        BI = 0,
        BD = 1,
        NONE = 1000
    };

    void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);

    struct NucMut {
        NucMut(MAT::nuc_mut mutation){
            position = mutation.position();

            if(mutation.has_gap_position()){
                gapPosition = mutation.gap_position();
                // std::cout << '0';
            } else {
                gapPosition = -1;
                // std::cout << '1';
            }
            condensed = mutation.condensed();
            nucs = mutation.nucs();
        }

        uint32_t position;
        int32_t gapPosition;
        uint32_t condensed;
        uint64_t nucs;
    };

    struct BlockMut {

        void loadFromProtobuf( MAT::block_mut mutation ){
            for(int i = 0; i < mutation.condensed_block_mut_size(); i++){
                condensedBlockMut.push_back(mutation.condensed_block_mut(i));
            }
        }

        std::vector< uint32_t > condensedBlockMut;
    };

    struct Block {
        Block(MAT::block b);

        uint32_t blockId;
        std::vector< uint32_t > consensusSeq;
        std::string chromosomeName;
    };

    struct GapList {

        std::vector< uint32_t > position;
        // std::vector< uint32_t > condensed;
        std::vector< uint32_t > blockId;
        std::vector< uint32_t > gapLength;
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
        std::vector< std::string > annotations;
    };

    class Tree {
    private:
        size_t m_currInternalNode{ 0 };
        size_t m_numLeaves{ 0 };
        size_t m_maxDepth{ 0 };
        float m_meanDepth{ 0 };

        Node* createTreeFromNewickString(std::string newick);
        void assignMutationsToNodes(Node* root, size_t& currentIndex, std::vector< MAT::node >& nodes);
        int getTotalParsimony(NucMutationType nucMutType, BlockMutationType blockMutType = NONE);
        int getTotalParsimonyParallel(NucMutationType nucMutType, BlockMutationType blockMutType = NONE);

        std::unordered_map<std::string, Node*> allNodes;
        
        std::vector< Node* > allLeaves; // Probably temporary for testing

        std::string newInternalNodeId() {
            return "node_" + std::to_string(++m_currInternalNode);
        }

    public:
        Tree(std::ifstream& fin);
        void printSummary();
        void printFASTA(std::ofstream& fout); // Old. Delete when done with testing
        void printFASTA_updated(std::ofstream& fout);

        void writeToFile(std::ofstream& fout);
        void printBfs(); // Temporary function. To be removed later;
        void sampleWriteToFile(std::ofstream& fout); // Temporary function. To be removed later;

        Node* root;
        GapList gaps;
        std::vector< Block > blocks;
    };

};

#endif // PANGENOME_MAP_HPP