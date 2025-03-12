#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <cstring>


#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"


static const int SANKOFF_INF = 100000001;

typedef std::vector< // block id
            std::pair< 
                std::vector< std::pair< char, std::vector< char > > >, // vector - nuc id, char - char at nuc id, vector <char> - nuc at gap id
                std::vector< 
                    std::vector< std::pair< char, std::vector< char > > > 
                > 
            > 
        > sequence_t;
// Individual block
typedef std::vector< std::pair< char, std::vector< char > > > block_t;

typedef  std::vector< std::pair< bool, std::vector< bool > > > blockExists_t;
// Forward or reverse strand
typedef  std::vector< std::pair< bool, std::vector< bool > > > blockStrand_t;

namespace panmanUtils {

enum FILE_TYPE {
    PANMAT = 0,
    GFA = 1,
    PANGRAPH=2,
    MSA = 3,
    MSA_OPTIMIZE = 4,
    // FASTA = 5
    };
};

namespace panmanUtils {

    enum NucMutationType {
        // Nucleotide Substutution
        NS = 0,
        // Nucleotide Deletion
        ND = 1,
        // Nucleotide Insertion
        NI = 2,
        // Single Nucleotide Substitution
        NSNPS = 3,
        // Single Nucleotide Insertion
        NSNPI = 4,
        // Single Nucleotide Deletion
        NSNPD = 5,
        // None
        NNONE = 2000
    };
    
    enum BlockMutationType {
        // Block Insertion
        BI = 1,
        // Block Deletion
        BD = 0,
        // Block Inversion
        BIn = 2,
        // None
        NONE = 1000
    };
    
    // Struct for representing Nucleotide Mutation
    struct NucMut {
        int32_t nucPosition;
        int32_t nucGapPosition;
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        uint8_t mutInfo;
        uint32_t nucs;
    
        // Create SNP mutation for MSA (optimized for memory)
        NucMut( const std::tuple< int, int8_t, int8_t>& mutationInfo ) {
            // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
            primaryBlockId = 0;
            secondaryBlockId = -1;
            nucPosition = std::get<0>(mutationInfo);
            nucGapPosition = -1;
            mutInfo = (int)std::get<1>(mutationInfo) + (1 << 4);
            nucs = ((int)std::get<2>(mutationInfo) << 20);
        }
        
        // Create SNP mutation
        NucMut( const std::tuple< int, int, int, int, int, int >& mutationInfo ) {
            // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
            primaryBlockId = std::get<0>(mutationInfo);
            secondaryBlockId = std::get<1>(mutationInfo);
            nucPosition = std::get<2>(mutationInfo);
            nucGapPosition = std::get<3>(mutationInfo);
            mutInfo = std::get<4>(mutationInfo) + (1 << 4);
            nucs = (std::get<5>(mutationInfo) << 20);
        }
    
        // Create non-SNP mutations from SNP mutations at consecutive positions for MSA
        NucMut(const std::vector< std::tuple< int, int8_t, int8_t > >& mutationArray,
               int start, int end) {
            primaryBlockId = 0;
            secondaryBlockId = -1;
    
            mutInfo = ((end - start) << 4);
            // type
            switch(std::get<1>(mutationArray[start])) {
            case panmanUtils::NucMutationType::NSNPS:
                mutInfo += panmanUtils::NucMutationType::NS;
                break;
            case panmanUtils::NucMutationType::NSNPI:
                mutInfo += panmanUtils::NucMutationType::NI;
                break;
            case panmanUtils::NucMutationType::NSNPD:
                mutInfo += panmanUtils::NucMutationType::ND;
                break;
            case panmanUtils::NucMutationType::NS:
                mutInfo += panmanUtils::NucMutationType::NS;
                break;
            case panmanUtils::NucMutationType::NI:
                mutInfo += panmanUtils::NucMutationType::NI;
                break;
            case panmanUtils::NucMutationType::ND:
                mutInfo += panmanUtils::NucMutationType::ND;
                break;
            }
    
            nucPosition = (int)std::get<0>(mutationArray[start]);
            nucGapPosition = -1;
    
            nucs = 0;
            for(int i = start; i < end; i++) {
                nucs += (std::get<2>(mutationArray[i]) << (4*(5-(i - start))));
            }
    
            // if (nucPosition == 0){
            //     std::cout << "\t Writing " << nucPosition << " " << 
            //                 (int)mutInfo << " " << 
            //                 nucs << " " <<
            //                 std::endl;
            // }
        }
    
        // Create non-SNP mutations from SNP mutations at consecutive positions
        NucMut(const std::vector< std::tuple< int, int, int, int, int, int > >& mutationArray,
               int start, int end) {
            primaryBlockId = std::get<0>(mutationArray[start]);
            secondaryBlockId = std::get<1>(mutationArray[start]);
    
            mutInfo = ((end - start) << 4);
            // type
            switch(std::get<4>(mutationArray[start])) {
            case panmanUtils::NucMutationType::NSNPS:
                mutInfo += panmanUtils::NucMutationType::NS;
                break;
            case panmanUtils::NucMutationType::NSNPI:
                mutInfo += panmanUtils::NucMutationType::NI;
                break;
            case panmanUtils::NucMutationType::NSNPD:
                mutInfo += panmanUtils::NucMutationType::ND;
                break;
            case panmanUtils::NucMutationType::NS:
                mutInfo += panmanUtils::NucMutationType::NS;
                break;
            case panmanUtils::NucMutationType::NI:
                mutInfo += panmanUtils::NucMutationType::NI;
                break;
            case panmanUtils::NucMutationType::ND:
                mutInfo += panmanUtils::NucMutationType::ND;
                break;
            }
    
            nucPosition = std::get<2>(mutationArray[start]);
            nucGapPosition = std::get<3>(mutationArray[start]);
    
            nucs = 0;
            for(int i = start; i < end; i++) {
                nucs += (std::get<5>(mutationArray[i]) << (4*(5-(i - start))));
            }
        }
    
        
    };
    
    // Struct for representing Block Mutations
    struct BlockMut {
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
    
        // Whether mutation is an insertion or deletion - Strand inversions are marked by
        // `blockMutInfo=false`, but they are not deletions
        bool  blockMutInfo;
    
        // Whether the block is being inverted or not. In case of insertion, whether the inserted
        // block is inverted or not
        bool inversion;
    
        BlockMut(size_t blockId, std::pair< BlockMutationType, bool > type, int secondaryBId = -1) {
            primaryBlockId = blockId;
            secondaryBlockId = secondaryBId;
            if(type.first == BlockMutationType::BI) {
                blockMutInfo = true;
            } else {
                // blockMutInfo is also set to false in the case of inversions. If the mutation
                // isn't an inversion, `blockMutInfo = false` indicates deletion
                blockMutInfo = false;
            }
    
            if(type.second == 1) {
                // If type.second == 1 (inversion)
                inversion = true;
            } else {
                inversion = false;
            }
        }
    
        BlockMut() {}
    
        
    };
    
    // List of default blocks in the global coordinate system of the PanMAT
    struct Block {
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
    
        std::vector< uint32_t > consensusSeq;
        std::string chromosomeName;
    
        Block(size_t primaryBlockId, std::string seq);
        // seq is a compressed form of the sequence where each nucleotide is stored in 4 bytes
        Block(int32_t primaryBlockId, int32_t secondaryBlockId, const std::vector< uint32_t >& seq);  
    };
    
    // List of gaps in the global coordinate system of the PanMAT
    struct GapList {
        std::vector< uint32_t > nucPosition;
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        std::vector< uint32_t > nucGapLength;
    
    };
    
    
    
    // @DEPRECATED. To be removed when secondary block ID is removed
    struct BlockGapList {
        std::vector< uint32_t > blockPosition;
        std::vector< uint32_t > blockGapLength;
    };
    
    
    // PanMAT tree node
    class Node {
      public:
        float branchLength = 0.0;
        size_t level;
        std::string identifier;
        Node* parent;
        std::vector< Node* > children;
        std::vector< NucMut > nucMutation;
        std::vector< BlockMut > blockMutation;
        std::vector< std::string > annotations;
        bool isComMutHead = false;
        int treeIndex = -1;
    
        Node(std::string id, float len);
        Node(std::string id, Node* par, float len);
    };
    
    
    // Data structure to represent a PangenomeMAT
    class Tree {
      private:
        Node* createTreeFromNewickString(std::string newick);
    
        // In the proto file, nodes are stored in preorder. Once the tree has been generated in
        // memory, assign mutations from the proto file to the tree nodes using preorder
        // traversal
        // Get the total number of mutations of given type
        int getTotalParsimonyParallel(NucMutationType nucMutType,
                                      BlockMutationType blockMutType = NONE);
        
        void getBlockMutationsParallel();
    
        // Run tree traversal to extract mutations in range
        panmanUtils::Node* extractPanMATSegmentHelper(panmanUtils::Node* node,
                const std::tuple< int, int, int, int >& start,
                const std::tuple< int, int, int, int >& end, const blockStrand_t& rootBlockStrand);
    
        // Tree traversal for FASTA writer
        void printFASTAHelper(panmanUtils::Node* root, sequence_t& sequence,
                              blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout,
                              bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
        void printFASTAHelperNew(panmanUtils::Node* root, 
                              std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                              std::vector<bool>& blockExists, 
                              std::vector<bool>& blockStrand, 
                              std::ostream& fout,
                              bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
        
        std::string printFASTAUltraFastHelper(
                              const std::vector<bool>& blockSequence,
                              std::unordered_map<int, int>& blockLengths,
                              const std::vector<panmanUtils::Node*>& nodesFromTipToRoot,  
                              std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                              std::vector<bool>& blockExists, 
                              std::vector<bool>& blockStrand, 
                              bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
        
        std::pair<std::vector<std::string>, std::vector<int>> extractSequenceHelper(
                              const std::vector<bool>& blockSequence,
                              std::unordered_map<int, int>& blockLengths,
                              const std::vector<panmanUtils::Node*>& nodesFromTipToRoot,  
                              std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                              std::vector<bool>& blockExists, 
                              std::vector<bool>& blockStrand, 
                              bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
        
        std::pair<std::vector<std::string>, std::vector<int>> extractSingleSequence(panmanUtils::Node* node, bool aligned=false, bool rootSeq=false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
        
        void printSingleNodeHelper(std::vector<panmanUtils::Node*> &nodeList, int nodeListIndex, sequence_t& sequence,
            blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart={-1,-1,-1,-1}, const std::tuple< int, int, int, int >& panMATEnd={-1,-1,-1,-1});
    
        // Merge parent and child nodes when compressing subtree
        void mergeNodes(Node* par, Node* chi);
    
        // Used to combine their mutations at corresponding positions when parent and child
        // nodes are combined
        std::pair< int, int > replaceMutation(std::pair<int,int> oldMutation,
                                              std::pair<int, int> newMutation);
    
        // Iterate through mutations and combine mutations at the same position
        std::vector< NucMut > consolidateNucMutations(const std::vector< NucMut >& nucMutation);
    
        // Used to confirm that consolidateNucMutations worked correctly. Can be removed in
        // production
        bool debugSimilarity(const std::vector< NucMut > array1,
                             const std::vector< NucMut > array2);
    
    
        // Used in rerooting
        void dfsExpansion(Node* node, std::vector< Node* >& vec);
        Node* transformHelper(Node* node);
        void adjustLevels(Node* node);
    
        // Check if tree is a polytomy
        bool hasPolytomy(Node* node);
    
        // Check if one PanMAT coordinate is greater than or equal to the other. Only the strand
        // information of the first block needs to be provided because if the block IDs are
        // different, the strand information does not change the result
        bool panMATCoordinateGeq(const std::tuple< int, int, int, int >& coor1,
                                 const std::tuple< int, int, int, int >& coor2, bool strand);
    
        // Check if one PanMAT coordinate is less than or equal to the other. Only the strand
        // information of the first block needs to be provided because if the block IDs are
        // different, the strand information does not change the result
        bool panMATCoordinateLeq(const std::tuple< int, int, int, int >& coor1,
                                 const std::tuple< int, int, int, int >& coor2, bool strand);
    
    
        
    
        std::unordered_map<std::string, std::vector< std::string > > annotationsToNodes;
      public:
        Node *root;
        std::vector< Block > blocks;
        std::vector< GapList > gaps;
    
        // @DEPRECATED: To be removed with secondary block ID
        BlockGapList blockGaps;
    
        // Specifies the circular offset required to print the original sequence
        std::unordered_map< std::string, int > circularSequences;
    
        // Specifies the block by which the rotation algorithm rotated the sequence
        std::unordered_map< std::string, int > rotationIndexes;
    
        // Specifies whether sequence is inverted or not by the rotation algorithm
        std::unordered_map< std::string, bool > sequenceInverted;
    
        std::unordered_map< std::string, Node* > allNodes;
    
        Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
        Tree(std::ifstream& fin, std::ifstream& secondFin,
             FILE_TYPE ftype = FILE_TYPE::GFA, std::string reference = "");    
        
        Tree(std::ifstream& newick_ifstream);

        
    
        size_t m_currInternalNode{ 0 };
        size_t m_numLeaves{ 0 };
        size_t m_maxDepth{ 0 };
        float m_meanDepth{ 0 };

        std::string newInternalNodeId() {
            return "node_" + std::to_string(++m_currInternalNode);
        }
    };

    
    
    
};