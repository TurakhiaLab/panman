#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <json/json.h>
#include "panman.capnp.h"
#include "common.hpp"

#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <kj/std/iostream.h>

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

    // Extract mutation from protobuf nucMut object
    NucMut(panman::NucMut::Reader mutation, int64_t blockId, bool blockGapExist) {
        nucPosition = mutation.getNucPosition();
        primaryBlockId = (blockId >> 32);
        mutInfo = (mutation.getMutInfo() & 0xFF);
        nucs = (mutation.getMutInfo() >> 8);
        nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

        if(blockGapExist) {
            secondaryBlockId = (blockId & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }

        if(mutation.getNucGapExist()) {
            nucGapPosition = mutation.getNucGapPosition();
        } else {
            nucGapPosition = -1;
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

    void loadFromProtobuf(panman::Mutation::Reader mutation) {
        primaryBlockId = (mutation.getBlockId() >> 32);
        if(mutation.getBlockGapExist()) {
            secondaryBlockId = (mutation.getBlockId() & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }
        blockMutInfo = mutation.getBlockMutInfo();
        // Whether the mutation is a block inversion or not. Inversion is marked by
        // `blockMutInfo = deletion` and `inversion = true`
        inversion = mutation.getBlockInversion();
    }

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
    void assignMutationsToNodes(Node* root, size_t& currentIndex,
                                std::vector<panman::Node::Reader>& storedNode);

    // Get the total number of mutations of given type
    int getTotalParsimonyParallel(NucMutationType nucMutType,
                                  BlockMutationType blockMutType = NONE);

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

    // Compress extracted subtree by combining parent and child nodes where parent has only
    // one child
    // void compressTreeParallel(Node* node, size_t level);
    void compressTreeParallel(Node* node, size_t level, const std::set< std::string >& nodeIdsToDefinitelyInclude);

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


    std::string newInternalNodeId() {
        return "node_" + std::to_string(++m_currInternalNode);
    }

    size_t m_currInternalNode{ 0 };
    size_t m_numLeaves{ 0 };
    size_t m_maxDepth{ 0 };
    float m_meanDepth{ 0 };

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

    Tree(const panman::Tree::Reader& mainTree);
    Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
    Tree(std::ifstream& fin, std::ifstream& secondFin,
         FILE_TYPE ftype = FILE_TYPE::GFA, std::string reference = "");

    // Copy blocks from current tree into new tree which is rooted at one of the internal
    // nodes of the current tree. Used in split for PanMAN
    Tree(Node* newRoot, const std::vector< Block >& b, const std::vector< GapList >& g,
         std::unordered_map< std::string, int >& cs,
         std::unordered_map< std::string, int >& ri,
         std::unordered_map< std::string, bool >& si,
         const BlockGapList& bgl);

    void protoMATToTree(const panman::Tree::Reader& mainTree);

    // Fitch Algorithm on Nucleotide mutations
    int nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states, int refState=-1);
    int nucFitchForwardPassOpt(Node* node, std::unordered_map< std::string, int >& states);
    // Default state is used in rerooting to a tip sequence. It is used to fix the state at
    // the root
    void nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states,
                              int parentState, int defaultState = (1<<28));
    void nucFitchBackwardPassOpt(Node* node, std::unordered_map< std::string, int >& states,
                                 int parentState, int defaultState = (1<<28));
    void nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states,
                                 std::unordered_map< std::string,
                                 std::pair< panmanUtils::NucMutationType, char > >& mutations,
                                 int parentState);
    void nucFitchAssignMutationsOpt(Node* node, std::unordered_map< std::string, int >& states,
                                    std::unordered_map< std::string,
                                    std::pair< panmanUtils::NucMutationType, char > >& mutations,
                                    int parentState);

    // Sankoff algorithm on Nucleotide Mutations
    std::vector< int > nucSankoffForwardPass(Node* node, std::unordered_map< std::string, std::vector< int > >& stateSets);
    std::vector< int > nucSankoffForwardPassOpt(Node* node, std::unordered_map< std::string, std::vector< int > >& stateSets);
    void nucSankoffBackwardPass(Node* node,
                                std::unordered_map< std::string, std::vector< int > >& stateSets,
                                std::unordered_map< std::string, int >& states, int parentPtr,
                                int defaultValue = (1<<28));
    void nucSankoffBackwardPassOpt(Node* node,
                                   std::unordered_map< std::string, std::vector< int > >& stateSets,
                                   std::unordered_map< std::string, int >& states, int parentPtr,
                                   int defaultValue = (1<<28));
    void nucSankoffAssignMutations(Node* node,
                                   std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                   std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState);
    void nucSankoffAssignMutationsOpt(Node* node,
                                      std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                      std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState);

    // Fitch algorithm on Block Mutations
    int blockFitchForwardPassNew(Node* node,
                                 std::unordered_map< std::string, int >& states);
    void blockFitchBackwardPassNew(Node* node,
                                   std::unordered_map< std::string, int >& states, int parentState,
                                   int defaultValue = (1<<28));
    void blockFitchAssignMutationsNew(Node* node,
                                      std::unordered_map< std::string, int >& states,
                                      std::unordered_map< std::string,
                                      std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState);

    // Sankoff algorithm on Block Mutations
    std::vector< int > blockSankoffForwardPass(Node* node, std::unordered_map< std::string,
            std::vector< int > >& stateSets);
    void blockSankoffBackwardPass(Node* node,
                                  std::unordered_map< std::string, std::vector< int > >& stateSets,
                                  std::unordered_map< std::string, int >& states, int parentPtr,
                                  int defaultValue = (1<<28));
    void blockSankoffAssignMutations(Node* node,
                                     std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                     std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState);

    // void printSummary();
    void printSummary(std::ostream &out);
    void printBfs(Node* node = nullptr);
    void printFASTA(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);
    void printFASTANew(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);
    void printSingleNode(std::ostream& fout, const sequence_t& sequence,
                                         const blockExists_t& blockExists, const blockStrand_t& blockStrand,
                                         std::string nodeIdentifier, std::tuple< int, int, int, int > &panMATStart, std::tuple< int, int, int, int > &panMATEnd);
    void printFASTAParallel(std::ostream& fout, bool aligned = false);
    void printMAF(std::ostream& fout);

    void printMAFNew(std::ostream& fout);
    void generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout);
    void printVCFParallel(std::string reference, std::ostream& fout);
    void extractAminoAcidTranslations(std::ostream& fout, int64_t start, int64_t end);

    // Extract PanMAT representing a segment of the genome. The start and end coordinates
    // are with respect to the root sequence. The strands of the terminal blocks in all
    // sequences are assumed to be the same as their strands in the root sequence for the
    // purpose of splitting the terminal blocks during extraction
    void extractPanMATSegment(kj::std::StdOutputStream& fout, int64_t start, int64_t end);
    void extractPanMATIndex(std::ostream& fout, int64_t start, int64_t end, std::string nodeIdentifier, bool single=true);

    Node* subtreeExtractParallel(std::vector< std::string > nodeIds, const std::set< std::string >& nodeIdsToDefinitelyInclude = {});
    // Node* subtreeExtractParallel(std::vector< std::string > nodeIds);
    void writeToFile(kj::std::StdOutputStream& fout, Node* node = nullptr);
    std::string getNewickString(Node* node);
    std::string getStringFromReference(std::string reference, bool aligned = true,
                                       bool incorporateInversions=true);
    void getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists,
                                  blockStrand_t& blockStrand, std::string reference, bool rotateSequence = false,
                                  int* rotIndex = nullptr);

    // For each node in the tree, print mutations with respect to the root node to the
    // output file
    void printMutations(std::ostream& fout);
    void printMutationsNew(std::ostream& fout);
    void printNodePaths(std::ostream& fout);

    void getBlockSequenceFromReference(block_t& sequence, bool& blockExists,
                                       bool& blockStrand, std::string reference, int64_t primaryBlockId,
                                       int64_t secondaryBlockId);

    // Split file provided as input.
    std::pair< Tree, Tree > splitByComplexMutations(const std::string& nodeId3);

    // get unaligned global coordinate
    int32_t getUnalignedGlobalCoordinate(int32_t primaryBlockId, int32_t secondaryBlockId,
                                         int32_t pos, int32_t gapPos, const sequence_t& sequence,
                                         const blockExists_t& blockExists, const blockStrand_t& blockStrand,
                                         int circularOffset = 0, bool * check = nullptr);
    std::tuple< int, int, int, int > globalCoordinateToBlockCoordinate(
        int64_t globalCoordinate,
        const sequence_t& sequence,
        const blockExists_t& blockExists,
        const blockStrand_t& blockStrand, int64_t circularOffset = 0);

    std::string getSequenceFromVCF(std::string sequenceId, std::ifstream& fin);
    bool verifyVCFFile(std::ifstream& fin);
    void vcfToFASTA(std::ifstream& fin, std::ofstream& fout);
    void annotate(std::ifstream& fin);
    std::vector< std::string > searchByAnnotation(std::string annotation);
    void convertToGFA(std::ostream& fout);
    void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);
    void getNodesPreorder(panmanUtils::Node* root, capnp::List<panman::Node>::Builder& nodesBuilder, size_t& nodeIndex);
    size_t getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition,
                               int nucGapPosition);

    // Transforms tree such that given node becomes child of new root
    void transform(Node* node);
    void reroot(std::string sequenceName);

    

};

// Represents complex mutations like Horizontal Gene Transfer or Recombinations
struct ComplexMutation {
    char mutationType;
    size_t treeIndex1, treeIndex2, treeIndex3;
    std::string sequenceId1, sequenceId2, sequenceId3;

    // coordinates of start in parent 1
    int32_t primaryBlockIdStart1;
    int32_t secondaryBlockIdStart1;
    int32_t nucPositionStart1;
    int32_t nucGapPositionStart1;

    // coordinates of end in parent 1
    int32_t primaryBlockIdEnd1;
    int32_t secondaryBlockIdEnd1;
    int32_t nucPositionEnd1;
    int32_t nucGapPositionEnd1;

    // coordinates of start in parent 2
    int32_t primaryBlockIdStart2;
    int32_t secondaryBlockIdStart2;
    int32_t nucPositionStart2;
    int32_t nucGapPositionStart2;

    // coordinates of end in parent 2
    int32_t primaryBlockIdEnd2;
    int32_t secondaryBlockIdEnd2;
    int32_t nucPositionEnd2;
    int32_t nucGapPositionEnd2;

    ComplexMutation(char mutType, int tIndex1, int tIndex2, int tIndex3, std::string sId1,
                    std::string sId2, std::string sId3, std::tuple< int,int,int,int > t1,
                    std::tuple< int,int,int,int > t2, std::tuple< int,int,int,int > t3,
                    std::tuple< int,int,int,int > t4) {
        mutationType = mutType;
        treeIndex1 = tIndex1;
        treeIndex2 = tIndex2;
        treeIndex3 = tIndex3;

        sequenceId1 = sId1;
        sequenceId2 = sId2;
        sequenceId3 = sId3;

        primaryBlockIdStart1 = std::get<0>(t1);
        secondaryBlockIdStart1 = std::get<1>(t1);
        nucPositionStart1 = std::get<2>(t1);
        nucGapPositionStart1 = std::get<3>(t1);

        primaryBlockIdEnd1 = std::get<0>(t2);
        secondaryBlockIdEnd1 = std::get<1>(t2);
        nucPositionEnd1 = std::get<2>(t2);
        nucGapPositionEnd1 = std::get<3>(t2);

        primaryBlockIdStart2 = std::get<0>(t3);
        secondaryBlockIdStart2 = std::get<1>(t3);
        nucPositionStart2 = std::get<2>(t3);
        nucGapPositionStart2 = std::get<3>(t3);

        primaryBlockIdEnd2 = std::get<0>(t4);
        secondaryBlockIdEnd2 = std::get<1>(t4);
        nucPositionEnd2 = std::get<2>(t4);
        nucGapPositionEnd2 = std::get<3>(t4);
    }

    ComplexMutation(panman::ComplexMutation::Reader cm) {
        mutationType = (cm.getMutationType()? 'H': 'R');
        treeIndex1 = cm.getTreeIndex1();
        treeIndex2 = cm.getTreeIndex2();
        treeIndex3 = cm.getTreeIndex3();
        sequenceId1 = cm.getSequenceId1();
        sequenceId2 = cm.getSequenceId2();
        sequenceId3 = cm.getSequenceId3();

        primaryBlockIdStart1 = (cm.getBlockIdStart1() >> 32);
        secondaryBlockIdStart1 = (cm.getBlockGapExistEnd1()?
                                  (cm.getBlockIdStart1()&(0xFFFFFFFF)): -1);
        nucPositionStart1 = cm.getNucPositionStart1();
        nucGapPositionStart1 = (cm.getNucGapExistStart1()? (cm.getNucGapPositionStart1()) : -1);

        primaryBlockIdStart2 = (cm.getBlockIdStart2() >> 32);
        secondaryBlockIdStart2 = (cm.getNucGapExistStart2()?
                                  (cm.getBlockIdStart2()&(0xFFFFFFFF)): -1);
        nucPositionStart2 = cm.getNucPositionStart2();
        nucGapPositionStart2 = (cm.getNucGapExistStart2()? (cm.getNucGapPositionStart2()) : -1);

        primaryBlockIdEnd1 = (cm.getBlockIdEnd1() >> 32);
        secondaryBlockIdEnd1 = (cm.getBlockGapExistEnd1()? (cm.getBlockIdEnd1()&(0xFFFFFFFF)): -1);
        nucPositionEnd1 = cm.getNucPositionEnd1();
        nucGapPositionEnd1 = (cm.getNucGapExistEnd1()? (cm.getNucGapPositionEnd1()) : -1);

        primaryBlockIdEnd2 = (cm.getBlockIdEnd2() >> 32);
        secondaryBlockIdEnd2 = (cm.getBlockGapExistEnd2()? (cm.getBlockIdEnd2()&(0xFFFFFFFF)): -1);
        nucPositionEnd2 = cm.getNucPositionEnd2();
        nucGapPositionEnd2 = (cm.getNucGapExistEnd2()? (cm.getNucGapPositionEnd2()) : -1);
    }

    panman::ComplexMutation::Builder toCapnProto() {
        panman::ComplexMutation::Builder cm(nullptr);
        cm.setMutationType(mutationType == 'H');
        cm.setTreeIndex1(treeIndex1);
        cm.setTreeIndex2(treeIndex2);
        cm.setTreeIndex3(treeIndex3);
        cm.setSequenceId1(sequenceId1);
        cm.setSequenceId2(sequenceId2);
        cm.setSequenceId3(sequenceId3);

        if(secondaryBlockIdStart1 != -1) {
            cm.setBlockGapExistStart1(true);
            cm.setBlockIdStart1(((int64_t)primaryBlockIdStart1 << 32)+secondaryBlockIdStart1);
        } else {
            cm.setBlockGapExistStart1(false);
            cm.setBlockIdStart1(((int64_t)primaryBlockIdStart1 << 32));
        }
        cm.setNucPositionStart1(nucPositionStart1);

        if(nucGapPositionStart1 != -1) {
            cm.setNucGapExistStart1(true);
            cm.setNucGapPositionStart1(nucGapPositionStart1);
        }

        if(secondaryBlockIdStart2 != -1) {
            cm.setBlockGapExistStart2(true);
            cm.setBlockIdStart2(((int64_t)primaryBlockIdStart2 << 32)+secondaryBlockIdStart2);
        } else {
            cm.setBlockGapExistStart2(false);
            cm.setBlockIdStart2(((int64_t)primaryBlockIdStart2 << 32));
        }
        cm.setNucPositionStart2(nucPositionStart2);

        if(nucGapPositionStart2 != -1) {
            cm.setNucGapExistStart2(true);
            cm.setNucGapPositionStart2(nucGapPositionStart2);
        }

        if(secondaryBlockIdEnd1 != -1) {
            cm.setBlockGapExistEnd1(true);
            cm.setBlockIdEnd1(((int64_t)primaryBlockIdEnd1 << 32)+secondaryBlockIdEnd1);
        } else {
            cm.setBlockGapExistEnd1(false);
            cm.setBlockIdEnd1(((int64_t)primaryBlockIdEnd1 << 32));
        }
        cm.setNucPositionEnd1(nucPositionEnd1);

        if(nucGapPositionEnd1 != -1) {
            cm.setNucGapExistEnd1(true);
            cm.setNucGapPositionEnd1(nucGapPositionEnd1);
        }

        if(secondaryBlockIdEnd2 != -1) {
            cm.setBlockGapExistEnd2(true);
            cm.setBlockIdEnd2(((int64_t)primaryBlockIdEnd2 << 32)+secondaryBlockIdEnd2);
        } else {
            cm.setBlockGapExistEnd2(false);
            cm.setBlockIdEnd2(((int64_t)primaryBlockIdEnd2 << 32));
        }
        cm.setNucPositionEnd2(nucPositionEnd2);

        if(nucGapPositionEnd2 != -1) {
            cm.setNucGapExistEnd2(true);
            cm.setNucGapPositionEnd2(nucGapPositionEnd2);
        }

        return cm;
    }

};

// Data structure to represent PanMAN
class TreeGroup {
  public:
    // List of PanMATs in PanMAN
    std::vector< Tree > trees;
    // List of complex mutations linking PanMATs
    std::vector< ComplexMutation > complexMutations;

    TreeGroup(std::istream& fin);
    // List of PanMAT files and a file with all the complex mutations relating these files
    TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile);
    TreeGroup(std::vector< Tree* >& t);
    TreeGroup(std::vector< Tree* >& tg, std::ifstream& mutationFile);

    TreeGroup* subnetworkExtract(std::unordered_map< int, std::vector< std::string > >& nodeIds);

    void printFASTA(std::ofstream& fout, bool rootSeq = false);
    void writeToFile(kj::std::StdOutputStream& fout);
    void printComplexMutations(std::ostream& fout);
};

};