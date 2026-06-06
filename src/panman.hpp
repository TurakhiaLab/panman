#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <json/json.h>
#include "panman.capnp.h"
#include "panman.pb.h"
#include "common.hpp"

#include <mutex>

#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <kj/std/iostream.h>

std::string reconstructNewick(const panman::Tree::Reader& tree);

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
        const auto& codec = getCodec(getActiveAlphabet());
        const int payloadShift = codec.bitsPerSymbol * ((24 / codec.bitsPerSymbol) - 1);
        // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
        primaryBlockId = 0;
        secondaryBlockId = -1;
        nucPosition = std::get<0>(mutationInfo);
        nucGapPosition = -1;
        mutInfo = (int)std::get<1>(mutationInfo) + (1 << 4);
        nucs = ((int)std::get<2>(mutationInfo) << payloadShift);
    }
    
    // Create SNP mutation
    NucMut( const std::tuple< int, int, int, int, int, int >& mutationInfo ) {
        const auto& codec = getCodec(getActiveAlphabet());
        const int payloadShift = codec.bitsPerSymbol * ((24 / codec.bitsPerSymbol) - 1);
        // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
        primaryBlockId = std::get<0>(mutationInfo);
        secondaryBlockId = std::get<1>(mutationInfo);
        nucPosition = std::get<2>(mutationInfo);
        nucGapPosition = std::get<3>(mutationInfo);
        mutInfo = std::get<4>(mutationInfo) + (1 << 4);
        nucs = (std::get<5>(mutationInfo) << payloadShift);
    }

    // Create non-SNP mutations from SNP mutations at consecutive positions for MSA
    NucMut(const std::vector< std::tuple< int, int8_t, int8_t > >& mutationArray,
           int start, int end) {
        const auto& codec = getCodec(getActiveAlphabet());
        const int payloadCapacity = 24 / codec.bitsPerSymbol;
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
            nucs += (std::get<2>(mutationArray[i])
                    << (codec.bitsPerSymbol * (payloadCapacity - 1 - (i - start))));
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
        const auto& codec = getCodec(getActiveAlphabet());
        const int payloadCapacity = 24 / codec.bitsPerSymbol;
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
            nucs += (std::get<5>(mutationArray[i])
                    << (codec.bitsPerSymbol * (payloadCapacity - 1 - (i - start))));
        }
    }

    // Extract mutation from protobuf nucMut object
    NucMut(panman::NucMut::Reader mutation, int64_t blockId, bool blockGapExist) {
        const uint32_t rawMutInfo = mutation.getMutInfo();
        nucPosition = mutation.getNucPosition();
        primaryBlockId = (blockId >> 32);
        mutInfo = (rawMutInfo & 0xFF);
        nucs = (rawMutInfo >> 8);
        nucs = expandMutationPayloadFromWire(nucs, static_cast<uint8_t>(mutInfo >> 4),
                                            getActiveAlphabet());

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

    NucMut(panmanOld::nucMut mutation, int64_t blockId, bool blockGapExist) {
        nucPosition = mutation.nucposition();
        primaryBlockId = (blockId >> 32);
        mutInfo = (mutation.mutinfo() & 0xFF);
        nucs = (mutation.mutinfo() >> 8);
        nucs = expandMutationPayloadFromWire(nucs, static_cast<uint8_t>(mutInfo >> 4),
                                            getActiveAlphabet());

        if(blockGapExist) {
            secondaryBlockId = (blockId & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }

        if(mutation.nucgapexist()) {
            nucGapPosition = mutation.nucgapposition();
        } else {
            nucGapPosition = -1;
        }
    }

    
};

// Struct for representing Block Mutations
struct BlockMut {
    int32_t primaryBlockId;
    int32_t secondaryBlockId;
    int chrIdx=0;

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
        /* if chrIdx exists, set it */
        chrIdx = mutation.getChrIdx();

    }

    void loadFromProtobuf(panmanOld::mutation mutation) {
        primaryBlockId = (mutation.blockid() >> 32);
        if(mutation.blockgapexist()) {
            secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }
        blockMutInfo = mutation.blockmutinfo();
        // Whether the mutation is a block inversion or not. Inversion is marked by
        // `blockMutInfo = deletion` and `inversion = true`
        inversion = mutation.blockinversion();
    }

    BlockMut(size_t blockId, std::pair< BlockMutationType, bool > type, int secondaryBId = -1, int chrIdxInput = -1) {
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

        if (chrIdxInput!=-1) chrIdx=chrIdxInput;
    }

    BlockMut() {}

    
};

// List of default blocks in the global coordinate system of the PanMAT
struct Block {
    int32_t primaryBlockId;
    int32_t secondaryBlockId;
    
    int64_t blockLength;
    
    std::vector< uint32_t > consensusSeq;
    std::string chromosomeName;

    Block(size_t primaryBlockId, std::string seq, Alphabet alphabet = getActiveAlphabet(), int64_t blockLength_ = 0);
    // seq is a compressed form of the sequence where each nucleotide is stored in 4 bytes
    Block(int32_t primaryBlockId, int32_t secondaryBlockId, const std::vector< uint32_t >& seq, int64_t blockLength_ = 0);  
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

// Block graph derived from a MAF. Block homology comes directly from the MAF:
// each alignment record is grouped (by identical gapped consensus) into a block
// "group". For every sequence we keep its block occurrences in MAF-coordinate
// order; the PanMAN consensus/pseudo-root is then built by progressively
// aligning each sequence's occurrences to a growing reference backbone.
struct MAFBlockGraph {
    // group id -> gapped consensus shared by all occurrences of the block
    std::vector<std::string> consensusByGroup;
    // One occurrence of a block in a sequence (one MAF row).
    struct Occ {
        int group;            // block group id (index into consensusByGroup)
        std::string row;      // aligned, gapped row text
        bool strand;          // true = '+', false = '-'
        int64_t start;        // MAF start coordinate (for ordering)
    };
    // sequence name -> its occurrences sorted by MAF start coordinate
    std::unordered_map<std::string, std::vector<Occ>> seqOccs;
    // sequence names in first-seen order (deterministic reference selection)
    std::vector<std::string> seqOrder;
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

struct Chr {
    int64_t chrIdx=0;
    std::string chrName;
    std::vector< int64_t > blockIds;

    void loadFromProtobuf(panman::ChrList::Reader chrList) {
        chrIdx = chrList.getChrIdx();
        chrName = chrList.getChrName();
        for(auto blockId: chrList.getBlockIds()) {
            blockIds.push_back(blockId);
        }
    }

    Chr(int64_t chrIdxInput, std::vector< int64_t >& blockIdsInput){
        chrIdx=chrIdxInput;
        for(int i=0; i<blockIdsInput.size(); i++){
            blockIds.push_back(blockIdsInput[i]);
        }
    }

    Chr(){};

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

    void assignMutationsToNodes(Node* root, size_t& currentIndex,
                                std::vector< panmanOld::node >& nodes);

    // Get the total number of mutations of given type
    size_t getTotalParsimonyParallel(NucMutationType nucMutType,
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
                          bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false, int chrID = -1);
    
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
    Alphabet alphabet = Alphabet::DNA;
    std::vector< Chr > chrList;
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
    Tree(const panmanOld::tree& mainTree);
    Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
    Tree(std::ifstream& fin, std::ifstream& secondFin,
         FILE_TYPE ftype = FILE_TYPE::GFA, std::string reference = "", std::string refSeqFile = "",
         Alphabet alphabetInput = getActiveAlphabet());

    // Copy blocks from current tree into new tree which is rooted at one of the internal
    // nodes of the current tree. Used in split for PanMAN
    Tree(Node* newRoot, const std::vector< Block >& b, const std::vector< GapList >& g,
         std::unordered_map< std::string, int >& cs,
         std::unordered_map< std::string, int >& ri,
         std::unordered_map< std::string, bool >& si,
         const BlockGapList& bgl);
    

    void protoMATToTree(const panman::Tree::Reader& mainTree);
    void protoMATToTree(const panmanOld::tree& mainTree);
    void test();

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
    void printFASTAUltraFast(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false, int chrID = -1);
    void printSingleNode(std::ostream& fout, const sequence_t& sequence,
                                         const blockExists_t& blockExists, const blockStrand_t& blockStrand,
                                         std::string nodeIdentifier, std::tuple< int, int, int, int > &panMATStart, std::tuple< int, int, int, int > &panMATEnd);
    void printFASTAParallel(std::ostream& fout, bool aligned = false);
    void printMAF(std::ostream& fout);

    void printMAFNew(std::ostream& fout);
    void generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout);
    // Parse a MAF into a block graph (see MAFBlockGraph). Each alignment record
    // becomes a column whose rows are the homologous segments; within-record
    // duplications are split into parallel instance columns.
    void readMAF(std::ifstream& fin, MAFBlockGraph& graph);
    void printVCFParallel(std::string reference, std::ostream& fout);
    void printVCFParallel(panmanUtils::Node* node, std::string& fileName);
    void extractAminoAcidTranslations(std::ostream& fout, int64_t start, int64_t end);
    void printConsensus(std::ostream& fout);
    void printPseduoRoot(std::ostream& fout);
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
    void printMutationsNew(std::ostream& fout, std::string& referenceString);
    void printMutationsNew(std::ostream& fout, std::vector<std::string>& nodes, std::string& referenceString);
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
    void convertToGFAEfficient(std::ostream& fout);
    void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);
    void getNodesPreorder(panmanUtils::Node* root, capnp::List<panman::Node>::Builder& nodesBuilder, size_t& nodeIndex);
    size_t getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition,
                               int nucGapPosition);

    // Transforms tree such that given node becomes child of new root
    void transform(Node* node);
    void reroot(std::string sequenceName);
    void ratioTest(std::ostream& fout);

    

};

// Directed edge in a PanMAN Network (in addition to tree edges stored via Node::parent/children)
struct NetworkEdge {
    Node* parent;
    Node* child;
};

// Phylogenetic network built from a sequence of trees related by single-SPR moves.
// Tree 1 provides the base topology (Node::parent/children); each subsequent SPR move
// adds a single directed "network edge" from the new grafting parent to the moved clade
// root. Clade identity across trees is determined by leaf-set equality.
class Network {
  public:
    // Root of the base (tree 1) topology.
    Node* root = nullptr;

    // All nodes owned by this network (keyed by identifier).
    std::unordered_map< std::string, Node* > allNodes;

    // Auxiliary, non-tree parent->child edges introduced by SPR moves.
    std::vector< NetworkEdge > networkEdges;

    // Depth of each node in tree 1 (root has depth 0). Nodes that were added later
    // (graft points whose clade did not exist in tree 1) are absent from this map.
    std::unordered_map< std::string, size_t > tree1Depth;

    // Map from clade (leaf set) -> Node. Used to identify "the same" internal node
    // across successive trees.
    std::map< std::set< std::string >, Node* > leafsetToNode;

    // Build the network by reading a text file containing a sequence of Newick trees,
    // each optionally preceded by a header line of the form:
    //     Step N: SPR move MCC {a,b,c} -> sibling of {d,e} ...
    // The first tree provides the base topology; each subsequent tree contributes
    // one SPR-derived network edge.
    explicit Network(std::istream& fin);
    ~Network();

    // Number of tree-1 leaves.
    size_t numLeaves() const;

    // Summary: #nodes, #tree edges, #network edges, #leaves.
    void printSummary(std::ostream& out) const;

    // Print the base tree in Newick format plus an extra section listing the
    // network edges as "parent -> child" pairs (one per line).
    void printNetwork(std::ostream& out) const;

  private:
    // Lightweight parse-only node used while parsing each Newick tree.
    struct PNode {
        std::string label;                    // empty for internal nodes
        PNode* parent = nullptr;
        std::vector< PNode* > children;
        std::set< std::string > leafSet;
    };

    // Monotonic counter used to synthesize internal node identifiers.
    size_t m_nextInternalId = 0;
    std::string newInternalId();

    // ---- Newick parsing helpers ----
    PNode* parseNewick(const std::string& newick);
    PNode* parseNewickHelper(const std::string& s, size_t& pos);
    void computeLeafSets(PNode* node);
    void freePNode(PNode* node);
    PNode* findByLeafSet(PNode* root, const std::set< std::string >& ls);

    // ---- Network construction helpers ----
    // Materialize tree 1 from its parse tree, creating Node* objects with
    // parent/children set, and populating `leafsetToNode`, `allNodes`, `tree1Depth`.
    void buildTree1(PNode* pRoot);

    // Return the network node representing this leaf set, creating a new one if
    // necessary. Single-element sets yield leaf nodes; otherwise, an internal
    // network node is created.
    Node* getOrCreateNode(const std::set< std::string >& ls);

    // Parse the "MCC {...}" portion of a step-header line into a leaf-set.
    std::set< std::string > parseMCCFromHeader(const std::string& header);

    // Add a directed network edge parent->child, guarding against duplicates,
    // redundancy with tree-1 edges, and cycles. In case of a cycle, resolve by
    // keeping the node that is deeper in tree 1 as the child and discarding the
    // edge whose direction contradicts that.
    void addNetworkEdgeSafe(Node* par, Node* ch);

    // Is `anc` an ancestor of `desc` in the current DAG (tree + network edges)?
    bool isAncestor(Node* anc, Node* desc) const;

    // ---- Newick printing ----
    void printNewick(std::ostream& out, Node* node) const;
};

// Coordinates of a parent sequence segment contributing to a complex mutation
struct ParentContribution {
    size_t treeIndex;
    std::string sequenceId;

    int32_t primaryBlockIdStart;
    int32_t secondaryBlockIdStart;
    int32_t nucPositionStart;
    int32_t nucGapPositionStart;

    int32_t primaryBlockIdEnd;
    int32_t secondaryBlockIdEnd;
    int32_t nucPositionEnd;
    int32_t nucGapPositionEnd;

    ParentContribution() = default;

    ParentContribution(size_t tIndex, std::string sId,
                       std::tuple< int,int,int,int > startCoords,
                       std::tuple< int,int,int,int > endCoords) {
        treeIndex = tIndex;
        sequenceId = std::move(sId);
        primaryBlockIdStart = std::get<0>(startCoords);
        secondaryBlockIdStart = std::get<1>(startCoords);
        nucPositionStart = std::get<2>(startCoords);
        nucGapPositionStart = std::get<3>(startCoords);
        primaryBlockIdEnd = std::get<0>(endCoords);
        secondaryBlockIdEnd = std::get<1>(endCoords);
        nucPositionEnd = std::get<2>(endCoords);
        nucGapPositionEnd = std::get<3>(endCoords);
    }

    static ParentContribution fromCapnProto(panman::ParentContribution::Reader pc) {
        ParentContribution parent;
        parent.treeIndex = pc.getTreeIndex();
        parent.sequenceId = pc.getSequenceId();

        parent.primaryBlockIdStart = (pc.getBlockIdStart() >> 32);
        parent.secondaryBlockIdStart = (pc.getBlockGapExistStart()?
                                        (pc.getBlockIdStart() & 0xFFFFFFFF): -1);
        parent.nucPositionStart = pc.getNucPositionStart();
        parent.nucGapPositionStart = (pc.getNucGapExistStart()?
                                      pc.getNucGapPositionStart() : -1);

        parent.primaryBlockIdEnd = (pc.getBlockIdEnd() >> 32);
        parent.secondaryBlockIdEnd = (pc.getBlockGapExistEnd()?
                                      (pc.getBlockIdEnd() & 0xFFFFFFFF) : -1);
        parent.nucPositionEnd = pc.getNucPositionEnd();
        parent.nucGapPositionEnd = (pc.getNucGapExistEnd()?
                                    pc.getNucGapPositionEnd() : -1);
        return parent;
    }

    static ParentContribution fromProtobuf(const panmanOld::parentContribution& pc) {
        ParentContribution parent;
        parent.treeIndex = pc.treeindex();
        parent.sequenceId = pc.sequenceid();

        parent.primaryBlockIdStart = (pc.blockidstart() >> 32);
        parent.secondaryBlockIdStart = (pc.blockgapexiststart()?
                                      (pc.blockidstart() & 0xFFFFFFFF) : -1);
        parent.nucPositionStart = pc.nucpositionstart();
        parent.nucGapPositionStart = (pc.nucgapexiststart()?
                                      pc.nucgappositionstart() : -1);

        parent.primaryBlockIdEnd = (pc.blockidend() >> 32);
        parent.secondaryBlockIdEnd = (pc.blockgapexistend()?
                                      (pc.blockidend() & 0xFFFFFFFF) : -1);
        parent.nucPositionEnd = pc.nucpositionend();
        parent.nucGapPositionEnd = (pc.nucgapexistend()?
                                    pc.nucgappositionend() : -1);
        return parent;
    }

    void toCapnProto(panman::ParentContribution::Builder& pc) const {
        pc.setTreeIndex(treeIndex);
        pc.setSequenceId(sequenceId);

        if(secondaryBlockIdStart != -1) {
            pc.setBlockGapExistStart(true);
            pc.setBlockIdStart(((int64_t)primaryBlockIdStart << 32) + secondaryBlockIdStart);
        } else {
            pc.setBlockGapExistStart(false);
            pc.setBlockIdStart(((int64_t)primaryBlockIdStart << 32));
        }
        pc.setNucPositionStart(nucPositionStart);
        if(nucGapPositionStart != -1) {
            pc.setNucGapExistStart(true);
            pc.setNucGapPositionStart(nucGapPositionStart);
        }

        if(secondaryBlockIdEnd != -1) {
            pc.setBlockGapExistEnd(true);
            pc.setBlockIdEnd(((int64_t)primaryBlockIdEnd << 32) + secondaryBlockIdEnd);
        } else {
            pc.setBlockGapExistEnd(false);
            pc.setBlockIdEnd(((int64_t)primaryBlockIdEnd << 32));
        }
        pc.setNucPositionEnd(nucPositionEnd);
        if(nucGapPositionEnd != -1) {
            pc.setNucGapExistEnd(true);
            pc.setNucGapPositionEnd(nucGapPositionEnd);
        }
    }

    panmanOld::parentContribution toProtobuf() const {
        panmanOld::parentContribution pc;
        pc.set_treeindex(treeIndex);
        pc.set_sequenceid(sequenceId);

        if(secondaryBlockIdStart != -1) {
            pc.set_blockgapexiststart(true);
            pc.set_blockidstart(((int64_t)primaryBlockIdStart << 32) + secondaryBlockIdStart);
        } else {
            pc.set_blockgapexiststart(false);
            pc.set_blockidstart(((int64_t)primaryBlockIdStart << 32));
        }
        pc.set_nucpositionstart(nucPositionStart);
        if(nucGapPositionStart != -1) {
            pc.set_nucgapexiststart(true);
            pc.set_nucgappositionstart(nucGapPositionStart);
        }

        if(secondaryBlockIdEnd != -1) {
            pc.set_blockgapexistend(true);
            pc.set_blockidend(((int64_t)primaryBlockIdEnd << 32) + secondaryBlockIdEnd);
        } else {
            pc.set_blockgapexistend(false);
            pc.set_blockidend(((int64_t)primaryBlockIdEnd << 32));
        }
        pc.set_nucpositionend(nucPositionEnd);
        if(nucGapPositionEnd != -1) {
            pc.set_nucgapexistend(true);
            pc.set_nucgappositionend(nucGapPositionEnd);
        }
        return pc;
    }
};

// Represents complex mutations like Horizontal Gene Transfer or Recombinations
struct ComplexMutation {
    char mutationType;
    size_t childTreeIndex;
    std::string childSequenceId;
    std::vector< ParentContribution > parents;

    ComplexMutation(char mutType, size_t childTreeIdx, std::string childSeqId,
                    std::vector< ParentContribution > parentContribs) {
        mutationType = mutType;
        childTreeIndex = childTreeIdx;
        childSequenceId = std::move(childSeqId);
        parents = std::move(parentContribs);
    }

    static ParentContribution legacyParentFromCapnProto(
        size_t treeIndex, const std::string& sequenceId,
        int64_t blockIdStart, bool blockGapExistStart,
        int32_t nucPositionStart, int32_t nucGapPositionStart, bool nucGapExistStart,
        int64_t blockIdEnd, bool blockGapExistEnd,
        int32_t nucPositionEnd, int32_t nucGapPositionEnd, bool nucGapExistEnd) {
        ParentContribution parent;
        parent.treeIndex = treeIndex;
        parent.sequenceId = sequenceId;
        parent.primaryBlockIdStart = (blockIdStart >> 32);
        parent.secondaryBlockIdStart = (blockGapExistStart?
                                        (blockIdStart & 0xFFFFFFFF) : -1);
        parent.nucPositionStart = nucPositionStart;
        parent.nucGapPositionStart = (nucGapExistStart? nucGapPositionStart : -1);
        parent.primaryBlockIdEnd = (blockIdEnd >> 32);
        parent.secondaryBlockIdEnd = (blockGapExistEnd?
                                      (blockIdEnd & 0xFFFFFFFF) : -1);
        parent.nucPositionEnd = nucPositionEnd;
        parent.nucGapPositionEnd = (nucGapExistEnd? nucGapPositionEnd : -1);
        return parent;
    }

    static void readLegacyInlineParents(panman::ComplexMutation::Reader cm,
                                        std::vector< ParentContribution >& out) {
        out.push_back(legacyParentFromCapnProto(
            cm.getTreeIndex1(), cm.getSequenceId1(),
            cm.getBlockIdStart1(), cm.getBlockGapExistStart1(),
            cm.getNucPositionStart1(), cm.getNucGapPositionStart1(), cm.getNucGapExistStart1(),
            cm.getBlockIdEnd1(), cm.getBlockGapExistEnd1(),
            cm.getNucPositionEnd1(), cm.getNucGapPositionEnd1(), cm.getNucGapExistEnd1()));
        out.push_back(legacyParentFromCapnProto(
            cm.getTreeIndex2(), cm.getSequenceId2(),
            cm.getBlockIdStart2(), cm.getBlockGapExistStart2(),
            cm.getNucPositionStart2(), cm.getNucGapPositionStart2(), cm.getNucGapExistStart2(),
            cm.getBlockIdEnd2(), cm.getBlockGapExistEnd2(),
            cm.getNucPositionEnd2(), cm.getNucGapPositionEnd2(), cm.getNucGapExistEnd2()));
    }

    static void readLegacyInlineParents(const panmanOld::complexMutation& cm,
                                        std::vector< ParentContribution >& out) {
        out.push_back(legacyParentFromCapnProto(
            cm.treeindex1(), cm.sequenceid1(),
            cm.blockidstart1(), cm.blockgapexiststart1(),
            cm.nucpositionstart1(), cm.nucgappositionstart1(), cm.nucgapexiststart1(),
            cm.blockidend1(), cm.blockgapexistend1(),
            cm.nucpositionend1(), cm.nucgappositionend1(), cm.nucgapexistend1()));
        out.push_back(legacyParentFromCapnProto(
            cm.treeindex2(), cm.sequenceid2(),
            cm.blockidstart2(), cm.blockgapexiststart2(),
            cm.nucpositionstart2(), cm.nucgappositionstart2(), cm.nucgapexiststart2(),
            cm.blockidend2(), cm.blockgapexistend2(),
            cm.nucpositionend2(), cm.nucgappositionend2(), cm.nucgapexistend2()));
    }

    ComplexMutation(panman::ComplexMutation::Reader cm) {
        mutationType = (cm.getMutationType()? 'H': 'R');
        childTreeIndex = cm.getTreeIndex3();
        childSequenceId = cm.getSequenceId3();

        if(cm.hasParents()) {
            try {
                for(auto parentReader: cm.getParents()) {
                    parents.push_back(ParentContribution::fromCapnProto(parentReader));
                }
            } catch (const kj::Exception&) {
                parents.clear();
                readLegacyInlineParents(cm, parents);
            }
        } else {
            readLegacyInlineParents(cm, parents);
        }
    }

    void toCapnProto(panman::ComplexMutation::Builder& cm) const {
        cm.setMutationType(mutationType == 'H');
        cm.setTreeIndex3(childTreeIndex);
        cm.setSequenceId3(childSequenceId);

        auto parentsBuilder = cm.initParents(parents.size());
        for(size_t i = 0; i < parents.size(); i++) {
            auto parentBuilder = parentsBuilder[i];
            parents[i].toCapnProto(parentBuilder);
        }
    }

    ComplexMutation(panmanOld::complexMutation cm) {
        mutationType = (cm.mutationtype()? 'H': 'R');
        childTreeIndex = cm.treeindex3();
        childSequenceId = cm.sequenceid3();

        if(cm.parents_size() > 0) {
            for(int i = 0; i < cm.parents_size(); i++) {
                parents.push_back(ParentContribution::fromProtobuf(cm.parents(i)));
            }
        } else {
            readLegacyInlineParents(cm, parents);
        }
    }

    panmanOld::complexMutation toProtobuf() const {
        panmanOld::complexMutation cm;
        cm.set_mutationtype(mutationType == 'H');
        cm.set_treeindex3(childTreeIndex);
        cm.set_sequenceid3(childSequenceId);

        for(const auto& parent: parents) {
            *cm.add_parents() = parent.toProtobuf();
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

    TreeGroup(std::istream& fin, bool isOld = false);
    // List of PanMAT files and a file with all the complex mutations relating these files
    TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile);
    TreeGroup(std::vector< Tree* >& t);
    TreeGroup(std::vector< Tree* >& tg, std::ifstream& mutationFile);

    TreeGroup* subnetworkExtract(std::unordered_map< int, std::vector< std::string > >& nodeIds);

    void printFASTA(std::ofstream& fout, bool rootSeq = false);
    void writeToFile(kj::std::StdOutputStream& fout);
    void printComplexMutations(std::ostream& fout);
};

// Build a PanMAN (TreeGroup) from TreeKnit output. Inputs:
//   treeknitDir  - directory of chrXX_window_N_N+1 subdirectories, each containing
//                  MCCs.json (exactly 2 MCCs -> one SPR move) and per-window
//                  <chrXX_A_B>_resolved.nwk files.
//   alignmentDir - directory of per-window MSA files named chrXX_A_B.fa.
//   refFile      - reference FASTA (.fa / .fa.gz / .fa.xz), one record per chromosome.
// Each window becomes a PanMAT (one block per window, uncovered reference gaps merged
// into the left neighbour block). Adjacent windows are linked by one SPR-derived
// recombination complex mutation (smaller MCC = moved clade; breakpoints = MSA
// coordinates of the child window block; parents identified by node id).
TreeGroup* buildNetworkFromTreeknit(const std::string& treeknitDir,
                                    const std::string& alignmentDir,
                                    const std::string& refFile);

};