#ifndef PANGENOME_MAT_HPP
#define PANGENOME_MAT_HPP

#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>

#include <json/json.h>
#include "mutation_annotation_test_proto3_optional_new.pb.h"
#include "kseq.h"
#include "syncmer.hpp"
#include "statsgenotype.hpp"

#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"

typedef std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >,
            std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence_t;
// Individual block
typedef std::vector< std::pair< char, std::vector< char > > > block_t;

typedef  std::vector< std::pair< bool, std::vector< bool > > > blockExists_t;
// Forward or reverse strand
typedef  std::vector< std::pair< bool, std::vector< bool > > > blockStrand_t;

namespace PangenomeMAT {

    static inline void printError(std::string e) {
        std::cout << "\033[1;31m" << "Error: " << "\033[0m" << e << std::endl;
    }

    // Get nucleotide character from 4-bit PanMAT code
    char getNucleotideFromCode(int code);

    // Get 4-bit PanMAT code from nucleotide character
    char getCodeFromNucleotide(char nuc);

    // Get complement Nucleotide character of given nucleotide character. Used to compute reverse
    // complement
    char getComplementCharacter(char nuc);

    // Given a sequence and block presence/strand information, print the sequence in FASTA format
    // where each line has length lineSize
    void printSequenceLines(const sequence_t& sequence,
        const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize,
        bool aligned, std::ofstream& fout, int offset = 0, bool debug = false);

    // Remove '-' character from sequence string
    std::string stripGaps(std::string sequenceString);
    std::string getDate();
    std::string stripString(std::string s);

    enum FILE_TYPE {
        PANMAT = 0,
        GFA = 1,
        PANGRAPH=2,
        MSA = 3,
        MSA_OPTIMIZE = 4,
        // FASTA = 5
    };

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

    void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);

    // Struct for representing Nucleotide Mutation
    struct NucMut {
        // For creating SNP mutation
        NucMut( const std::tuple< int, int, int, int, int, int >& mutationInfo ) {
            // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
            primaryBlockId = std::get<0>(mutationInfo);
            secondaryBlockId = std::get<1>(mutationInfo);
            nucPosition = std::get<2>(mutationInfo);
            nucGapPosition = std::get<3>(mutationInfo);
            mutInfo = std::get<4>(mutationInfo) + (1 << 4);
            nucs = (std::get<5>(mutationInfo) << 20);
        }

        // For creating non-SNP mutations from SNP mutations at consecutive positions
        NucMut(const std::vector< std::tuple< int, int, int, int, int, int > >& mutationArray,
            int start, int end) {
            primaryBlockId = std::get<0>(mutationArray[start]);
            secondaryBlockId = std::get<1>(mutationArray[start]);

            mutInfo = ((end - start) << 4);
            // type
            switch(std::get<4>(mutationArray[start])) {
                case PangenomeMAT::NucMutationType::NSNPS:
                    mutInfo += PangenomeMAT::NucMutationType::NS;
                    break;
                case PangenomeMAT::NucMutationType::NSNPI:
                    mutInfo += PangenomeMAT::NucMutationType::NI;
                    break;
                case PangenomeMAT::NucMutationType::NSNPD:
                    mutInfo += PangenomeMAT::NucMutationType::ND;
                    break;
                case PangenomeMAT::NucMutationType::NS:
                    mutInfo += PangenomeMAT::NucMutationType::NS;
                    break;
                case PangenomeMAT::NucMutationType::NI:
                    mutInfo += PangenomeMAT::NucMutationType::NI;
                    break;
                case PangenomeMAT::NucMutationType::ND:
                    mutInfo += PangenomeMAT::NucMutationType::ND;
                    break;
            }

            nucPosition = std::get<2>(mutationArray[start]);
            nucGapPosition = std::get<3>(mutationArray[start]);

            nucs = 0;
            for(int i = start; i < end; i++) {
                nucs += (std::get<5>(mutationArray[i]) << (4*(5-(i - start))));
            }
        }

        NucMut(PanMAT::nucMut mutation, int64_t blockId, bool blockGapExist) {
            nucPosition = mutation.nucposition();
            primaryBlockId = (blockId >> 32);
            mutInfo = (mutation.mutinfo() & 0xFF);
            nucs = (mutation.mutinfo() >> 8);
            nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

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

        int32_t nucPosition;
        int32_t nucGapPosition;
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        uint8_t mutInfo;
        uint32_t nucs;
    };

    // Struct for representing Block Mutations
    struct BlockMut {
        void loadFromProtobuf(PanMAT::mutation mutation) {
            primaryBlockId = (mutation.blockid() >> 32);
            if(mutation.blockgapexist()) {
                secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
            } else {
                secondaryBlockId = -1;
            }
            blockMutInfo = mutation.blockmutinfo();
            // If the mutation is a block inversion or not. Inversion is marked by
            // `blockMutInfo = deletion` and `inversion = true`
            inversion = mutation.blockinversion();
        }

        BlockMut(size_t blockId, bool type, int secondaryBId = -1) {
            primaryBlockId = blockId;
            secondaryBlockId = secondaryBId;
            blockMutInfo = type;

            // @TODO (Harsh): Update GFA interoperability to account for block strands. Remove this
            // constructor and just use the previous constructor
            inversion = false;
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

        int32_t primaryBlockId;
        int32_t secondaryBlockId;

        // If mutation is an insertion or deletion - Strand inversions are marked by
        // `blockMutInfo=false`, but they are not deletions
        bool  blockMutInfo;

        // Whether the block is being inverted or not
        bool inversion;
    };

    // List of default blocks in the global coordinate system of the PanMAT
    struct Block {
        Block(size_t primaryBlockId, std::string seq);
        // seq is a compressed form of the sequence where each nucleotide is stored in 4 bytes
        Block(int32_t primaryBlockId, int32_t secondaryBlockId, const std::vector< uint32_t >& seq);

        int32_t primaryBlockId;
        int32_t secondaryBlockId;

        std::vector< uint32_t > consensusSeq;
        std::string chromosomeName;
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
        Node(std::string id, float len);
        Node(std::string id, Node* par, float len);

        float branchLength;

        size_t level;

        std::string identifier;
        Node* parent;
        std::vector< Node* > children;

        std::vector< NucMut > nucMutation;
        std::vector< BlockMut > blockMutation;

        std::vector< std::string > annotations;
    };

    class Pangraph{
    private:
        void topologicalSortHelper(size_t nodeId, std::vector< size_t >& topoArray,
            std::vector< bool >& visited);
    public:
        Pangraph(Json::Value& pangraphData);
        std::vector< size_t > getTopologicalSort();
        std::unordered_map< std::string,std::vector< int > >
            getAlignedSequences(const std::vector< size_t >& topoArray);
        // Patch to incorporate strands
        std::unordered_map< std::string,std::vector< int > >
            getAlignedStrandSequences(const std::vector< size_t >& topoArray);

        // Graph adjacency list
        size_t numNodes;
        std::vector< std::vector< size_t > > adj;
        std::unordered_map< std::string, std::vector< size_t > > intSequences;
        std::vector< size_t > topoSortedIntSequences;

        std::unordered_map< std::string, std::vector< std::string > > paths;
        
        // Added as a patch to incorporate strands
        std::unordered_map< std::string, std::vector< int > > strandPaths;

        // Represents the "number" parameter of each block
        std::unordered_map< std::string, std::vector< size_t > > blockNumbers;

        // Circular offset of given sequence. Zero if not circular
        std::unordered_map< std::string, int > circularSequences;

        // If a sequence was rotated, which block was it rotated by
        std::unordered_map< std::string, int > rotationIndexes;

        // Specifies whether sequence is inverted by the rotation algorithm or not
        std::unordered_map< std::string, bool > sequenceInverted;

        std::unordered_map< std::string, std::string > stringIdToConsensusSeq;

        // block identifier to list of gaps - < position, size > pairs
        std::unordered_map< std::string,
            std::vector< std::pair< size_t, size_t > > > stringIdToGaps;

        // Mapping from new integer ID of block to old string ID. Duplicated blocks have same string
        // ID
        std::unordered_map< size_t, std::string > intIdToStringId;

        // List of substitutions for each block
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
            tbb::concurrent_unordered_map< size_t,
            std::vector< std::pair< size_t, std::string > > > > > substitutions;
        // List of insertions for each block
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
            tbb::concurrent_unordered_map< size_t, std::vector< std::tuple< size_t, size_t,
            std::string > > > > > insertions;
        // List of deletions for each block
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
            tbb::concurrent_unordered_map< size_t,
            std::vector< std::pair< size_t, size_t > > > > > deletions;
    };

    class GFAGraph {
    private:
        bool checkForCyclesHelper(size_t nodeId, std::vector< int >& color);
        void topologicalSortHelper(size_t nodeId, std::vector< size_t >& topoArray, std::vector< bool >& visited);
    public:
        size_t numNodes;
        // Graph adjacency list
        std::vector< std::vector< size_t > > adj;

        // Names of the sequences
        std::vector< std::string > pathIds;

        // The sequences themselves where each node has an integer ID
        std::vector< std::vector< size_t > > intSequences;

        // Raw sequence corresponding to each node
        std::vector< std::string > intNodeToSequence;

        GFAGraph(const std::vector< std::string >& pathNames, const std::vector< std::vector< std::string > >& sequences, std::map< std::string, std::string >& nodes);

        bool pathExists(size_t nId1, size_t nId2, std::vector< bool >& visited);
        bool checkForCycles();
        std::vector< size_t > getTopologicalSort();
        std::vector< std::vector< int64_t > > getAlignedSequences(const std::vector< size_t >& topoArray);
    };

    class Tree {
        private:
            Node* createTreeFromNewickString(std::string newick);
            
            // In the proto file, nodes are stored in preorder. Once the tree has been generated in
            // memory, assign mutations from the proto file to the tree nodes using preorder
            // traversal
            void assignMutationsToNodes(Node* root, size_t& currentIndex,
                std::vector< PanMAT::node >& nodes);

            // Get the total number of mutations of given type
            int getTotalParsimonyParallel(NucMutationType nucMutType,
                BlockMutationType blockMutType = NONE);

            // Tree traversal for FASTA writer
            void printFASTAHelper(PangenomeMAT::Node* root, sequence_t& sequence,
                blockExists_t& blockExists, blockStrand_t& blockStrand, std::ofstream& fout,
                bool aligned = false);

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
            void compressTreeParallel(Node* node, size_t level);

            // Used in rerooting
            void dfsExpansion(Node* node, std::vector< Node* >& vec);
            Node* transformHelper(Node* node);
            void adjustLevels(Node* node);

            // Read mutation matrices from file. Infer mutation matrices from tree if no input file provided.
            void fillMutationMats(statsgenotype::mutationMatrices& mutmat, Node* node, std::ifstream* min = nullptr);
            // Used to get substitutions for inferring mutation matrices
            std::pair<std::string, std::string> get_substitution(const std::string& nid, const PangenomeMAT::NucMut& nucmut);

            std::string newInternalNodeId() {
                return "node_" + std::to_string(++m_currInternalNode);
            }

            size_t m_currInternalNode{ 0 };
            size_t m_numLeaves{ 0 };
            size_t m_maxDepth{ 0 };
            float m_meanDepth{ 0 };

            std::unordered_map<std::string, std::vector< std::string > > annotationsToNodes;
        public:
            Tree(const PanMAT::tree& mainTree);
            Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
            Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype = FILE_TYPE::GFA);

            // Copy blocks from current tree into new tree which is rooted at one of the internal nodes of the current tree. Used in split for PanMAN
            Tree(Node* newRoot, const std::vector< Block >& b, const std::vector< GapList >& g, const std::unordered_map< std::string, int >& c, const BlockGapList& bgl);

            void protoMATToTree(const PanMAT::tree& mainTree);

            int nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states);
            void nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState, int defaultState = (1<<28));
            void nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT::NucMutationType, char > >& mutations, int parentState);

            // To account for strand
            int blockFitchForwardPassNew(Node* node, std::unordered_map< std::string, int >& states);
            void blockFitchBackwardPassNew(Node* node, std::unordered_map< std::string, int >& states, int parentState,  int defaultValue = (1<<28));
            void blockFitchAssignMutationsNew(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT::BlockMutationType, bool > >& mutations, int parentState);

            int blockFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states);
            // The defaultValue parameter is used in rerooting.
            void blockFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState,  int defaultValue = 2);
            void blockFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, bool >& mutations, int parentState);

            void printSummary();
            void printBfs(Node* node = nullptr);
            void printFASTA(std::ofstream& fout, bool aligned = false);
            void printFASTAParallel(std::ofstream& fout, bool aligned = false);

            void printMAF(std::ofstream& fout);
            void generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout);

            void printVCFParallel(std::string reference, std::ofstream& fout);

            Node* subtreeExtractParallel(std::vector< std::string > nodeIds);
            void writeToFile(std::ostream& fout, Node* node = nullptr);
            std::string getNewickString(Node* node);
            std::string getStringFromReference(std::string reference, bool aligned = true, bool incorporateInversions=true);
            void getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists, blockStrand_t& blockStrand, std::string reference, int* rotIndex = nullptr);
            void getBlockSequenceFromReference(block_t& sequence, bool& blockExists, bool& blockStrand, std::string reference, int64_t primaryBlockId, int64_t secondaryBlockId);

            // Split file provided as input.
            std::pair< Tree, Tree > splitByComplexMutations(const std::string& nodeId3);
            std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > > globalCoordinates;

            // get unaligned global coordinate
            int32_t getUnalignedGlobalCoordinate(int32_t primaryBlockId, int32_t secondaryBlockId, int32_t pos, int32_t gapPos, const sequence_t& sequence, const blockExists_t& blockExists);
            std::tuple< int, int, int, int > globalCoordinateToBlockCoordinate(int64_t globalCoordinate, const sequence_t& sequence, const blockExists_t& blockExists, const blockStrand_t& blockStrand);

            std::string getSequenceFromVCF(std::string sequenceId, std::ifstream& fin);
            bool verifyVCFFile(std::ifstream& fin);
            void vcfToFASTA(std::ifstream& fin, std::ofstream& fout);
            void annotate(std::ifstream& fin);
            std::vector< std::string > searchByAnnotation(std::string annotation);
            void convertToGFA(std::ofstream& fout);
            void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);
            void getNodesPreorder(PangenomeMAT::Node* root, PanMAT::tree& treeToWrite);
            void setupGlobalCoordinates();
            size_t getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition, int nucGapPosition);
            
            pair< vector<statsgenotype::variationSite>, pair<size_t, size_t> > getVariantSites(std::ifstream& fin, std::ifstream* min = nullptr);
            int printSamplePlacementVCF(std::ifstream& fin, std::ifstream* min = nullptr);
            
            // Transforms tree such that given node becomes child of new root
            void transform(Node* node);
            void reroot(std::string sequenceName);

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

    };

    struct ComplexMutation{
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

        ComplexMutation(char mutType, int tIndex1, int tIndex2, int tIndex3, std::string sId1, std::string sId2, std::string sId3, std::tuple< int,int,int,int > t1, std::tuple< int,int,int,int > t2, std::tuple< int,int,int,int > t3, std::tuple< int,int,int,int > t4) {
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

        ComplexMutation(PanMAT::complexMutation cm) {
            mutationType = (cm.mutationtype()? 'H': 'R');
            treeIndex1 = cm.treeindex1();
            treeIndex2 = cm.treeindex2();
            treeIndex3 = cm.treeindex3();
            sequenceId1 = cm.sequenceid1();
            sequenceId2 = cm.sequenceid2();
            sequenceId3 = cm.sequenceid3();

            primaryBlockIdStart1 = (cm.blockidstart1() >> 32);
            secondaryBlockIdStart1 = (cm.blockgapexiststart1()? (cm.blockidstart1()&(0xFFFFFFFF)): -1);
            nucPositionStart1 = cm.nucpositionstart1();
            nucGapPositionStart1 = (cm.nucgapexiststart1()? (cm.nucgappositionstart1()) : -1);

            primaryBlockIdStart2 = (cm.blockidstart2() >> 32);
            secondaryBlockIdStart2 = (cm.blockgapexiststart2()? (cm.blockidstart2()&(0xFFFFFFFF)): -1);
            nucPositionStart2 = cm.nucpositionstart2();
            nucGapPositionStart2 = (cm.nucgapexiststart2()? (cm.nucgappositionstart2()) : -1);

            primaryBlockIdEnd1 = (cm.blockidend1() >> 32);
            secondaryBlockIdEnd1 = (cm.blockgapexistend1()? (cm.blockidend1()&(0xFFFFFFFF)): -1);
            nucPositionEnd1 = cm.nucpositionend1();
            nucGapPositionEnd1 = (cm.nucgapexistend1()? (cm.nucgappositionend1()) : -1);

            primaryBlockIdEnd2 = (cm.blockidend2() >> 32);
            secondaryBlockIdEnd2 = (cm.blockgapexistend2()? (cm.blockidend2()&(0xFFFFFFFF)): -1);
            nucPositionEnd2 = cm.nucpositionend2();
            nucGapPositionEnd2 = (cm.nucgapexistend2()? (cm.nucgappositionend2()) : -1);
        }

        PanMAT::complexMutation toProtobuf() {
            PanMAT::complexMutation cm;
            cm.set_mutationtype(mutationType == 'H');
            cm.set_treeindex1(treeIndex1);
            cm.set_treeindex2(treeIndex2);
            cm.set_treeindex3(treeIndex3);
            cm.set_sequenceid1(sequenceId1);
            cm.set_sequenceid2(sequenceId2);
            cm.set_sequenceid3(sequenceId3);

            if(secondaryBlockIdStart1 != -1) {
                cm.set_blockgapexiststart1(true);
                cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32)+secondaryBlockIdStart1);
            } else {
                cm.set_blockgapexiststart1(false);
                cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32));
            }
            cm.set_nucpositionstart1(nucPositionStart1);

            if(nucGapPositionStart1 != -1) {
                cm.set_nucgapexiststart1(true);
                cm.set_nucgappositionstart1(nucGapPositionStart1);
            }

            if(secondaryBlockIdStart2 != -1) {
                cm.set_blockgapexiststart2(true);
                cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32)+secondaryBlockIdStart2);
            } else {
                cm.set_blockgapexiststart2(false);
                cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32));
            }
            cm.set_nucpositionstart2(nucPositionStart2);

            if(nucGapPositionStart2 != -1) {
                cm.set_nucgapexiststart2(true);
                cm.set_nucgappositionstart2(nucGapPositionStart2);
            }

            if(secondaryBlockIdEnd1 != -1) {
                cm.set_blockgapexistend1(true);
                cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32)+secondaryBlockIdEnd1);
            } else {
                cm.set_blockgapexistend1(false);
                cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32));
            }
            cm.set_nucpositionend1(nucPositionEnd1);

            if(nucGapPositionEnd1 != -1) {
                cm.set_nucgapexistend1(true);
                cm.set_nucgappositionend1(nucGapPositionEnd1);
            }

            if(secondaryBlockIdEnd2 != -1) {
                cm.set_blockgapexistend2(true);
                cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32)+secondaryBlockIdEnd2);
            } else {
                cm.set_blockgapexistend2(false);
                cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32));
            }
            cm.set_nucpositionend2(nucPositionEnd2);

            if(nucGapPositionEnd2 != -1) {
                cm.set_nucgapexistend2(true);
                cm.set_nucgappositionend2(nucGapPositionEnd2);
            }

            return cm;
        }

    };

    class TreeGroup {
    public:
        std::vector< Tree > trees;
        std::vector< ComplexMutation > complexMutations;
        
        TreeGroup(std::ifstream& fin);
        TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile);
        TreeGroup(const std::vector< Tree >& t);

        void printFASTA(std::ofstream& fout);
        void writeToFile(std::ofstream& fout);
        void printComplexMutations();
    };

    /* seed indexing */

    struct Counter
    {
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    int32_t count = 0;
    };

    template<typename T1, typename T2> int32_t intersection_size(const T1& s1, const T2& s2)
    {
    Counter c;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
    return c.count;
    }

    template <typename INT, typename T>
    void removeIndices(std::vector<T>& v, std::stack<INT>& rm);
    void discardSyncmers(std::vector<kmer_t> &inSyncmers, const std::vector<std::pair<int32_t, int32_t>>& B, std::string &gappedSequence, std::unordered_map<std::string, kmer_t> &to_insert, std::unordered_map<std::string, bool> &variable_syncmers, seedIndex &index, std::string nid, size_t k);
    void shaveDFS(Node *currNode, std::vector<kmer_t> &currNodeSyncmers, seedIndex &index_orig, std::unordered_map<std::string, std::vector<kmer_t>> &deletions, std::unordered_map<std::string, bool> &invariants);
    void placeDFS(Node *currNode, std::vector<kmer_t> &currNodeSyncmers, std::unordered_map<std::string, bool> &querySyncmers, seedIndex &index, dynamicJaccard dj, std::unordered_map<std::string, float> &scores);

    void writeIndexDFS(Node *currNode, seedIndex &index, std::stringstream &ss, std::vector<kmer_t> &seeds);
    void writeIndex(Node *node, std::ofstream &fout, seedIndex &index);
    void indexSyncmersHelper(Tree *T, Node *root, sequence_t &sequence, blockExists_t &blockExists, blockStrand_t &blockStrand, seedIndex &index, std::vector<kmer_t> &syncmers, std::unordered_map<std::string, int32_t> &counts, std::unordered_map<std::string, bool> &variable_syncmers, size_t k, size_t s);
    void indexSyncmers(Tree *t, std::ofstream& fout, size_t k, size_t s);
    seedIndex shaveIndex(Node *root, seedIndex &index, std::unordered_map<std::string, bool> &invariants, std::vector<kmer_t> &initialSyncmers);
    std::pair<int32_t, int32_t> getRecomputePositions(std::pair<int32_t, int32_t> p, std::string &gappedSequence, int32_t k);
    int32_t alignedEndPos(int32_t pos, int32_t k, std::string &gappedSequence);
    void loadIndex(Node *root, std::ifstream &indexFile, seedIndex &index);
    void placeSample(Tree *T, std::string fastqPath, seedIndex &index, size_t k, size_t s, string& best_match);
    std::set<kmer_t> syncmersFromFastq(std::string fastqPath,  std::vector<read_t> &reads, size_t k, size_t s);
    std::string getCurrentFastaSequence(const sequence_t& sequence, const blockExists_t& blockExists, blockStrand_t& blockStrand, bool aligned);


};

#endif // PANGENOME_MAT_HPP