#ifndef PANGENOME_MAT_NEW_HPP
#define PANGENOME_MAT_NEW_HPP

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
#include "spoa/spoa.hpp"

#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"

typedef std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence_t;
typedef  std::vector< std::pair< bool, std::vector< bool > > > blockExists_t;
// Forward or reverse strand
typedef  std::vector< std::pair< bool, std::vector< bool > > > blockStrand_t;

namespace PangenomeMAT {

    char getNucleotideFromCode(int code);
    char getCodeFromNucleotide(char nuc);
    char getComplementCharacter(char nuc);
    void printSequenceLines(const std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
        const std::vector< std::pair< bool, std::vector< bool > > >& blockExists, blockStrand_t& blockStrand, size_t lineSize, bool aligned, std::ofstream& fout, bool debug = false);
    std::pair< int, int > replaceMutation(std::pair<int,int> oldMutation, std::pair<int, int> newMutation);
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
        NS = 0,
        ND = 1,
        NI = 2,
        NSNPS = 3,
        NSNPI = 4,
        NSNPD = 5,
        NNONE = 2000
    };

    enum BlockMutationType {
        BI = 1,
        BD = 0,
        BIn = 2,
        NONE = 1000
    };

    void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);

    struct NucMut {
        // For creating SNP mutation
        NucMut( const std::tuple< int, int, int, int, int, int >& mutationInfo ){
            // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
            primaryBlockId = std::get<0>(mutationInfo);
            secondaryBlockId = std::get<1>(mutationInfo);
            nucPosition = std::get<2>(mutationInfo);
            nucGapPosition = std::get<3>(mutationInfo);
            mutInfo = std::get<4>(mutationInfo) + (1 << 4);
            nucs = (std::get<5>(mutationInfo) << 20);
        }

        // For creating non-SNP mutations from SNP mutations
        NucMut(const std::vector< std::tuple< int, int, int, int, int, int > >& mutationArray, int start, int end){
            primaryBlockId = std::get<0>(mutationArray[start]);
            secondaryBlockId = std::get<1>(mutationArray[start]);

            mutInfo = ((end - start) << 4);
            // type
            switch(std::get<4>(mutationArray[start])){
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
            for(int i = start; i < end; i++){
                nucs += (std::get<5>(mutationArray[i]) << (4*(5-(i - start))));
            }
        }

        NucMut(MATNew::nucMut mutation, int64_t blockId, bool blockGapExist){
            nucPosition = mutation.nucposition();
            primaryBlockId = (blockId >> 32);
            mutInfo = (mutation.mutinfo() & 0xFF);
            nucs = (mutation.mutinfo() >> 8);
            nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

            if(blockGapExist){
                secondaryBlockId = (blockId & 0xFFFFFFFF);
            } else {
                secondaryBlockId = -1;
            }

            if(mutation.nucgapexist()){
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

    struct BlockMut {

        void loadFromProtobuf(MATNew::mutation mutation){
            primaryBlockId = (mutation.blockid() >> 32);
            if(mutation.blockgapexist()){
                secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
            } else {
                secondaryBlockId = -1;
            }
            blockMutInfo = mutation.blockmutinfo();
            
            inversion = mutation.blockinversion();
        }

        BlockMut(size_t blockId, bool type, int secondaryBId = -1){
            primaryBlockId = blockId;
            secondaryBlockId = secondaryBId;
            blockMutInfo = type;

            // Change
            inversion = false;
        }

        BlockMut(size_t blockId, std::pair< BlockMutationType, bool > type, int secondaryBId = -1){
            primaryBlockId = blockId;
            secondaryBlockId = secondaryBId;
            if(type.first == BlockMutationType::BI){
                blockMutInfo = true;
            } else {
                // blockMutInfo is also set to false in the case of inversions. Generally it indicates deletions
                blockMutInfo = false;
            }

            if(type.second == 1){
                // If type.second == 1 (inversion)
                inversion = true;
            } else {
                inversion = false;
            }
        }

        BlockMut(){
            
        }

        int32_t primaryBlockId;
        int32_t secondaryBlockId;

        // The following two variables are separate because the strand patch was added later

        // If mutation is an insertion or deletion - Strand inversions are marked by blockMutInfo=false, however, they are not deletions
        bool  blockMutInfo;
        
        // Whether the block is being inverted or not
        bool inversion;
    };

    struct Block {
        Block(MATNew::block b, const std::vector< uint32_t >& blockConsensusSeq);
        Block(size_t blockId, std::string seq);

        int32_t primaryBlockId;
        int32_t secondaryBlockId;

        std::vector< uint32_t > consensusSeq;
        std::string chromosomeName;
    };

    struct GapList {

        std::vector< uint32_t > nucPosition;
        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        std::vector< uint32_t > nucGapLength;

    };

    struct BlockGapList {
        std::vector< uint32_t > blockPosition;
        std::vector< uint32_t > blockGapLength;
    };

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
        bool checkForCyclesHelper(size_t nodeId, std::vector< int >& color);
        void topologicalSortHelper(size_t nodeId, std::vector< size_t >& topoArray, std::vector< bool >& visited);
    public:
        // Graph adjacency list
        size_t numNodes;
        std::vector< std::vector< size_t > > adj;
        std::unordered_map< std::string, std::vector< size_t > > intSequences;

        std::unordered_map< std::string, std::vector< std::string > > paths;
        
        // Added as a patch to incorporate strands
        std::unordered_map< std::string, std::vector< int > > strandPaths;

        std::unordered_map< std::string, std::string > stringIdToConsensusSeq;
        std::unordered_map< std::string, std::vector< std::pair< size_t, size_t > > > stringIdToGaps;
        std::unordered_map< size_t, std::string > intIdToStringId;
        
        // substitutions
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< size_t , std::vector< std::pair< size_t, std::string > > > > > substitutions;
        // insertions
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< size_t , std::vector< std::tuple< size_t, size_t, std::string > > > > > insertions;
        // deletions
        tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< size_t , std::vector< std::pair< size_t, size_t > > > > > deletions;

        Pangraph(Json::Value& pangraphData);
        bool pathExists(size_t nId1, size_t nId2, std::vector< bool >& visited);
        std::vector< size_t > getTopologicalSort();
        std::unordered_map< std::string,std::vector< int > > getAlignedSequences(const std::vector< size_t >& topoArray);
        
        // Patch to incorporate strands
        std::unordered_map< std::string,std::vector< int > > getAlignedStrandSequences(const std::vector< size_t >& topoArray);

        bool checkForCycles();
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
            size_t m_currInternalNode{ 0 };
            size_t m_numLeaves{ 0 };
            size_t m_maxDepth{ 0 };
            float m_meanDepth{ 0 };

            std::unordered_map<std::string, std::vector< std::string > > annotationsToNodes;
            std::unordered_map<std::string, Node*> allNodes;

            Node* createTreeFromNewickString(std::string newick);
            void assignMutationsToNodes(Node* root, size_t& currentIndex, std::vector< MATNew::node >& nodes);
            int getTotalParsimonyParallel(NucMutationType nucMutType, BlockMutationType blockMutType = NONE);
            void invertTree(Node* root);
            void compressTreeParallel(Node* node, size_t level);
            void mergeNodes(Node* par, Node* chi);
            bool debugSimilarity(const std::vector< NucMut > array1, const std::vector< NucMut > array2);
            void dfsExpansion(Node* node, std::vector< Node* >& vec);
            Node* transformHelper(Node* node);
            void adjustLevels(Node* node);
            void setupGlobalCoordinates();

            std::vector< Node* > allLeaves;

            std::string newInternalNodeId() {
                return "node_" + std::to_string(++m_currInternalNode);
            }

        public:
            Tree(const MATNew::tree& mainTree);
            Tree(std::ifstream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
            Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype = FILE_TYPE::GFA);

            void protoMATToTree(const MATNew::tree& mainTree);

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
            void printVCFParallel(std::string reference, std::ofstream& fout);

            Node* subtreeExtractParallel(std::vector< std::string > nodeIds);
            void writeToFile(std::ofstream& fout, Node* node = nullptr);
            std::string getNewickString(Node* node);
            std::string getStringFromReference(std::string reference, bool aligned = true, bool incorporateInversions=true);
            void getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists, std::string reference);
            

            // get unaligned global coordinate
            int32_t getUnalignedGlobalCoordinate(int32_t primaryBlockId, int32_t secondaryBlockId, int32_t pos, int32_t gapPos, const sequence_t& sequence, const blockExists_t& blockExists);
            std::tuple< int, int, int, int > globalCoordinateToBlockCoordinate(int64_t globalCoordinate, const sequence_t& sequence, const blockExists_t& blockExists);

            std::string getSequenceFromVCF(std::string sequenceId, std::ifstream& fin);
            bool verifyVCFFile(std::ifstream& fin);
            void vcfToFASTA(std::ifstream& fin, std::ofstream& fout);
            void annotate(std::ifstream& fin);
            std::vector< std::string > searchByAnnotation(std::string annotation);
            void convertToGFA(std::ofstream& fout);
            void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);
            void getNodesPreorder(PangenomeMAT::Node* root, MATNew::tree& treeToWrite);
            size_t getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition, int nucGapPosition);

            // Transforms tree such that given node becomes child of new root
            void transform(Node* node);
            void reroot(std::string sequenceName);

            Node *root;
            std::vector< Block > blocks;
            std::vector< GapList > gaps;
            BlockGapList blockGaps;

    };

    struct ComplexMutation{
        char mutationType;
        size_t treeIndex1, treeIndex2, treeIndex3;
        std::string sequenceId1, sequenceId2;
        
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

        ComplexMutation(char mutType, int tIndex1, int tIndex2, int tIndex3, std::string sId1, std::string sId2, std::tuple< int,int,int,int > t1, std::tuple< int,int,int,int > t2, std::tuple< int,int,int,int > t3, std::tuple< int,int,int,int > t4){
            mutationType = mutType;
            treeIndex1 = tIndex1;
            treeIndex2 = tIndex2;
            treeIndex3 = tIndex3;

            sequenceId1 = sId1;
            sequenceId2 = sId2;

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

        ComplexMutation(MATNew::complexMutation cm){
            mutationType = (cm.mutationtype()? 'H': 'R');
            treeIndex1 = cm.treeindex1();
            treeIndex2 = cm.treeindex2();
            treeIndex3 = cm.treeindex3();
            sequenceId1 = cm.sequenceid1();
            sequenceId2 = cm.sequenceid2();

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

        MATNew::complexMutation toProtobuf(){
            MATNew::complexMutation cm;
            cm.set_mutationtype(mutationType == 'H');
            cm.set_treeindex1(treeIndex1);
            cm.set_treeindex2(treeIndex2);
            cm.set_treeindex3(treeIndex3);
            cm.set_sequenceid1(sequenceId1);
            cm.set_sequenceid2(sequenceId2);

            if(secondaryBlockIdStart1 != -1){
                cm.set_blockgapexiststart1(true);
                cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32)+secondaryBlockIdStart1);
            } else {
                cm.set_blockgapexiststart1(false);
                cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32));
            }
            cm.set_nucpositionstart1(nucPositionStart1);

            if(nucGapPositionStart1 != -1){
                cm.set_nucgapexiststart1(true);
                cm.set_nucgappositionstart1(nucGapPositionStart1);
            }

            if(secondaryBlockIdStart2 != -1){
                cm.set_blockgapexiststart2(true);
                cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32)+secondaryBlockIdStart2);
            } else {
                cm.set_blockgapexiststart2(false);
                cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32));
            }
            cm.set_nucpositionstart2(nucPositionStart2);

            if(nucGapPositionStart2 != -1){
                cm.set_nucgapexiststart2(true);
                cm.set_nucgappositionstart2(nucGapPositionStart2);
            }

            if(secondaryBlockIdEnd1 != -1){
                cm.set_blockgapexistend1(true);
                cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32)+secondaryBlockIdEnd1);
            } else {
                cm.set_blockgapexistend1(false);
                cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32));
            }
            cm.set_nucpositionend1(nucPositionEnd1);

            if(nucGapPositionEnd1 != -1){
                cm.set_nucgapexistend1(true);
                cm.set_nucgappositionend1(nucGapPositionEnd1);
            }

            if(secondaryBlockIdEnd2 != -1){
                cm.set_blockgapexistend2(true);
                cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32)+secondaryBlockIdEnd2);
            } else {
                cm.set_blockgapexistend2(false);
                cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32));
            }
            cm.set_nucpositionend2(nucPositionEnd2);

            if(nucGapPositionEnd2 != -1){
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
        
        void printFASTA(std::ofstream& fout);
        void writeToFile(std::ofstream& fout);
        void printComplexMutations();
    };

};

#endif // PANGENOME_MAT_NEW_HPP