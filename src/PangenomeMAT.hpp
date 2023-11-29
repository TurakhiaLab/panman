#ifndef PANGENOME_MAT_HPP
#define PANGENOME_MAT_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include "mutation_annotation_test_proto3_optional.pb.h"

#include <json/json.h>
#include "mutation_annotation_test_proto3_optional_new.pb.h"
#include "kseq.hpp"
#include "syncmer.hpp"

#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"

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
        // For creating SNP mutation
        NucMut( const std::tuple< int, int, int, int, int >& mutationInfo ){
            // bid, pos, gapPos, type, char
            position = std::get<1>(mutationInfo);
            gapPosition = std::get<2>(mutationInfo);
            condensed = (std::get<0>(mutationInfo) << 8) + (1<<7) + (std::get<4>(mutationInfo) << 3) + (std::get<3>(mutationInfo));
        }

        // For creating non-SNP mutations from SNP mutations
        NucMut(const std::vector< std::tuple< int, int, int, int, int > >& mutationArray, int start, int end){
            condensed = ((std::get<0>(mutationArray[start])) << 8) + ((end - start) << 3);
            // type
            switch(std::get<3>(mutationArray[start])){
                case PangenomeMAT::NucMutationType::NSNPS:
                    condensed += PangenomeMAT::NucMutationType::NS;
                    break;
                case PangenomeMAT::NucMutationType::NSNPI:
                    condensed += PangenomeMAT::NucMutationType::NI;
                    break;
                case PangenomeMAT::NucMutationType::NSNPD:
                    condensed += PangenomeMAT::NucMutationType::ND;
                    break;
            }

            position = std::get<1>(mutationArray[start]);
            gapPosition = std::get<2>(mutationArray[start]);

            nucs = 0;
            for(int i = start; i < end; i++){
                nucs += ((uint64_t)std::get<4>(mutationArray[i]) << (4*(15-(i - start))));
            }
        }

        NucMut(MAT::nuc_mut mutation){
            position = mutation.position();

            if(mutation.gap_exist()){
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

            std::unordered_map<std::string, std::vector< std::string > > annotationsToNodes;
        public:
            Tree(const PanMAT::tree& mainTree);
            Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
            Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype = FILE_TYPE::GFA);

            // Copy blocks from current tree into new tree which is rooted at one of the internal nodes of the current tree. Used in split for PanMAN
            Tree(Node* newRoot, const std::vector< Block >& b, const std::vector< GapList >& g, const std::unordered_map< std::string, int >& c, const BlockGapList& bgl);

            void protoMATToTree(const PanMAT::tree& mainTree);

            void condenseTree(Node *node);

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
    void placeSample(Tree *T, std::string fastqPath, seedIndex &index, size_t k, size_t s);
    std::set<kmer_t> syncmersFromFastq(std::string fastqPath,  std::vector<read_t> &reads, size_t k, size_t s);
    std::string getCurrentFastaSequence(const sequence_t& sequence, const blockExists_t& blockExists, blockStrand_t& blockStrand, bool aligned);

};

#endif // PANGENOME_MAP_HPP