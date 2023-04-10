#ifndef PANGENOME_MAT_NEW_HPP
#define PANGENOME_MAT_NEW_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <json/json.h>
#include "mutation_annotation_test_proto3_optional_new.pb.h"
#include "vg.pb.h"
#include "AuxilaryMAT.hpp"
#include "spoa/spoa.hpp"

#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"

namespace PangenomeMAT2 {

    char getNucleotideFromCode(int code);
    char getCodeFromNucleotide(char nuc);
    void printSequenceLines(const std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
        const std::vector< std::pair< bool, std::vector< bool > > >& blockExists, size_t lineSize, bool aligned, std::ofstream& fout);
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
        FASTA = 5
    };

    enum NucMutationType {
        NS = 0,
        ND = 1,
        NI = 2,
        NSNPS = 3,
        NSNPI = 4,
        NSNPD = 5,
    };

    enum BlockMutationType {
        BI = 1,
        BD = 0,
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
                case PangenomeMAT2::NucMutationType::NSNPS:
                    mutInfo += PangenomeMAT2::NucMutationType::NS;
                    break;
                case PangenomeMAT2::NucMutationType::NSNPI:
                    mutInfo += PangenomeMAT2::NucMutationType::NI;
                    break;
                case PangenomeMAT2::NucMutationType::NSNPD:
                    mutInfo += PangenomeMAT2::NucMutationType::ND;
                    break;
                case PangenomeMAT2::NucMutationType::NS:
                    mutInfo += PangenomeMAT2::NucMutationType::NS;
                    break;
                case PangenomeMAT2::NucMutationType::NI:
                    mutInfo += PangenomeMAT2::NucMutationType::NI;
                    break;
                case PangenomeMAT2::NucMutationType::ND:
                    mutInfo += PangenomeMAT2::NucMutationType::ND;
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
        }

        BlockMut(size_t blockId, bool type){
            primaryBlockId = blockId;
            secondaryBlockId = -1;
            blockMutInfo = type;
        }
        BlockMut(){
            
        }

        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        bool  blockMutInfo;
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
            void getNodesPreorder(PangenomeMAT2::Node* root, MATNew::tree& treeToWrite);
            AuxilaryMAT::Node* convertToAuxMatHelper(PangenomeMAT2::Node* currentNode, std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
                std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > >& coordinates,\
                std::vector< std::pair< bool, std::vector< bool > > >& blockExists
            );

            std::vector< Node* > allLeaves;

            std::string newInternalNodeId() {
                return "node_" + std::to_string(++m_currInternalNode);
            }

        public:
            Tree(std::ifstream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
            Tree(std::ifstream& fin, std::ifstream& secondFin, FILE_TYPE ftype = FILE_TYPE::GFA);

            int nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states);
            void nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState);
            void nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, std::pair< PangenomeMAT2::NucMutationType, char > >& mutations, int parentState);

            int blockFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states);
            void blockFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states, int parentState);
            void blockFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states, std::unordered_map< std::string, bool >& mutations, int parentState);

            void printSummary();
            void printBfs(Node* node = nullptr);
            void printFASTA(std::ofstream& fout, bool aligned = false);
            void printFASTAParallel(std::ofstream& fout, bool aligned = false);
            Node* subtreeExtractParallel(std::vector< std::string > nodeIds);
            void writeToFile(std::ofstream& fout, Node* node = nullptr);
            std::string getNewickString(Node* node);
            std::string getStringFromReference(std::string reference, bool aligned = true);
            void printVCFParallel(std::string reference, std::ofstream& fout);
            std::string getSequenceFromVCF(std::string sequenceId, std::ifstream& fin);
            bool verifyVCFFile(std::ifstream& fin);
            void vcfToFASTA(std::ifstream& fin, std::ofstream& fout);
            void annotate(std::ifstream& fin);
            std::vector< std::string > searchByAnnotation(std::string annotation);
            void convertToVG(std::ofstream& fout);
            void convertToGFA(std::ofstream& fout);
            void printFASTAFromVG(std::ifstream& fin, std::ofstream& fout);
            void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);

            AuxilaryMAT::Tree* convertToAuxMat();


            Node *root;
            std::vector< Block > blocks;
            std::vector< GapList > gaps;
            BlockGapList blockGaps;

    };

};

#endif // PANGENOME_MAT_NEW_HPP