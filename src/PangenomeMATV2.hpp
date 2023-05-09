#ifndef PANGENOME_MAT_NEW_HPP
#define PANGENOME_MAT_NEW_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <iterator>
#include <algorithm>
#include <atomic>
#include "mutation_annotation_test_proto3_optional_new.pb.h"
#include "AuxilaryMAT.hpp"
#include "syncmer.hpp"
#include <math.h>

#define PMAT_VERSION "1.0-beta"
#define VCF_VERSION "4.2"

namespace PangenomeMAT2 {

    char getNucleotideFromCode(int code);
    void printSequenceLines(const std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > >& sequence,\
        const std::vector< std::pair< bool, std::vector< bool > > >& blockExists, size_t lineSize, bool aligned, std::ofstream& fout);
    std::pair< int, int > replaceMutation(std::pair<int,int> oldMutation, std::pair<int, int> newMutation);
    std::string stripGaps(std::string sequenceString);
    std::string getDate();
    std::string stripString(std::string s);

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
            }

            nucPosition = std::get<2>(mutationArray[start]);
            nucGapPosition = std::get<3>(mutationArray[start]);

            nucs = 0;
            for(int i = start; i < end; i++){
                nucs += (std::get<5>(mutationArray[i]) << (4*(5-(i - start))));
            }
        }

        NucMut(MATNew::nucMut mutation){
            nucPosition = mutation.nucposition();
            primaryBlockId = (mutation.blockid() >> 32);
            mutInfo = (mutation.mutinfo() & 0xFF);
            nucs = (mutation.mutinfo() >> 8);
            nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

            if(mutation.blockgapexist()){
                secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
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

        void loadFromProtobuf(MATNew::blockMut mutation){
            primaryBlockId = (mutation.blockid() >> 32);
            if(mutation.blockgapexist()){
                secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
            } else {
                secondaryBlockId = -1;
            }
            blockMutInfo = mutation.blockmutinfo();
        }

        int32_t primaryBlockId;
        int32_t secondaryBlockId;
        bool  blockMutInfo;
    };

    struct Block {
        Block(MATNew::block b);

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
            //void processNode(std::string referenceSequence, std::string altSequence, std::string nid, std::string parentId, std::ofstream& fout);


        public:
            Tree(std::ifstream& fin);
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
            AuxilaryMAT::Tree* convertToAuxMat();
            void indexSyncmers(std::ofstream& fout);
            void placeSample(std::string fastqPath, seedIndex &index);
            void placeDFS(Node *currNode, std::set<kmer_t> currNodeSyncmers, std::set<kmer_t> &querySyncmers, seedIndex &index, dynamicJaccard dj, std::unordered_map<std::string, float> &scores);
            void loadIndex(std::ifstream &indexFile, seedIndex &index);

            Node *root;
            std::vector< Block > blocks;
            std::vector< GapList > gaps;
            BlockGapList blockGaps;

    };

};

#endif // PANGENOME_MAT_NEW_HPP