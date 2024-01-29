#ifndef CORE_LIKELIHOOD_FP64
#define CORE_LIKELIHOOD_FP64

#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <queue>
#include <unordered_map>
#include <string>
#include <utility>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <bits/stdc++.h>


#define NORM 1
#define STATES 5
#define DEBUG 0
#define PRINT_INFERENCE 0
#define PRINT_RATE_MATRIX 0
#define PRINT_CHILDREN 0
#define SCALE_THRESHOLD 0.0000000000000000000542101086242752217003726400434970855712890625

namespace utility
{

    size_t length;
    std::vector<double> pi;
    std::vector<double> subs_param;
    std::vector<std::vector<double>> rate_matrix;
    tbb::concurrent_unordered_map< std::string, std::pair<int, std::string> > seqs;
    void msa_seq(std::string input_file);
    void rate_matrix_calc();
    void matrix_exp(double bl, std::vector<std::vector<double>>& mat_out);
    
    enum SCALE_TYPE { MAX, SUM1, NONE};
    static SCALE_TYPE SCALE_TYPE_ASSIGN = NONE;

    class Node 
    {
    public:
        Node(std::string id, double len);
        Node(std::string id, Node* par, double len);

        double branchLength;
        size_t level;

        std::string identifier;
        Node* parent;
        std::vector< Node* > children;

        // For Parsimony
        std::vector<int8_t> fitch_inference;

        // For Likelihood
        std::vector<std::vector<double>> bottom;
        std::vector<std::vector<double>> up;
        std::vector<std::vector<double>> marginal;
        std::vector<std::vector<int8_t>> inference;

        // Scaling Factor
        std::vector<double> scale_vector;
    };
    
    class Tree
    {
    private:
        size_t m_currInternalNode{ 0 };
        size_t m_numLeaves{ 0 };
        size_t m_maxDepth{ 0 };
        double m_meanDepth{ 0 };
        std::string newInternalNodeId() {
            return "node_" + std::to_string(++m_currInternalNode);
        }

    public:
        Node* root;
        std::unordered_map< std::string, Node* > allNodes;
        Tree(std::string newick);
    };

    void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);
    std::string stripString(std::string s);
    void bottom_up(utility::Tree& tree);
    void top_down(utility::Tree& tree);
    void marginal(utility::Tree& tree);
    void felsenstein_pruning(utility::Tree& tree);
    void fitch(Tree& tree);
    void printLikelihoodInference(Tree& tree, Node* node);
    void printParsimonyInference(Tree& tree, Node* node);
};

#endif