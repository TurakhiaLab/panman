#include "core_likelihood_fp32.hpp"
#include "core_likelihood_fp32.cpp"
#include "tree_fp32.cpp"

int main(int argc, char** argv)
{
    if (argc != 3)
        std::cerr << "Usage: "<< argv[0] <<" [alignment file] [newick tree file]" << std::endl;

    utility::msa_seq(argv[1]);    // Store MSA data into data-structure
    utility::subs_param.resize(10); // Set Subs parameter
    for (size_t i = 0; i < 10; i++) utility::subs_param[i] = 1;
    utility::rate_matrix_calc(); // Find the rate matrix

    
    if (PRINT_RATE_MATRIX) // Print matrix_exp Matrix
    {
        for (size_t i = 0; i < utility::rate_matrix.size(); i++)
        {
            for (size_t j = 0; j < utility::rate_matrix[0].size(); j++) std::cout << utility::rate_matrix[i][j] << "\t";
            std::cout << "\n";
        }
        for (auto &p: utility::pi) std::cout << p << "\n";
    }


    std::ifstream newickTree(argv[2]);
    std::string newickString;
    newickTree >> newickString;
    utility::Tree tree(newickString);

    utility::felsenstein_pruning(tree);

    // utility::bottom_up(tree);
    // std::cout << "\n";
    // utility::top_down(tree);
    // std::cout << "\n";
    // utility::marginal(tree);


    // utility::fitch(tree);

    // if (PRINT_INFERENCE)
    // {
    //     std::string tip="Mon";
    //     auto search = tree.allNodes.find(tip);
    //     if (search != tree.allNodes.end()) utility::printLikelihoodInference(tree, search->second);
    //     else {std::cerr<<"Tip not found";}
    // }

    return 0;
}

