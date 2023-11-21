#include "core_likelihood.hpp"
#include "core_likelihood.cpp"
#include "tree.cpp"

int main(int argc, char** argv)
{
    if (argc != 3)
        std::cerr << "Usage: "<< argv[0] <<" [alignment file] [newick tree file]" << std::endl;

    
    // Store MSA data into data-structure
    utility::msa_seq(argv[1]);

    // Set Subs parameter
    utility::subs_param.resize(10);
    for (size_t i = 0; i < 10; i++)
    {
        utility::subs_param[i] = 1;
    }
    
    // Find the rate matrix
    utility::rate_matrix_calc();

    std::vector<std::vector<double>> mat_out;

    // Print matrix_exp Matrix
    if (DEBUG)
    {
        for (size_t i = 0; i < utility::rate_matrix.size(); i++)
        {
            for (size_t j = 0; j < utility::rate_matrix[0].size(); j++)
            {
                std::cout << utility::rate_matrix[i][j] << "\t";
            }
            std::cout << "\n";
        }

        // Print Pi
        for (auto &p: utility::pi)
        {
            std::cout << p << "\n";
        }
    }


    std::ifstream newickTree(argv[2]);
    std::string newickString;
    newickTree >> newickString;
    utility::Tree tree(newickString);

    // for (auto s: tree.allNodes)
    // {
    //     std::cout << s.first << std::endl;
    // }

    utility::bottom_up(tree);
    utility::top_down(tree);
    utility::marginal(tree);

    return 0;
}

