// #include "core_likelihood.hpp"

void utility::msa_seq(std::string input_file)
{
    std::ifstream in(input_file);

    if(!in.good()){
        std::cerr << "Error opening " << input_file << "." << std::endl;
    }

    std::string line, name, content;
    length = 0;
    int count = 0;
    while( std::getline( in, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                seqs[name] = std::make_pair(count,content);
                count++;
                if (length != 0 && length != content.size())
                {
                    std::cerr << "Length of aligned sequences is unequal." << std::endl;
                }
                else 
                {
                    length = content.size();
                }
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }

    if( !name.empty() ){ // Print out what we read from the last entry
        seqs[name] = std::make_pair(count,content);
        count++;
    }

}

void utility::rate_matrix_calc()
{
    // Calculate Rate Matrix
    tbb::concurrent_vector<tbb::concurrent_vector<int>> freq_vec( STATES, tbb::concurrent_vector<int> (seqs.size(), 0));
    std::vector<double> freq;
    freq.resize(STATES);
    double total;

    tbb::parallel_for_each(seqs, [&](auto& seq){
        std::string seq_ = seq.second.second;
        int index = seq.second.first;

        for (size_t i = 0; i < length; i++)
        {
            if (seq_[i] == 'A' || seq_[i] == 'a' )
            {
                freq_vec[0][index]++;
            }
            else if (seq_[i] == 'C' || seq_[i] == 'c' )
            {
                freq_vec[1][index]++;
            }
            else if (seq_[i] == 'G' || seq_[i] == 'g' )
            {
                freq_vec[2][index]++;
            }
            else if (seq_[i] == 'T' || seq_[i] == 't' )
            {
                freq_vec[3][index]++;
            }
            else
            {
                if (STATES == 5)
                {
                    freq_vec[4][index]++;
                }
            }
        }
    });

    for (int i = 0; i < STATES; i++)
    {
        freq[i] = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, freq_vec[0].size()), 
            0, [&](tbb::blocked_range<int>& r, int init) -> int
            {
                for(int j = r.begin(); j != r.end(); j++)
                {
                    init += freq_vec[i][j];
                }
            return init;
            },
            std::plus<int>()
        );
        total += freq[i];
    }
    
    // Handle corner cases: 1) What if freq[i] = 0?
    double total_pi = 0;
    for (int i = 0; i < STATES; i++)
    {
        freq[i] /= total;
        total_pi += freq[i];
        // std::cout << freq[i] << std::endl;
    }

    pi.resize(STATES);
    pi = freq;

    // Normalize Subs Param
    std::vector<double> subs_param_normalized(subs_param.size());
    for (size_t i = 0; i < subs_param.size(); i++)
    {
        subs_param_normalized[i] = subs_param[i]/subs_param[subs_param.size() - 1];
    }

    // for (size_t i = 0; i < subs_param.size(); i++)
    // {
    //     std::cout << "Factor: " << subs_param_normalized[i] << "\n";
    // } 

    rate_matrix.resize(STATES);
    for (int i = 0; i < STATES; i++)
    {   
        rate_matrix[i].resize(STATES);
        rate_matrix[i][i] = 0;
    }

    int k = 0;
    for (int i = 0; i < STATES; i++)
    {
        for (int j = i+1; j < STATES; j++)
        {
            double alpha = subs_param_normalized[k++];
            rate_matrix[i][j] = rate_matrix[j][i] = alpha * sqrt(pi[i] * pi[j]);
            rate_matrix[i][i] -= alpha * pi[j];
            rate_matrix[j][j] -= alpha * pi[i];
        }
    }

    double mean = 0;
    for (int i = 0; i < STATES; ++i) mean += pi[i] * (-rate_matrix[i][i]);
    for (int i = 0; i < STATES; ++i)
        for (int j = 0; j < STATES; ++j)
            rate_matrix[i][j] /= mean;
    
    

    // double norm = 0;
    // for (int i = 0; i < STATES; i++)
    // {
    //     rate_matrix[i][i] = rate_matrix[i][i] - total_pi;  
    //     norm += std::abs(rate_matrix[i][i])*pi[i];
    // }


}

void utility::matrix_exp(double bl, std::vector<std::vector<double>>& mat_out)
{
    // Approximating e^(l*x) = I + l*x + (l*x)^2/2;
    mat_out.resize(STATES);
    for (size_t i = 0; i < rate_matrix.size(); i++)
    {
        mat_out[i].resize(STATES);
        for (size_t j = 0; j < rate_matrix[i].size(); j++)
        {
            mat_out[i][j] = bl*rate_matrix[i][j];
            if (i == j)
            {
                mat_out[i][j] += 1;
            }
        }
    }

    for (size_t i = 0; i < rate_matrix.size(); i++)
    { 
        for (size_t j = 0; j < rate_matrix[i].size(); j++)
        { 
            double value = 0;
            for (size_t k = 0; k < rate_matrix[i].size(); k++)
            { 
                value += rate_matrix[i][k] * rate_matrix[k][j]; 
            } 
            mat_out[i][j] += value*bl*bl*(1/2);
        }   
    } 
    // printf ("P-matrix for branch length %f\n", bl);

    // for (size_t i = 0; i < mat_out.size(); i++)
    // {
    //     for (size_t j = 0; j < mat_out[i].size(); j++)
    //     {
    //         std::cout << mat_out[i][j] << "\t";
    //     }
    //     printf ("\n");
    // }

}


void postorder_traversal(utility::Node* node, utility::Tree& tree)
{
    std::vector<std::vector<double>> bottom;
    std::string identifier = node->identifier;
    bottom.resize(utility::length);

    // std::cout << node->identifier << "\t";
    // for (auto &child:node->children)
    // {
    //     std::cout << child->identifier << "\t";
    // }
    // std::cout << "\n";
    
    if (node->children.size() == 0)
    {
        std::string seq = utility::seqs[identifier].second;
        for (size_t i = 0; i < bottom.size(); i++)
        {
            bottom[i].resize(STATES);
            for (size_t j= 0; j < STATES; j++)
            {
                bottom[i][j] = 0;
            }
        }    
        for (size_t i = 0; i < bottom.size(); i++)
        {
            if (seq[i] == 'A')
            {
                bottom[i][0] = 1;
            }   
            else if (seq[i] == 'C') 
            {
                bottom[i][1] = 1;
            }
            else if (seq[i] == 'G') 
            {
                bottom[i][2] = 1;
            }
            else if (seq[i] == 'T') 
            {
                bottom[i][3] = 1;
            }
            else // N are considered as -
            {
                if (STATES == 5)
                {
                    bottom[i][4] = 1;
                }
            }
        }
        node->bottom = bottom;
        return;
    }
    
    for (size_t i = 0; i < bottom.size(); i++)
    {
        bottom[i].resize(STATES);
        for (int j= 0; j < STATES; j++)
        {
            // Initialize to 1 as child prob are multiplied
            bottom[i][j] = 1;
        }
    }

    for (auto child: node->children)
    {
        postorder_traversal(child, tree);
        double bl = child->branchLength;
        std::vector<std::vector<double>> mat_out;
        utility::matrix_exp(bl, mat_out);

        for (int i = 0; i < STATES; i++)
        {
            for (int j = 0; j < STATES; j++)
            {
                if (mat_out[i][j] < 0)
                {
                    std::cout << "Issue: " << mat_out[i][j] << std::endl;
                }
            }
        }

        std::vector<std::vector<double>> child_bottom = child->bottom;

        for (size_t i = 0; i < bottom.size(); i++)
        {
            for (size_t j = 0; j < STATES; j++)
            {
                double prob = 0;
                for (size_t k = 0; k < STATES; k++)
                {
                    // std::cout << i << " " << j << " " << k << std::endl;
                    prob += child_bottom[i][k]*mat_out[j][k];
                }
                bottom[i][j] *= prob;
            }
        }
    }
    
    node->bottom = bottom;
    
}

void utility::bottom_up(utility::Tree& tree)
{
    utility::Node* root = tree.root;
    postorder_traversal(root, tree);
    return;
}

void utility::felsenstein_pruning(utility::Tree& tree)
{
    utility::bottom_up(tree);

    // Multiply conditional prob with stationary prob at root 
    utility::Node* root = tree.root;
    std::vector<std::vector<double>> bottom = root->bottom;
    double lk = 0;

    for (size_t i = 0; i < bottom.size(); i++)
    {
        double lk_node = 0;
        for (size_t j = 0; j < STATES; j++)
        {
            bottom[i][j] *= pi[j];
            lk_node += bottom[i][j];
        }
        // if (lk_node < 0)
        // {
            // std::cout << log(lk_node) << " " << lk_node << std::endl;
        // }
        lk += log(lk_node);
    }

    std::cout << "Likelihood Value:" << lk << std::endl;

}

void top_down_preorder_traversal(utility::Node* node, utility::Tree& tree, std::vector<utility::Node*> siblings)
{
    size_t i,j,k,m;
    std::vector<std::vector<double>> up;
    up.resize(utility::length);

    std::vector<std::vector<double>> bottom = node->bottom;

    for (i=0; i<siblings.size(); i++)
    {
        if (siblings[i] == node)
        {
            siblings[i] = siblings.back();
            siblings.pop_back();
            break;
        }
    }

    std::cout << node->identifier << "\t";
    for (i=0; i<siblings.size(); i++)
    {
        std::cout << siblings[i]->identifier << " ";
    }
    std::cout << "\n";

    for (i=0; i<utility::length; i++) up[i].resize(STATES,1); // Initializing with 1;

    if (node != tree.root)
    {
        std::vector<std::vector<double>> mat_out_curr;
        utility::matrix_exp(node->branchLength, mat_out_curr);

        std::vector<std::vector<std::vector<double>>> mat_out;
        mat_out.resize(siblings.size());
        for (i=0; i<siblings.size(); i++)
        {
            utility::matrix_exp(siblings[i]->branchLength, mat_out[i]);
        }

        std::vector<std::vector<double>> parent_up = node->parent->up;
        for (m=0; m<utility::length; m++)
        {
            for (i=0; i<STATES; i++)
            {
                double prob = 0;
                for (j=0; j<STATES; j++)
                {
                    double prob_sib = 0;
                    int sib_count = 0;
                    for(auto &sib: siblings)
                    {
                        double prob_sib_local = 0;
                        std::vector<std::vector<double>> sib_bottom = sib->bottom;
                        for (k=0; k<STATES; k++)
                        {
                            prob_sib_local += sib_bottom[m][k]*mat_out[sib_count][j][k];
                        }
                        prob_sib *= prob_sib_local;
                        sib_count++;
                    }
                    prob += parent_up[m][j]*mat_out_curr[i][j]*prob_sib;
                }
                up[m][i] *= prob;
            }
        }
    }

    node->up = up;
    for (auto &child: node->children)
    {
        std::vector<utility::Node*> node_siblings = node->children;
        top_down_preorder_traversal(child, tree, node_siblings);
    }

    return;
}


void utility::top_down(utility::Tree& tree)
{

    utility::Node* root = tree.root;
    std::vector<utility::Node*> siblings;
    top_down_preorder_traversal(root, tree, siblings);

}

void marginal_preorder_traversal(utility::Node* node, utility::Tree& tree)
{
    size_t i,j;
    std::vector<std::vector<double>> marginal;
    marginal.resize(utility::length);
    std::vector<std::vector<double>> up = node->up;
    std::vector<std::vector<double>> bottom = node->bottom;

    for (i=0; i<utility::length; i++) marginal[i].resize(STATES,1); // Initializing with 1;

    for (i=0; i<utility::length; i++)
    {
        double total_prob = 0;
        for (j=0; j<STATES; j++) total_prob += utility::pi[j]*up[i][j]*bottom[i][j];
        for (j=0; j<STATES; j++) marginal[i][j] = (utility::pi[j]*up[i][j]*bottom[i][j])/total_prob;

        // infer characters
        double max = 0;
        for (j=0; j<STATES; j++)
        {
            if (marginal[i][j] > max) max = marginal[i][j];
        }
        for (j=0; j<STATES; j++)
        {
            if (marginal[i][j] == max) node->inference.push_back(j);
        }
    }
    node->marginal = marginal;

    for (auto &child: node->children)
    {
        marginal_preorder_traversal(child, tree);
    }

    return;
}


void utility::marginal(utility::Tree& tree)
{

    utility::Node* root = tree.root;
    marginal_preorder_traversal(root, tree);

}


