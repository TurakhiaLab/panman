// #ifndef CORE_LIKELIHOOD_FP64
// #include "core_likelihood_fp64.hpp"
// #endif

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
    std::vector<float> freq;
    freq.resize(STATES);
    float total = 0.0000000001;

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
    float total_pi = 0;
    for (int i = 0; i < STATES; i++)
    {
        freq[i] /= total;
        total_pi += freq[i];
        // std::cout << freq[i] << std::endl;
    }

    pi.resize(STATES);
    pi = freq;

    // Normalize Subs Param
    std::vector<float> subs_param_normalized(subs_param.size());
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
            float alpha = subs_param_normalized[k++];
            rate_matrix[i][j] = rate_matrix[j][i] = alpha * sqrt(pi[i] * pi[j]);
            rate_matrix[i][i] -= alpha * pi[j];
            rate_matrix[j][j] -= alpha * pi[i];
        }
    }

    float mean = 0;
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

void utility::matrix_exp(float bl, std::vector<std::vector<float>>& mat_out)
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
            float value = 0;
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

void scale(size_t i, std::vector<std::vector<float>>& bottom, utility::SCALE_TYPE SCALE_TYPE_ASSIGN, float max_prob, float sum_prob, utility::Node* node)
{
    for (size_t j=0; j<STATES; j++)
    {
        if (SCALE_TYPE_ASSIGN == utility::MAX) bottom[i][j] /= max_prob;
        else if (SCALE_TYPE_ASSIGN == utility::SUM1) bottom[i][j] /= sum_prob;
        else if (SCALE_TYPE_ASSIGN == utility::NONE) continue;
    }
    if (SCALE_TYPE_ASSIGN == utility::MAX) node->scale_vector[i] *= max_prob;
    else if (SCALE_TYPE_ASSIGN == utility::SUM1) node->scale_vector[i] *= sum_prob;
   

}


void postorder_traversal(utility::Node* node, utility::Tree& tree)
{
    std::vector<std::vector<float>> bottom;
    std::string identifier = node->identifier;
    bottom.resize(utility::length);

    node->scale_vector.resize(utility::length, 1.0);

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
        float bl = child->branchLength;
        std::vector<std::vector<float>> mat_out;
        utility::matrix_exp(bl, mat_out);

        // for (int i = 0; i < STATES; i++)
        // {
        //     for (int j = 0; j < STATES; j++)
        //     {
        //         if (mat_out[i][j] < 0)
        //         {
        //             std::cout << "Issue: " << mat_out[i][j] << std::endl;
        //         }
        //     }
        // }

        std::vector<std::vector<float>> child_bottom = child->bottom;

        for (size_t i = 0; i < bottom.size(); i++)
        {
            float max_prob = 0.;
            float sum_prob = 0.;
            for (size_t j = 0; j < STATES; j++)
            {
                float prob = 0;
                for (size_t k = 0; k < STATES; k++)
                {
                    // std::cout << i << " " << j << " " << k << std::endl;
                    prob += child_bottom[i][k]*mat_out[j][k];
                }
                bottom[i][j] *= prob;
                sum_prob += bottom[i][j];
                if (max_prob<bottom[i][j]) max_prob = bottom[i][j];
            }
            scale(i, bottom, utility::SCALE_TYPE_ASSIGN, max_prob, sum_prob, node);
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

double postorder_traversal_scaling(utility::Node* node)
{
    double lk = 0.0;
    for (auto &child: node->children) lk += postorder_traversal_scaling(child);

    for (size_t i=0; i<node->scale_vector.size(); i++)
    {
        lk += log(node->scale_vector[i]);
    }
    return lk;
}

void utility::felsenstein_pruning(utility::Tree& tree)
{

    utility::bottom_up(tree);

    // Multiply conditional prob with stationary prob at root 
    utility::Node* root = tree.root;
    std::vector<std::vector<float>> bottom = root->bottom;
    double lk = 0;

    for (size_t i = 0; i < bottom.size(); i++)
    {
        float lk_node = 0;
        for (size_t j = 0; j < STATES; j++)
        {
            // std::cout << bottom[i][j] << "\t" << pi[j] << "\t";
            bottom[i][j] *= pi[j];
            lk_node += bottom[i][j];
        }
        // std::cout << "\n";

        // if (lk_node < 0)
        // {
        //     std::cout << log(lk_node) << " " << lk_node << std::endl;
        // }
        lk += (double)log(lk_node);
    }

    // Adding scaling factor
    lk += postorder_traversal_scaling(root);


    printf("%lf\n", lk);

}

void top_down_preorder_traversal(utility::Node* node, utility::Tree& tree, std::vector<utility::Node*> siblings)
{
    size_t i,j,k,m,l;
    std::vector<std::vector<float>> up;
    up.resize(utility::length);
    // std::cout << node->identifier << "\t";
    std::vector<std::vector<float>> bottom = node->bottom;

    for (i=0; i<siblings.size(); i++)
    {
        if (siblings[i] == node)
        {
            siblings[i] = siblings.back();
            siblings.pop_back();
            break;
        }
    }

    if (PRINT_CHILDREN)
    {
        std::cout << node->identifier << "\t";
        for (i=0; i<siblings.size(); i++)
        {
            std::cout << siblings[i]->identifier << " ";
        }
        std::cout << "\n";
    }

    for (i=0; i<utility::length; i++) up[i].resize(STATES,1); // Initializing with 1;

    if (node != tree.root)
    {
        std::vector<std::vector<float>> mat_out_curr;
        utility::matrix_exp(node->branchLength, mat_out_curr);

        std::vector<std::vector<std::vector<float>>> mat_out;
        std::vector<std::vector<std::vector<float>>> sib_bottom;
        mat_out.resize(siblings.size());
        sib_bottom.resize(siblings.size());
        for (i=0; i<siblings.size(); i++)
        {
            utility::matrix_exp(siblings[i]->branchLength, mat_out[i]);
            sib_bottom[i] = siblings[i]->bottom;
        }

        std::vector<std::vector<float>> parent_up = node->parent->up;
        for (m=0; m<utility::length; m++)
        {
            for (i=0; i<STATES; i++)
            {
                float prob = 0;
                for (j=0; j<STATES; j++)
                {
                    float prob_sib = 1;
                    int sib_count = 0;
                    for(l=0; l<siblings.size(); l++)
                    {
                        float prob_sib_local = 0;
                        for (k=0; k<STATES; k++)
                        {
                            prob_sib_local += sib_bottom[sib_count][m][k]*mat_out[sib_count][j][k];
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
    std::vector<std::vector<float>> marginal;
    marginal.resize(utility::length);
    std::vector<std::vector<float>> up = node->up;
    std::vector<std::vector<float>> bottom = node->bottom;

    for (i=0; i<utility::length; i++) marginal[i].resize(STATES,1); // Initializing with 1;
    node->inference.resize(utility::length);
    for (i=0; i<utility::length; i++)
    {
        float total_prob = 0;
        for (j=0; j<STATES; j++) total_prob += utility::pi[j]*up[i][j]*bottom[i][j];
        for (j=0; j<STATES; j++) marginal[i][j] = (utility::pi[j]*up[i][j]*bottom[i][j])/total_prob;
        // infer characters
        float max = 0;
        for (j=0; j<STATES; j++) { if (marginal[i][j] > max) max = marginal[i][j]; }
        for (j=0; j<STATES; j++) { if (marginal[i][j] == max) node->inference[i].push_back(int8_t(pow(2,j)));}

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

char int2char(int v)
{
    switch (v)
    {
    case 1:
        return 'A';
    case 2:
        return 'C';
    case 4:
        return 'G';
    case 8:
        return 'T';
    case 16:
        return '-';
    default:
        return '-';
    }

}

void utility::printLikelihoodInference(utility::Tree& tree, utility::Node* node)
{
    
    size_t index = 0;
    std::cout << node->identifier << "\t";
    for (index=0; index<utility::length; index++)
    {
        char tip_char = utility::seqs.find(node->identifier)->second.second[index];
        std::cout << tip_char;
    }
    std::cout << "\n";


    while (node->parent != nullptr)
    {   
        node = node->parent;
        std::cout << node->identifier << "\t";
        for (index=0; index<utility::length; index++)
        {
            char tip_char = int2char(node->inference[index][0]);
            std::cout << tip_char;
        }
        std::cout << "\n";
    }
}

int8_t char2int(char& c)
{
    switch (c)
    {
    case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 4;
    case 'T':
        return 8;
    case '-':
        return 16;
    default:
        return 16;
    }
}

int8_t fitch_op(int8_t& a, int8_t&b)
{
    int8_t c = ((a&b) == 0) ? (a|b) : (a&b);
    return c;
}

void bottom_up_fitch(utility::Tree& tree, utility::Node* node)
{
    size_t i;
    std::vector<int8_t> inference(utility::length, 0);
    if (!node->children.size()) //tips
    {
        std::string seq = utility::seqs[node->identifier].second;
        for (i=0; i<utility::length; i++)
        {
            inference[i]= char2int(seq[i]);
        }
    }
    else //internal nodes
    {
        bool child_count = false;
        for (auto &child: node->children)
        {
            bottom_up_fitch(tree, child);
            std::vector<int8_t> child_inference = child->fitch_inference;
            if (!child_count)
            {
                child_count = true;
                for (i=0; i<utility::length; i++)
                {
                    inference[i] = child_inference[i];
                }
                continue;
            }
            for (i=0; i<utility::length; i++)
            {
                inference[i] = fitch_op(inference[i], child_inference[i]);
            }
        }
    }
    node->fitch_inference = inference;
}

int8_t accInt82Int8(int8_t& acc)
{
    if (acc >= 16) return 16;
    else if (acc >= 8) return 8;
    else if (acc >= 4) return 4;
    else if (acc >= 2) return 2;
    else return 1;
}

void top_dowm_fitch(utility::Tree& tree, utility::Node* node)
{
    size_t i;
    std::vector<int8_t> inference = node->fitch_inference;

    if (node == tree.root)
    {
        for (i=0; i<utility::length; i++)
        {
            inference[i] = accInt82Int8(inference[i]);
        }
    }
    else
    {
        std::vector<int8_t> parent_inference = node->parent->fitch_inference;
        std::vector<int8_t> inference = node->fitch_inference;
        for (i=0; i<utility::length; i++)
        {
            inference[i] = ((parent_inference[i]&inference[i]) == 0) ? accInt82Int8(inference[i]) : parent_inference[i];
        }
    }
    node->fitch_inference = inference;

    for (auto &child: node->children)
    {
        top_dowm_fitch(tree, child);
    }
    return;
}




void utility::fitch(utility::Tree& tree)
{
    utility::Node* root = tree.root;
    bottom_up_fitch(tree, root);
    top_dowm_fitch(tree, root);
}


void utility::printParsimonyInference(utility::Tree& tree, utility::Node* node)
{
    size_t index = 0;
    std::cout << node->identifier << "\t";
    for (index=0; index<utility::length; index++)
    {
        char tip_char = utility::seqs.find(node->identifier)->second.second[index];
        std::cout << tip_char;
    }
    std::cout << "\n";


    while (node->parent != nullptr)
    {   
        node = node->parent;
        std::cout << node->identifier << "\t";
        for (index=0; index<utility::length; index++)
        {
            char tip_char = int2char(node->fitch_inference[index]);
            std::cout << tip_char;
        }
        std::cout << "\n";
    }
}