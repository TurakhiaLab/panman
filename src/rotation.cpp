#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <limits>
#include <chrono>

using namespace std;


std::pair<int, int> rotate_alignment(std::vector<std::string> consensus, std::vector<std::string> sample)
{
    std::vector<std::pair<int, int>> score(sample.size(), std::make_pair(-1,-1));
    std::vector<std::pair<int, int>> next_score(sample.size(), std::make_pair(-1,-1));
    std::pair<int,int>max_(0,0);
    int match = 5;
    int gap = 1;
    int mismatch = 2;
    for(auto i = 0; i < consensus.size(); i++)
    {
        for (auto j = 0; j < sample.size(); j++)
        {
            int up_idx = (j == 0) ? sample.size() - 1 : j - 1;
            int diag_idx = (j == 0) ? sample.size() - 1 : j - 1;
            int left_idx = j;

            int left_value = score[left_idx].first;
            int diag_value = score[diag_idx].first;
            int up_value = next_score[up_idx].first;

            left_value -= gap;
            up_value = (j==0)?-1: up_value - gap;
            diag_value = (consensus[i]==sample[j]) ? diag_value + match: diag_value - mismatch;

            if (diag_value >= left_value)
            {
                if (diag_value >= up_value)
                {
                    next_score[j] = make_pair(diag_value ,(score[diag_idx].second == -1) ? j:score[diag_idx].second);
                }
                else
                {
                   next_score[j] = make_pair(up_value ,(j==0)?-1:next_score[up_idx].second); 
                }
            }
            else
            {
                if (left_value >= up_value)
                {
                    next_score[j] = make_pair(left_value,score[left_idx].second);
                }
                else
                {
                    next_score[j] = make_pair(up_value ,(j==0)?-1:next_score[up_idx].second);
                }
            }

            if (next_score[j].first > max_.first)
            {
                max_.first = next_score[j].first;
                max_.second = next_score[j].second;
            }
        }

        for (auto z = 0; z < sample.size(); z++)
        {
            score[z] = next_score[z];
        }
    }

    
    // cout << max_.first << " " << max_.second << " ";

    return max_;

}

std::vector<std::string>rotate_sample(std::vector<std::string> consensus, std::vector<std::string> sample, std::unordered_map<std::string, int> &blockSizeMap, int &rotation_index, bool &invert)
{
    std::vector<std::string> rotated_sample;
    std::pair<int,int> front_rotate = rotate_alignment(consensus, sample);

    reverse(sample.begin(), sample.end());
    
    std::pair<int,int> back_rotate = rotate_alignment(consensus, sample);

    int rotate;
    invert = back_rotate.first > front_rotate.first;
    if(!invert)
    {
        reverse(sample.begin(), sample.end());
        rotate = front_rotate.second;
        cout << "Front\n";
    }
    else {
        rotate = back_rotate.second;
        cout << "Back\n";

    }

    // for (auto i=0;i<rotate;i++){
    //     rotation_index += blockSizeMap[sample[i]];
    // }
    rotation_index = rotate;
    int index;
    for (auto i=0;i < sample.size(); i++)
    {
        index = (i+rotate)%sample.size();
        rotated_sample.push_back(sample[index]);
        // output_file << sample[index];
    }

    return rotated_sample;
}


// int main(int argc, char* argv[]) {

//     std::ifstream consensus_file (argv[1]);
//     std::ifstream sample_file (argv[2]);

//     std::vector<std::string> consensus = {};
//     std::string line, colname;
//     int val;

//     std::vector<int> intSequenceConsensus={};
//     std::vector<int> intSequenceSample={};
//     std::vector<int> intSequenceConsensus_new={};


//     if(consensus_file.good())
//     {
//         std::getline(consensus_file, line);

//         std::stringstream ss(line);

//         while(std::getline(ss, colname, ',')){  
//             if (colname == "\n")
//             {
//                 break;
//             }
//             consensus.push_back(colname);
//         }
//     }
//     consensus_file.close();

//     std::vector<std::string> sample = {};

//     if(sample_file.good())
//     {
//         std::getline(sample_file, line);

//         std::stringstream ss(line);

//         while(std::getline(ss, colname, ',')){  
//             if (colname == "\n")
//             {
//                 break;
//             }
//             sample.push_back(colname);
//         }
//     }
//     sample_file.close();    
//     // Sample points

//     size_t numBlocks = 0;
//     std::unordered_map<int,std::string> intToString;
//     // std::unordered_map<string,int> stringToint;
//     std::unordered_map< std::string, std::vector< int > > intSequences; 
//     std::vector<std::string> consensus_new;
    
//     std::string consensus_name = "consensus";
//     std::string sample_name = "sample";
//     intSequences[consensus_name]={};
//     for (auto s: consensus)
//     {
//         intToString[numBlocks] = s;
//         intSequenceConsensus.push_back(numBlocks);
//         intSequences[consensus_name].push_back(numBlocks);
//         // stringToint[s] += 1;
//         numBlocks++;

//     }
    
//     std::vector<std::string> rotated_sample;
//     rotated_sample = rotate_sample(consensus, sample);

//     // cout << rotate << " " << sample.size() << endl;

    
//     return 0;
// }