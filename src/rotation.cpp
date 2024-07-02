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


std::pair<int, int> rotate_alignment(const std::vector<std::string>& consensus, const std::vector<std::string>& sample) {
    std::vector<std::pair<int, int>> score(sample.size(), std::make_pair(-1,-1));
    std::vector<std::pair<int, int>> next_score(sample.size(), std::make_pair(-1,-1));
    std::pair<int,int>max_(0,0);
    int match = 5;
    int gap = 1;
    int mismatch = 2;
    for(size_t i = 0; i < consensus.size(); i++) {
        for (size_t j = 0; j < sample.size(); j++) {
            int up_idx = (j == 0) ? sample.size() - 1 : j - 1;
            int diag_idx = (j == 0) ? sample.size() - 1 : j - 1;
            int left_idx = j;

            int left_value = score[left_idx].first;
            int diag_value = score[diag_idx].first;
            int up_value = next_score[up_idx].first;

            left_value -= gap;
            up_value = (j==0)?-1: up_value - gap;
            diag_value = (consensus[i]==sample[j]) ? diag_value + match: diag_value - mismatch;

            if (diag_value >= left_value) {
                if (diag_value >= up_value) {
                    next_score[j] = make_pair(diag_value,(score[diag_idx].second == -1) ? j:score[diag_idx].second);
                } else {
                    next_score[j] = make_pair(up_value,(j==0)?-1:next_score[up_idx].second);
                }
            } else {
                if (left_value >= up_value) {
                    next_score[j] = make_pair(left_value,score[left_idx].second);
                } else {
                    next_score[j] = make_pair(up_value,(j==0)?-1:next_score[up_idx].second);
                }
            }

            if (next_score[j].first > max_.first) {
                max_.first = next_score[j].first;
                max_.second = next_score[j].second;
            }
        }

        for (size_t z = 0; z < sample.size(); z++) {
            score[z] = next_score[z];
        }
    }


    cout << max_.first << " " << max_.second << " ";

    return max_;

}

std::vector<std::string>rotate_sample(const std::vector<std::string>& consensus, std::vector<std::string>& sample, std::vector<int>& blockStrand, std::vector< size_t >& blockNumbers, std::unordered_map<std::string, int> &blockSizeMap, int &rotation_index, bool &invert) {
    std::vector<std::string> rotated_sample;
    std::pair<int,int> front_rotate = rotate_alignment(consensus, sample);

    reverse(sample.begin(), sample.end());
    reverse(blockStrand.begin(), blockStrand.end());
    reverse(blockNumbers.begin(), blockNumbers.end());

    // Rotated block index
    int rotate = 0;

#ifdef ALLOW_INVERSIONS
    std::pair<int,int> back_rotate = rotate_alignment(consensus, sample);
    invert = back_rotate.first > front_rotate.first;
    if(!invert) {
#endif
        reverse(sample.begin(), sample.end());
        reverse(blockStrand.begin(), blockStrand.end());
        reverse(blockNumbers.begin(), blockNumbers.end());
        rotate = front_rotate.second;
#ifdef ALLOW_INVERSIONS
        cout << "Front\n";
    } else {
        rotate = back_rotate.second;
        cout << "Back\n";
    }
#endif

    rotation_index = (sample.size() - rotate) % sample.size();
    int index;
    std::vector<int> newBlockStrand = {};
    std::vector<size_t> newBlockNumbers = {};

    for (size_t i = 0; i < sample.size(); i++) {
        index = (i+rotate)%sample.size();
        rotated_sample.push_back(sample[index]);
        newBlockStrand.push_back(blockStrand[index]);
        newBlockNumbers.push_back(blockNumbers[index]);
    }

    blockStrand = newBlockStrand;
    blockNumbers = newBlockNumbers;
    return rotated_sample;
}
