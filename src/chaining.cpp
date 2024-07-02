#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <limits>
#include <chrono>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/concurrent_vector.h>

std::pair <int,int> origin (-1,-1);
std::pair <int,int> base (0,0);


using namespace std;



struct hashPair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        auto hash1 = hash<T1> {}(p.first);
        auto hash2 = hash<T2> {}(p.second);

        if (hash1 != hash2) {
            return hash1 ^ hash2;
        }

        // If hash1 == hash2, their XOR is zero.
        return hash1;
    }
};
// Structure to represent a node in the range tree
struct Node {
    std::pair<int,int> point;
    Node* left;
    Node* right;
    Node* parent;
    int score;
};

Node* createNode(std::pair<int,int> point) {
    Node* newNode = new Node;
    newNode->point = point;
    newNode->left = newNode->right = nullptr;
    return newNode;
}

// Function to compare points based on x-coordinate
bool compareX(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    return a.first < b.first;
}

bool comparePoint(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    if (a.first == b.first)
        return (a.second < b.second);
    return a.first < b.first;
}

// Function to compare points based on y-coordinate
bool compareY(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    return a.second < b.second;
}


Node* constructRangeTree(tbb::concurrent_vector<std::pair<int,int>>& points, int start, int end) {
    if (start > end)
        return nullptr;

    sort(points.begin() + start, points.begin() + end + 1, compareX);

    int mid = (start + end) / 2;
    Node* root = createNode(points[mid]);

    root->left = constructRangeTree(points, start, mid - 1);
    root->right = constructRangeTree(points, mid + 1, end);

    return root;
}



// Function to perform range query on the 2D range tree
void queryRange(Node* root, std::pair<int,int> rangeStart, std::pair<int,int> rangeEnd, vector<std::pair<int,int>>& result) {
    if (root == nullptr)
        return;

    if (root->point.first >= rangeStart.first && root->point.first <= rangeEnd.first &&
            root->point.second >= rangeStart.second && root->point.second <= rangeEnd.second) {
        result.push_back(root->point);
    }

    if (root->left != nullptr && rangeStart.first <= root->point.first)
        queryRange(root->left, rangeStart, rangeEnd, result);

    if (root->right != nullptr && rangeEnd.first >= root->point.first)
        queryRange(root->right, rangeStart, rangeEnd, result);
}

void find_chain(Node* root, std::pair<int,int> point, std::unordered_map<std::pair<int,int>, std::pair<int, std::pair<int,int>>,hashPair>&map, int K, pair<int,int> &curr_base,pair<int,int> &max_score_point) {

    std::vector<std::pair<int,int>> result;
    std::pair<int,int> new_base ((point.first - K > 0 ? point.first - K: 0), (point.second - K > 0 ? point.second - K: 0));
    std::pair<int,int> new_point (point.first - 1, point.second - 1);
    int match = 50;
    if (point.first == 0 and point.second == 0) {
        std::pair<int,int> p(-1,-1);
        std::pair<int,std::pair<int,int>> temp_map (match, p);
        map[point] = temp_map;
        return;
    }
    queryRange(root, new_base, new_point, result);
    // std::cout << new_base.first << " " << new_base.second << " " << point.first << " " << point.second  << " " << result.size()  << endl;

    int temp_score = 10;
    std::pair<int,int> temp_node = origin;

    // Barrier
    int x_b = -1;
    int y_b = -1;
    for (vector<std::pair<int,int>>::reverse_iterator i = result.rbegin(); i != result.rend(); ++i ) {
        std::pair<int,int> p = *i;

        if (p.first <= x_b and p.second <= y_b) {
            continue;
        }
        // int cost = (std::min(point.first-p.first, point.second-p.second)) + (std::max(point.second-p.second,point.first-p.first));
        int cost = -(point.first - p.first + point.second - p.second);
        if (cost + map[p].first + match > temp_score) {
            temp_score = cost + map[p].first + match ;
            temp_node = p;
        }
        if (x_b < p.first) {
            x_b = p.first - 1;
        }
        if (y_b < p.second) {
            y_b = p.second - 1;
        }

    }

    std::pair<int,std::pair<int,int>> temp_map (temp_score, temp_node);
    map[point] = temp_map;

}


std::vector<std::pair<int,int>> chaining (std::vector<std::string> &consensus, std::vector<std::string> &sample) {
    std::vector<std::pair<int,int>> chain;
    int K = 500;

    std::vector<std::pair<int,int>>  points;
    // std::cout << "Finding seeds sequencial ";
    auto start = std::chrono::high_resolution_clock::now();

    // Concurrent Vector
    tbb::concurrent_vector<std::pair<int,int>> points_conc;
    tbb::parallel_for((size_t)0, consensus.size(),[&](size_t i) {
        tbb::parallel_for((size_t)0, sample.size(),[&](size_t j) {
            if (consensus[i]==sample[j]) {
                std::pair<int,int> p (i,j);
                points_conc.push_back(p);
            }
        });
    });


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds time_ = (end-start);
    // std::cout << time_.count() << "\n";

    // Sort points by x-coordinate
    std::vector<std::pair<int,int>> pointsX = points;
    // cout << "Sorting ";
    start = std::chrono::high_resolution_clock::now();
    tbb::parallel_sort(points_conc.begin(), points_conc.end(), comparePoint);
    end = std::chrono::high_resolution_clock::now();
    time_ = (end-start);
    // cout << time_.count() << "\n";


    // Constructing the 2D range tree
    // std::cout << "Range Tree Construction Function ";
    start = std::chrono::high_resolution_clock::now();
    Node* root = constructRangeTree(points_conc, 0, points_conc.size() - 1);
    end = std::chrono::high_resolution_clock::now();
    time_ = (end-start);
    // cout << time_.count() << "\n";

    if (points_conc.size() == 0) {
        return chain;
    }

    std::unordered_map<std::pair<int,int>, std::pair<int,std::pair<int,int>>, hashPair> map;
    for (auto point: points_conc) {
        std::pair<int,std::pair<int,int>> h (-1, origin);
        map[point] = h;
    }

    pair<int,int> max_score_point = base;
    pair<int,int> curr_base = base;

    for (auto point: points_conc) {
        find_chain(root, point, map, K, curr_base, max_score_point);
    }

    int max_score = -1;
    std::pair<int,int> max_score_seed = {};
    for (auto m: map) {
        if (m.second.first>max_score) {
            max_score = m.second.first;
            max_score_seed = m.first;
        }
        // cout << "Score of (" << m.first.first << "," << m.first.second << ") = " << m.second.first << std::endl;
    }

    while (true) {
        chain.push_back(max_score_seed);
        // cout << max_score_seed.first << "," << max_score_seed.second << endl;
        max_score_seed = map[max_score_seed].second;
        if (max_score_seed == origin) {
            break;
        }
    }

    return chain;
}


void build_consensus (
    std::vector<std::pair<int,int>> &chain,
    std::vector<std::string> &consensus,
    std::vector<std::string> &sample,
    std::vector<int> &intSequenceConsensus,
    std::vector<int> &intSequenceSample,
    size_t &numBlocks,
    std::vector<std::string> &consensus_new,
    std::vector<int> &intSequenceConsensus_new,
    std::unordered_map<int,std::string> &intToString
) {

    int prev_consensus_coord = -1;
    int prev_sample_coord = -1;
    for (vector<std::pair<int,int>>::reverse_iterator i = chain.rbegin(); i != chain.rend(); ++i ) {
        int consensus_coord = i->first;
        int sample_coord = i->second;

        for (auto j = prev_consensus_coord + 1; j < consensus_coord; ++j) {
            consensus_new.push_back(consensus[j]);
            intSequenceConsensus_new.push_back(intSequenceConsensus[j]);
        }
        for (auto j = prev_sample_coord + 1; j < sample_coord; ++j) {
            consensus_new.push_back(sample[j]);
            intSequenceSample.push_back(numBlocks);
            intToString[numBlocks] = sample[j];
            intSequenceConsensus_new.push_back(numBlocks);
            numBlocks++;
        }
        consensus_new.push_back(consensus[consensus_coord]);
        intSequenceSample.push_back(intSequenceConsensus[consensus_coord]);
        intSequenceConsensus_new.push_back(intSequenceConsensus[consensus_coord]);
        prev_consensus_coord = consensus_coord;
        prev_sample_coord = sample_coord;
    }

    for (auto j = prev_consensus_coord + 1; j < (int)consensus.size(); ++j) {
        consensus_new.push_back(consensus[j]);
        intSequenceConsensus_new.push_back(intSequenceConsensus[j]);
    }

    for (auto j = prev_sample_coord + 1; j < (int)sample.size(); ++j) {
        consensus_new.push_back(sample[j]);
        intSequenceSample.push_back(numBlocks);
        intToString[numBlocks] = sample[j];
        intSequenceConsensus_new.push_back(numBlocks);
        numBlocks++;
    }
}

void chain_align (
    std::vector<std::string> &consensus,
    std::vector<std::string> &sample,
    std::vector<int> &intSequenceConsensus,
    std::vector<int> &intSequenceSample,
    size_t &numBlocks,
    std::vector<std::string> &consensus_new,
    std::vector<int> &intSequenceConsensus_new,
    std::unordered_map<int,std::string> &intToString
) {
    // cout << "Entering chain Alignment Function\n";
    // if (consensus.size() != 0)
    {
        std::vector<std::pair<int,int>> chain = chaining(consensus, sample);
        build_consensus (chain,
                         consensus,
                         sample,
                         intSequenceConsensus,
                         intSequenceSample,
                         numBlocks,
                         consensus_new,
                         intSequenceConsensus_new,
                         intToString);
    }

}
