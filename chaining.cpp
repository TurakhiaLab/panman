#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <limits>

std::pair <int,int> origin (-1,-1);
std::pair <int,int> base (0,0);


using namespace std;

// Structure to represent a point in 2D space
// struct Point {
//     int x, y;
// };

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
 
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

// Function to compare points based on x-coordinate
bool compareX(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    return a.first < b.first;
}

bool compareXG(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    return a.first > b.first;
}

// Function to compare points based on y-coordinate
bool compareY(const std::pair<int,int>& a, const std::pair<int,int>& b) {
    return a.second < b.second;
}

// Function to construct the 2D range tree
Node* constructRangeTree(vector<std::pair<int,int>>& points, int start, int end, bool isXLevel, Node* p) {
    if (start > end) {
        return nullptr;
    }

    // Sort the points based on the current level (x or y)
    if (isXLevel) {
        sort(points.begin() + start, points.begin() + end + 1, compareX);
    } else {
        sort(points.begin() + start, points.begin() + end + 1, compareY);
    }

    // Find the middle point
    int mid = (start + end) / 2;

    // Create a new node with the middle point
    Node* newNode = new Node;
    newNode->point = points[mid];
    newNode->parent = p;
    newNode->left = constructRangeTree(points, start, mid - 1, !isXLevel, newNode);
    newNode->right = constructRangeTree(points, mid + 1, end, !isXLevel, newNode);
    
    return newNode;
}

// Function to perform range query on the 2D range tree
void rangeQuery(Node* root, std::pair<int,int>& queryPoint1, std::pair<int,int>& queryPoint2, vector<std::pair<int,int>>& result) {
    if (root == nullptr) {
        return;
    }

    // Check if the current node lies within the given range
    if (queryPoint1.first <= root->point.first && root->point.first < queryPoint2.first &&
        queryPoint1.second <= root->point.second && root->point.second < queryPoint2.second) {
        result.push_back(root->point);
    }

    // Check if the current node can have points within the given range
    if ((root->left != nullptr) && (queryPoint1.first <= root->point.first)) {
        rangeQuery(root->left, queryPoint1, queryPoint2, result);
    }

    if ((root->right != nullptr) && (root->point.first <= queryPoint2.first)) {
        rangeQuery(root->right, queryPoint1, queryPoint2, result);
    }
}


// void find_range_start(Node* root, std::pair<int,int> point)
// {
//     if ()
// }

void find_chain(Node* root, std::pair<int,int> point, std::unordered_map<std::pair<int,int>, std::pair<int, std::pair<int,int>>,hash_pair>&map)
{
    int INF = std::numeric_limits<int>::max();

    std::vector<std::pair<int,int>> result;
    rangeQuery(root, base, point, result);
    // sort(result.begin(), result.begin(), compareXG);

    std::unordered_map<std::pair<int,int>, bool, hash_pair> visited;
    int temp_score = 10;
    int match = 50;
    std::pair<int,int> temp_node = origin; 

    int alpha = -1;
    int beta = 1;

    // Barrier
    int x_b = 0;
    int y_b = 0;
    for (vector<std::pair<int,int>>::reverse_iterator i = result.rbegin(); i != result.rend(); ++i ) 
    {
        std::pair<int,int> p = *i;
        if (p.first<x_b and p.second<y_b)
        {
            continue;
        }
        // if (visited[p])
        // {
        //     continue;
        // }
        visited[p]=1;
        // int cost = alpha*(abs(abs(point.second-p.second)-abs(point.first-p.first))) + beta*(std::min(point.first-p.first, point.second-p.second));
        int cost = alpha*(std::max(point.second-p.second,point.first-p.first)) + beta*(std::min(point.first-p.first, point.second-p.second));
        // int cost = abs(abs(point.second-p.second)+abs(point.first-p.first));
        // cout << cost << " " << map[p].first<< endl;
        if (cost + map[p].first + match > temp_score)
        {
            // cout << "Max score found at:" << p.first << "," << p.second << "=" << cost + map[p].first + match<< endl;
            temp_score = cost + map[p].first + match ;
            temp_node = p;
        }

        std::vector<std::pair<int,int>> inter;
        if (x_b < p.first)
        {
            x_b = p.first;
        }
        if (y_b < p.second)
        {
            y_b = p.second;
        }
        // rangeQuery(root, base, p, inter);
        // for (auto t: inter)
        // {
        //     // cout << "visting: " << t.first << "," << t.second << endl;
        //     visited[t] = 1;
        // }
    }
    
    std::pair<int,std::pair<int,int>> temp_map (temp_score, temp_node);
    map[point] = temp_map;

}

std::vector<std::pair<int,int>> chaining (std::vector<std::string> &consensus, std::vector<std::string> &sample)
{
    cout << "Entering chaining Function\n";
    std::vector<std::pair<int,int>> chain;
    int INF = std::numeric_limits<int>::max();

    // std::vector<Point> points = {Point(2, 3), Point(5, 7), Point(8, 1), Point(10, 12)};
    std::vector<std::pair<int,int>> points = {};

    for (int i =0; i<consensus.size(); i++)
    {
        for (int j = 0; j < sample.size(); j++)
        {
            if (consensus[i]==sample[j])
            {
                std::pair<int,int> p = {i,j};
                points.push_back(p);
            }
        }
    }

    // Constructing the 2D range tree
    cout << "Entering Range Tree Construction Function\n";
    Node* root = constructRangeTree(points, 0, points.size() - 1, true, nullptr);
    cout << "Exiting Range Tree Construction Function\n";
    // Sort points by x-coordinate
    std::vector<std::pair<int,int>> pointsX = points;
    std::sort(pointsX.begin(), pointsX.end(), compareX);

    std::unordered_map<std::pair<int,int>, std::pair<int,std::pair<int,int>>, hash_pair> map;
    for (auto point: pointsX)
    {

        std::pair<int,std::pair<int,int>> h (-1, origin);
        map[point] = h;
    }
    int count = 0;
    cout << "total Seeds:" << pointsX.size() << endl;
    for (auto point: pointsX)
    {   
        cout << (count++) << endl;
        find_chain(root, point, map);
        // cout << "Chained:" << point.first << "," << point.second << " with: " << map[point].second.first << "," << map[point].second.second << endl;
    } 
    
    int max_score = -1;
    std::pair<int,int> max_score_seed = {};
    for (auto m: map)
    {
        if (m.second.first>max_score)
        {
            max_score = m.second.first;
            max_score_seed = m.first;
        }
        // cout << "Score of (" << m.first.first << "," << m.first.second << ") = " << m.second.first << std::endl;
    }

    while (true)
    {
        chain.push_back(max_score_seed);
        cout << max_score_seed.first << "," << max_score_seed.second << endl;
        max_score_seed = map[max_score_seed].second;
        if (max_score_seed == origin)
        {
            break;
        }
    }
    cout << "Exiting chaining Function\n";
    return chain;
}


void build_consensus (
    std::vector<std::pair<int,int>> &chain,
    std::vector<std::string> &consensus, 
    std::vector<std::string> &sample,
    std::string consensus_name, 
    std::string sample_name, 
    std::unordered_map< std::string, std::vector< int > > &intSequences,
    int &numBlocks,
    std::vector<std::string> &consensus_new 
)
{
    

    int prev_consensus_coord = -1;    
    int prev_sample_coord = -1;
    for (vector<std::pair<int,int>>::reverse_iterator i = chain.rbegin(); i != chain.rend(); ++i ) 
    {
        int consensus_coord = i->first;    
        int sample_coord = i->second;

        for (auto j = prev_consensus_coord + 1; j < consensus_coord; ++j)
        {
            consensus_new.push_back(consensus[j]);
        } 
        for (auto j = prev_sample_coord + 1; j < sample_coord; ++j)
        {
            consensus_new.push_back(sample[j]);
            intSequences[sample_name].push_back(numBlocks);
            numBlocks++;
        }    
        consensus_new.push_back(consensus[consensus_coord]);
        intSequences[sample_name].push_back(intSequences[consensus_name][consensus_coord]);
        prev_consensus_coord = consensus_coord;
        prev_sample_coord = sample_coord;
    } 

    for (auto j = prev_consensus_coord + 1; j < consensus.size(); ++j)
    {
        consensus_new.push_back(consensus[j]);
    } 

    for (auto j = prev_sample_coord + 1; j < sample.size(); ++j)
    {
        consensus_new.push_back(sample[j]);
        intSequences[sample_name].push_back(numBlocks);
        numBlocks++;
    }
}

void chain_align (
    std::vector<std::string> &consensus, 
    std::vector<std::string> &sample,
    std::string consensus_name, 
    std::string sample_name, 
    std::unordered_map< std::string, std::vector< int > > &intSequences,
    int &numBlocks,
    std::vector<std::string> &consensus_new,
    std::unordered_map<int,std::string> &intToString 
)
{
    cout << "Entering chain Alignment Function\n";
    if (consensus.size() != 0)
    {
        std::vector<std::pair<int,int>> chain = chaining(consensus, sample);
        build_consensus (chain, 
                    consensus, 
                    sample, 
                    consensus_name, 
                    sample_name, 
                    intSequences, 
                    numBlocks, 
                    consensus_new); 
    }
    else {
        consensus_new = sample;
        intSequences[sample_name]={};
        for (auto s: sample)
        {
            intToString[numBlocks] = s;
            intSequences[sample_name].push_back(numBlocks);
            numBlocks++;
        }
    }
}

int main(int argc, char* argv[]) {

    std::ifstream consensus_file (argv[1]);
    std::ifstream sample_file (argv[2]);

    std::vector<std::string> consensus = {};
    std::string line, colname;
    int val;

    if(consensus_file.good())
    {
        std::getline(consensus_file, line);

        std::stringstream ss(line);

        while(std::getline(ss, colname, ',')){  
            if (colname == "\n")
            {
                break;
            }
            consensus.push_back(colname);
        }
    }
    consensus_file.close();

    std::vector<std::string> sample = {};

    if(sample_file.good())
    {
        std::getline(sample_file, line);

        std::stringstream ss(line);

        while(std::getline(ss, colname, ',')){  
            if (colname == "\n")
            {
                break;
            }
            sample.push_back(colname);
        }
    }
    sample_file.close();    
    // Sample points

    int numBlocks = 0;
    std::unordered_map<int,std::string> intToString;
    // std::unordered_map<string,int> stringToint;
    std::unordered_map< std::string, std::vector< int > > intSequences; 
    std::vector<std::string> consensus_new;
    
    string consensus_name = "consensus";
    string sample_name = "sample";
    intSequences[consensus_name]={};
    for (auto s: consensus)
    {
        intToString[numBlocks] = s;
        intSequences[consensus_name].push_back(numBlocks);
        // stringToint[s] += 1;
        numBlocks++;

    }

    chain_align (consensus, 
                sample, 
                consensus_name, 
                sample_name, 
                intSequences, 
                numBlocks, 
                consensus_new,
                intToString); 

    for (auto i = 0; i <  consensus_new.size(); ++i)
    {
        cout << consensus_new[i] << " ";
    }
    cout << "\n";

    for (auto i = 0; i <  intSequences[consensus_name].size(); ++i)
    {
        cout << intSequences[consensus_name][i] << " ";
    }
    cout << "\n";
    
    for (auto i = 0; i <  intSequences[sample_name].size(); ++i)
    {
        cout << intSequences[sample_name][i] << " ";
    }
    cout << "\n";

    cout << numBlocks;


    // // Querying for points within the range (2,2) and (7,7)
    // std::pair<int,int> queryPoint1 = {0, 0};
    // std::pair<int,int> queryPoint2 = {4, 4};
    // vector<std::pair<int,int>> result;
    // rangeQuery(root, queryPoint1, queryPoint2, result);

    // // Printing the result
    // cout << "Points within the range (" << queryPoint1.first << "," << queryPoint1.second << ") and (" << queryPoint2.first << "," << queryPoint2.second << "):" << endl;
    // for (const auto& point : result) {
    //     cout << "(" << point.first << "," << point.second << ")" << endl;
    // }

    return 0;
}