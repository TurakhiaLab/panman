#include <iostream>
#include <algorithm>
#include <stack>
#include "panmanUtils.cuh"
#include "fitchSankoff.cuh"
#ifndef UTILS
#include "utils.hpp"
#endif
panmanUtils::Node::Node(std::string id, float len) {
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

panmanUtils::Node::Node(std::string id, Node* par, float len) {
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

void stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0, temp_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        std::string sub;
        if (temp_pos == 0) {
            sub = s.substr(start_pos, end_pos-start_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 1) {
                temp_pos = start_pos;
            }
            else {
                words.emplace_back(sub);
            }
        }
        else {
            sub = s.substr(temp_pos, end_pos-temp_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 0) {
                temp_pos = 0;
                words.emplace_back(sub);
            }
        }
        // words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}


std::string stripString(std::string s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            return s.substr(i);
        }
    }
    return s;
}

panmanUtils::Node* panmanUtils::Tree::createTreeFromNewickString(std::string newickString) {
    newickString = stripString(newickString);

    Node* treeRoot = nullptr;

    std::vector<std::string> leaves;
    std::vector<size_t> numOpen;
    std::vector<size_t> numClose;
    std::vector<std::queue<float>> branchLen (128);  // will be resized later if needed
    size_t level = 0;
    // std::cout << newickString << std::endl;

    std::vector<std::string> s1;
    stringSplit(newickString, ',', s1);

    numOpen.reserve(s1.size());
    numClose.reserve(s1.size());
    

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        size_t leafDepth = 0;

        bool stop = false;
        bool branchStart = false;
        bool nameZone = false;
        bool hasApo = false;
        std::string leaf = "";
        std::string branch = "";

        for (auto c: s) {
            if (nameZone) {
                leaf += c;
                if (c == '\'') nameZone = false;
            } else if (c == '\'' && !nameZone) {
                nameZone = true;
                hasApo = true;
                leaf += c;
            } else if (c == ':') {
                stop = true;
                branch = "";
                branchStart = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branchLen.size() <= level) {
                    branchLen.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                // float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                float len = (branch.size() > 0) ? std::stof(branch) : 1.0;
                if (len == 0) len = 1.0;
                branchLen[level].push(len);
                level--;
                branchStart = false;
            } else if (!stop) {
                leaf += c;
                branchStart = false;
                leafDepth = level;

            } else if (branchStart) {
                if (isdigit(c)  || c == '.') {
                    branch += c;
                }
            }
        }
        if (hasApo && leaf[0] == '\'' && leaf[leaf.length()-1] == '\'') leaf = leaf.substr(1, leaf.length()-2);
        leaves.push_back(std::move(leaf));
        numOpen.push_back(no);
        numClose.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : 1.0;
        if (len == 0) len = 1.0;
        branchLen[level].push(len);

        // Adjusting max and mean depths
        m_maxDepth = std::max(m_maxDepth, leafDepth);
        m_meanDepth += leafDepth;

    }


    m_meanDepth /= leaves.size();

    // std::cout << m_meanDepth << " " << level << std::endl;
    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    m_numLeaves = leaves.size();

    std::stack<Node*> parentStack;
    int cc = 0;
    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = numOpen[i];
        auto nc = numClose[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = newInternalNodeId();
            Node* newNode = nullptr;
            if (parentStack.size() == 0) {
                newNode = new Node(nid, branchLen[level].front());
                treeRoot = newNode;
            } else {
                newNode = new Node(nid, parentStack.top(), branchLen[level].front());
        
            }
            branchLen[level].pop();
            level++;

            if (allNodes.find(nid) != allNodes.end()) {
                fprintf(stderr, "ERROR: Node with id %s already exists!\n", nid.c_str());
            }
            allNodes[nid] = newNode;
            parentStack.push(newNode);
            cc++;
        }
        if (allNodes.find(leaf) != allNodes.end()) {
            fprintf(stderr, "ERROR: Node with id %s already exists!\n", leaf.c_str());
        }
        Node* leafNode = new Node(leaf, parentStack.top(), branchLen[level].front());
        allNodes[leaf] = leafNode;

        branchLen[level].pop();
        for (size_t j=0; j<nc; j++) {
            parentStack.pop();
            level--;
        }
    }

    if (treeRoot == nullptr) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    treeRoot->branchLength = 0.0;
    std::cout << "Tree created with " << m_numLeaves << " leaves and " << allNodes.size() << " nodes\n";
    return treeRoot;
}


panmanUtils::Tree::Tree(std::ifstream& newick_ifstream) {
    std::string newickString;
    std::getline(newick_ifstream, newickString);
    root = createTreeFromNewickString(newickString);
}

int main(int argc, char**argv){
    std::ifstream newick_file_istream (argv[1]); // newick file
    std::string fname = argv[2]; // MSA fasta file
    std::string ref_file = argv[3]; // ref alignment fasta file

    auto start = std::chrono::high_resolution_clock::now();
    
    panmanUtils::Tree* T = new panmanUtils::Tree(newick_file_istream);
if (1) {
    printf("=== Tree information ===\n");
    printf("Number of nodes: %d\n", T->allNodes.size());
    printf("============================\n");
}
    std::unordered_map<std::string, std::string> seqs;
    std::unordered_map<std::string, std::string> refs;
    if(ref_file != "")
        read_seqs(ref_file, refs);

    utility::util *u = new utility::util(10000, 300, 0, 0);
    u->seq_file_name = fname;

    for (auto &a: refs)
    {
        u->ref_name = a.first;
        u->ref_seq = a.second;
        u->msa_len = a.second.size();
        u->consensus = a.second;
    }

    fitch_sankoff_on_gpu(T, seqs, u);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds total = end - start;
    printf("Total time: %lf mins\n", ((double)total.count()/1000000000)/60);
    return 0;
}
