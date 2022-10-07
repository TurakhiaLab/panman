#include <iostream>
#include <string>
#include <vector>
#include "PangenomeMAT.hpp"
#include "mutation_annotation.pb.h"

PangenomeMAT::Node::Node(std::string id, float len){
	identifier = id;
	level = 1;
	branchLength = len;
	parent = NULL;
}

PangenomeMAT::Node::Node(std::string id, Node* par, float len){
	identifier = id;
	parent = par;
	level = par->level + 1;
	branchLength = len;
}


void PangenomeMAT::stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    // TIMEIT();
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        // if ((end_pos == start_pos) || end_pos >= s.length()) {
        if (end_pos >= s.length()) {
            break;
        }
        words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}

PangenomeMAT::Node* PangenomeMAT::Tree::createTreeFromNewickString(std::string newick_string) {
    // TIMEIT();

    PangenomeMAT::Node* newTreeRoot;

    std::vector<std::string> leaves;
    std::vector<size_t> num_open;
    std::vector<size_t> num_close;
    std::vector<std::queue<float>> branch_len (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    stringSplit(newick_string, ',', s1);

    num_open.reserve(s1.size());
    num_close.reserve(s1.size());

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        bool stop = false;
        bool branch_start = false;
        std::string leaf = "";
        std::string branch = "";
        for (auto c: s) {
            if (c == ':') {
                stop = true;
                branch = "";
                branch_start = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branch_len.size() <= level) {
                    branch_len.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                branch_len[level].push(len);
                level--;
                branch_start = false;
            } else if (!stop) {
                leaf += c;
                branch_start = false;
            } else if (branch_start) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        leaves.push_back(std::move(leaf));
        num_open.push_back(no);
        num_close.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        branch_len[level].push(len);
    }

    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    std::stack<Node*> parent_stack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = num_open[i];
        auto nc = num_close[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = newInternalNodeId();
            Node* new_node = NULL;
            if (parent_stack.size() == 0) {
                newTreeRoot = new Node(nid, branch_len[level].front());
            } else {
                new_node = new Node(nid, parent_stack.top(), branch_len[level].front());
                (parent_stack.top()->children).push_back(new_node);

            }
            branch_len[level].pop();
            level++;
            parent_stack.push(new_node);
        }
        Node* tempNode = new Node(leaf, parent_stack.top(), branch_len[level].front());
        parent_stack.top()->children.push_back(tempNode);
        branch_len[level].pop();
        for (size_t j=0; j<nc; j++) {
            parent_stack.pop();
            level--;
        }
    }

    if (newTreeRoot == NULL) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    return newTreeRoot;
}

PangenomeMAT::Tree::Tree(std::ifstream& fin){
	MAT::tree mainTree;

	if(!mainTree.ParseFromIstream(&fin)){
		std::cerr << "Could not read tree from input file." << std::endl;
		return;
	}

	root = createTreeFromNewickString(mainTree.newick());

	// Traversal test
	std::queue<Node *> bfsQueue;
	int prevLev = 0;
	bfsQueue.push(root);
	while(!bfsQueue.empty()){
		Node* current = bfsQueue.front();
		bfsQueue.pop();

		if(current->level != prevLev){
			std::cout << '\n';
			prevLev = current->level;
		}
		std::cout << '(' << current->identifier << "," << current->branchLength << ") ";

		for(auto child: current->children){
			bfsQueue.push(child);
		}
	}
}

int main(int argc, char* argv[]){
	std::ifstream input(argv[1]);

	PangenomeMAT::Tree T(input);
}