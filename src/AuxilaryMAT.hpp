#ifndef AUXILARY_MAT_HPP
#define AUXILARY_MAT_HPP

#include <string>
#include <vector>
#include <fstream>

namespace AuxilaryMAT {

    struct Substitution {
        uint32_t position;
        char newNuc;
    };

    struct Node {
        std::string identifier;
        std::vector< Substitution > substitutions;
        std::vector< Node* > children;
    };

    struct Tree {
        std::string consensusSeq;
        Node* root;

        void printFASTA(std::ofstream& fout);

        private:
        void printFASTAHelper(Node* currentNode, std::ofstream& fout);

    };

};

#endif // AUXILARY_MAT_HPP