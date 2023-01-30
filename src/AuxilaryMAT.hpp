#ifndef AUXILARY_MAT_HPP
#define AUXILARY_MAT_HPP

#include <string>
#include <vector>
#include <fstream>

namespace AuxilaryMAT {

    struct Substitution {
        uint32_t position;
        char nuc;
    };

    struct Node {
        std::string identifier;
        std::vector< Substitution > substitutions;
        std::vector< Node* > children;
    };

    struct Tree {
        uint32_t consensusSeqLength;
        Node* root;

        void printFASTA(std::ofstream& fout);

        private:
        void printFASTAHelper(Node* currentNode, std::string& sequence, std::ofstream& fout, size_t lineLength=70);

    };

};

#endif // AUXILARY_MAT_HPP