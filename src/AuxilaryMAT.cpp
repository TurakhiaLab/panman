#include "AuxilaryMAT.hpp"
#include <vector>
#include <tuple>

std::string stripGaps(const std::string& s){
    std::string t;
    for(auto c: s){
        if(c != '-'){
            t+=c;
        }
    }
    return t;
}

void AuxilaryMAT::Tree::printFASTAHelper(AuxilaryMAT::Node* currentNode, std::string& sequence, std::ofstream& fout, size_t lineLength){
    // pos, oldVal, newVal
    std::vector< std::tuple<int, char, char> > appliedMutations;
    for(auto sub: currentNode->substitutions){
        int position = sub.position;
        char oldChar = sequence[position];
        char newChar = sub.nuc;
        sequence[position] = newChar;
        appliedMutations.push_back(std::make_tuple(position, oldChar, newChar));
    }

    if(currentNode->children.size() == 0){
        fout << '>' << (currentNode->identifier) << '\n';
        std::string sequenceToPrint = stripGaps(sequence);
        for(size_t i = 0; i < sequenceToPrint.length(); i+=lineLength){
            fout << sequenceToPrint.substr(i, std::min(lineLength, sequenceToPrint.length() - i)) << '\n';
        }
    } else {
        for(auto child: currentNode->children){
            printFASTAHelper(child, sequence, fout);
        }
    }

    for(auto mut: appliedMutations){
        sequence[std::get<0>(mut)] = std::get<1>(mut);
    }
}

void AuxilaryMAT::Tree::printFASTA(std::ofstream& fout){
    std::string sequence(consensusSeqLength, '-');
    printFASTAHelper(root, sequence, fout);
}