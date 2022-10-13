#include <iostream>
#include "PangenomeMAT.hpp"

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cout << "Please provide file name.\n";
        return -1;
    }

    try {
        std::ifstream input(argv[1]);
        PangenomeMAT::Tree T(input);
        T.printSummary();
    } catch(std::exception& e){
        std::cout << e.what() << std::endl;
        return -1;
    }

}