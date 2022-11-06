#include <iostream>
#include <chrono>
#include <filesystem>

#include <fstream>

#include "PangenomeMAT.hpp"

std::vector< std::string > splitString(const std::string& s){
    std::vector< std::string > res;
    std::string current;

    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            current += s[i];
        } else {
            if(current.length()){
                res.push_back(current);
                current = "";
            }
        }
    }

    if(current.length()){
        res.push_back(current);
        current = "";
    }

    return res;
}

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cout << "Please provide file name.\n";
        return -1;
    }

    try {

        if (std::string(argv[1]) == "--help") {
            std::cout << "./panmat-utils [FILENAME]\n";
            return 0;
        }

        std::ifstream input(argv[1]);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        
        PangenomeMAT::Tree T(input);

        auto treeBuiltEnd= std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << '\n';
        
        // T.printBfs();

        while(true){
            std::cout << "> ";

            std::string command;
            std::getline (std::cin, command);

            std::vector< std::string > splitCommand = splitString(command);

            if(splitCommand.size() == 1 && splitCommand[0] == "summary"){
                auto summaryStart = std::chrono::high_resolution_clock::now();
                T.printSummary();
                auto summaryEnd = std::chrono::high_resolution_clock::now();

                std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

                std::cout << "\nSummary creation time: " << summaryTime.count() << '\n';
            } else if(splitCommand.size() >= 2 && splitCommand[0] == "fasta"){
                std::string fileName = splitCommand[1];
                bool aligned = false;

                if(splitCommand.size() == 3 && splitCommand[1] == "--aligned"){
                    aligned = true;
                    fileName = splitCommand[2];
                }

                std::filesystem::create_directory("./fasta");
                std::ofstream fout("./fasta/" + fileName + ".fasta");

                auto fastaStart = std::chrono::high_resolution_clock::now();
                
                T.printFASTA(fout, aligned);

                auto fastaEnd = std::chrono::high_resolution_clock::now();
                
                std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                std::cout << "\nFASTA execution time: " << fastaTime.count() << '\n';

                fout.close();
            } else if(splitCommand.size() == 2 && splitCommand[0] == "write") {
                std::string fileName = splitCommand[1];
                std::filesystem::create_directory("./pmat");
                std::ofstream fout("./pmat/" + fileName + ".pmat");

                auto writeStart = std::chrono::high_resolution_clock::now();
                
                T.writeToFile(fout);

                auto writeEnd = std::chrono::high_resolution_clock::now();
                
                std::chrono::nanoseconds writeTime = writeEnd - writeStart;

                std::cout << "\nTree Write execution time: " << writeTime.count() << '\n';

                fout.close();
            } else if(splitCommand.size() > 2 && splitCommand[0] == "subtree"){
                if(splitCommand[1].substr(0,8) == "--input="){
                    std::string inputFileName = splitCommand[1].substr(8);

                    std::ifstream fin(inputFileName);

                    std::vector< std::string > nodeIds;
                    std::string nodeId;

                    while(fin >> nodeId){
                        nodeIds.push_back(nodeId);
                    }

                    fin.close();

                    std::string outputFileName = splitCommand[2];
                    std::filesystem::create_directory("./pmat");
                    std::ofstream fout("./pmat/" + outputFileName + ".pmat");

                    auto subtreeStart = std::chrono::high_resolution_clock::now();

                    T.writeToFile(fout, T.subtreeExtractParallel(nodeIds));

                    auto subtreeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << '\n';

                    fout.close();

                } else {
                    std::string fileName = splitCommand[1];
                    std::filesystem::create_directory("./pmat");
                    std::ofstream fout("./pmat/" + fileName + ".pmat");

                    std::vector< std::string > nodeIds;
                    for(size_t i = 2; i < splitCommand.size(); i++){
                        nodeIds.push_back(splitCommand[i]);
                    }
                    auto subtreeStart = std::chrono::high_resolution_clock::now();

                    T.writeToFile(fout, T.subtreeExtractParallel(nodeIds));

                    auto subtreeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << '\n';

                    fout.close();
                }
            } else if(splitCommand.size() > 2 && splitCommand[0] == "subtree-newick"){

                if(splitCommand[1].substr(0,8) == "--input="){
                    std::string inputFileName = splitCommand[1].substr(8);

                    std::ifstream fin(inputFileName);

                    std::vector< std::string > nodeIds;
                    std::string nodeId;

                    while(fin >> nodeId){
                        nodeIds.push_back(nodeId);
                    }

                    fin.close();

                    std::string outputFileName = splitCommand[2];
                    std::filesystem::create_directory("./newick");
                    std::ofstream fout("./newick/" + outputFileName + ".newick");

                    auto subtreeStart = std::chrono::high_resolution_clock::now();

                    fout << T.getNewickString(T.subtreeExtractParallel(nodeIds));

                    auto subtreeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << '\n';

                    fout.close();

                } else {
                    std::string fileName = splitCommand[1];
                    std::filesystem::create_directory("./newick");
                    std::ofstream fout("./newick/" + fileName + ".newick");

                    std::vector< std::string > nodeIds;
                    for(size_t i = 2; i < splitCommand.size(); i++){
                        nodeIds.push_back(splitCommand[i]);
                    }
                    auto subtreeStart = std::chrono::high_resolution_clock::now();

                    fout << T.getNewickString(T.subtreeExtractParallel(nodeIds));

                    auto subtreeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << '\n';

                    fout.close();
                }
            } else if(splitCommand.size() == 1 && splitCommand[0] == "newick"){
                std::cout << T.getNewickString(T.root) << std::endl;
            } else if(splitCommand.size() == 1 && splitCommand[0] == "exit"){
                return 0;
            }
        }

        // T.printBfs();

    } catch(std::exception& e){
        std::cout << e.what() << std::endl;
        return -1;
    }

}