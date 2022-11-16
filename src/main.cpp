#include <iostream>
#include <chrono>
#include <filesystem>

#include <fstream>

#include "PangenomeMAT.hpp"

std::vector< std::string > splitString(const std::string& s, char delim = ' '){
    std::vector< std::string > res;
    std::string current;

    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != delim){
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

void stripString(std::string& s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            s = s.substr(i);
            return;
        }
    }
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
            } else if(splitCommand.size() == 3 && splitCommand[0] == "vcf"){
                if(splitCommand[1].substr(0,12) != "--reference="){
                    std::cout << "Please enter reference sequence ID as arg1. Format: --reference=<id>\n";
                    return 0;
                }
                std::string reference = splitCommand[1].substr(12);
                std::string fileName = splitCommand[2];
                std::filesystem::create_directory("./vcf");
                std::ofstream fout("./vcf/" + fileName + ".vc");

                auto vcfStart = std::chrono::high_resolution_clock::now();

                T.printVCF(reference, fout);

                auto vcfEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;

                std::cout << "\nVCF execution time: " << vcfTime.count() << '\n';

                fout.close();

                // Debugging

                // std::ifstream fin("./vcf/" + fileName + ".vc");

                // T.getSequenceFromVCF("NC_007946.1", fin);

                // fin.close();

            } else if(splitCommand.size() == 1 && splitCommand[0] == "newick"){
                std::cout << T.getNewickString(T.root) << std::endl;
            } else if(splitCommand.size() == 2 && splitCommand[0] == "annotate"){
                std::string fileName = splitCommand[1];
                std::ifstream fin(fileName);

                auto annotateStart = std::chrono::high_resolution_clock::now();

                T.annotate(fin);

                auto annotateEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds annotateTime = annotateEnd - annotateStart;

                std::cout << "Annotate time: " << annotateTime.count() << std::endl;


            } else if(splitCommand.size() > 1 && splitCommand[0] == "search"){
                std::string word;
                std::string restOfCommand;
                for(size_t i = 0; i < command.length(); i++){
                    if(command[i] != ' '){
                        word+=command[i];
                    }
                    if(word == "search"){
                        restOfCommand = command.substr(i+1);
                        break;
                    }
                }

                std::vector< std::string > annotationVector = splitString(restOfCommand, ',');

                for(size_t i = 0; i < annotationVector.size(); i++){
                    stripString(annotationVector[i]);
                    std::cout << annotationVector[i] << ": ";
                    auto result = T.searchByAnnotation(annotationVector[i]);
                    for(auto r: result){
                        std::cout << r << ";";
                    }
                    std::cout << std::endl;
                }
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