#define NEW_MAT
#define NEW_PARSER

#include <iostream>
#include <chrono>
#include <filesystem>
#include <boost/program_options.hpp>

#include <fstream>

#ifdef NEW_MAT
#include "PangenomeMATNew.hpp"
#else
#include "PangenomeMAT.hpp"
#endif

namespace po = boost::program_options;

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

po::options_description globalDesc("Pangenome MAT Command Line Arguments");
po::positional_options_description globalPositionArgumentDesc;
po::options_description summaryDesc("Summary Command Line Arguments");
po::options_description fastaDesc("FASTA Command Line Arguments");
po::positional_options_description fastaPositionArgumentDesc;
po::options_description writeDesc("MAT Writer Command Line Arguments");
po::positional_options_description writePositionArgumentDesc;
po::options_description subtreeDesc("Subtree Extract Command Line Arguments");
po::positional_options_description subtreePositionArgumentDesc;
po::options_description vcfDesc("VCF writer Command Line Arguments");
po::positional_options_description vcfPositionArgumentDesc;
po::options_description annotateDesc("MAT Annotate Command Line Arguments");
po::positional_options_description annotatePositionArgumentDesc;
po::options_description searchDesc("Search by annotation Command Line Arguments");
po::positional_options_description searchPositionArgumentDesc;


void setupOptionDescriptions(){
    // Global option descriptions
    globalDesc.add_options()
        ("help", "produce help message")
        ("input-file,I", po::value< std::string >(), "pmat file path")
    ;

    // Adding input file as positional argument
    globalPositionArgumentDesc.add("input-file", -1);

    // FASTA option descriptions
    fastaDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
        ("aligned", "print in aligned format")
        ("num_cores", po::value< size_t >(), "number of cores if parallelized")
    ;

    // Adding output file as positional argument
    fastaPositionArgumentDesc.add("output-file", -1);

    // MAT Writer option descriptions
    writeDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    writePositionArgumentDesc.add("output-file", -1);

    // Subtree Extract option descriptions
    subtreeDesc.add_options()
        ("help", "produce help message")
        ("newick", po::value< bool >()->default_value(false), "just print newick string")
        ("input-file", po::value< std::string >(), "Input file name if reading node IDs from file")
        ("output-file", po::value< std::string >()->required(), "Output file name")
        ("node-ids", po::value< std::vector< std::string > >()->multitoken(), "Node IDs to extract")
    ;

    // Adding output file as positional argument
    subtreePositionArgumentDesc.add("output-file", -1);

    // VCF Writer option descriptions
    vcfDesc.add_options()
        ("help", "produce help message")
        ("reference", po::value< std::string >()->required(), "Sequence ID of the reference sequence")
        ("output-file", po::value< std::string >()->required(), "Output file name")
        ("fasta-file", po::value< std::string >(), "FASTA file name if it should also be created")
    ;

    // Adding output file as positional argument
    vcfPositionArgumentDesc.add("output-file", -1);

    // MAT Annotate option descriptions
    annotateDesc.add_options()
        ("help", "produce help message")
        ("input-file", po::value< std::string >(), "Name of the file containing annotation info")
    ;

    // Adding input file as positional argument
    annotatePositionArgumentDesc.add("input-file", -1);

    // Search by annotation option descriptions
    searchDesc.add_options()
        ("help", "produce help message")
        ("keywords", po::value< std::vector< std::string > >()->multitoken(), "list of keywords to search for")
    ;

    searchPositionArgumentDesc.add("keywords", -1);

}

void printError(std::string e){
    std::cout << "\033[1;31m" << "Error: " << "\033[0m" << e << "\n";
}

void updatedParser(int argc, char* argv[]){

    setupOptionDescriptions();

    po::variables_map globalVm;
    po::store(po::command_line_parser(argc, argv).options(globalDesc).positional(globalPositionArgumentDesc).run(), globalVm);
    po::notify(globalVm);

#ifdef NEW_MAT

    PangenomeMATNew::Tree *T;
    if(globalVm.count("help")){
        std::cout << globalDesc;
        return;
    } else if(globalVm.count("input-file")){
        std::string fileName = globalVm["input-file"].as< std::string >();
        std::ifstream inputStream(fileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new PangenomeMATNew::Tree(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        inputStream.close();
    } else {
        printError("Incorrect Format");
        // std::cout << "\033[1;31m" << "Error: " << "\033[0m" << "Incorrect Format\n";
        std::cout << globalDesc;
        return;
    }

    char** splitCommandArray;

    while(true){
        std::cout << "> ";

        std::string command;
        std::getline (std::cin, command);
        stripString(command);
        std::vector< std::string > splitCommand;
        PangenomeMATNew::stringSplit(command, ' ', splitCommand);

        splitCommandArray = new char*[splitCommand.size()];
        for(size_t i = 0; i < splitCommand.size(); i++){
            splitCommandArray[i] = new char[splitCommand[i].length() + 1];
            strcpy(splitCommandArray[i], splitCommand[i].c_str());
        }
        
        try{
            if(strcmp(splitCommandArray[0], "summary") == 0){
                // If command was summary
                auto summaryStart = std::chrono::high_resolution_clock::now();
                T->printSummary();
                auto summaryEnd = std::chrono::high_resolution_clock::now();

                std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

                std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
            } else if(strcmp(splitCommandArray[0], "fasta") == 0){
                // If command was fasta
                po::variables_map fastaVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(fastaDesc).positional(fastaPositionArgumentDesc).run(), fastaVm);

                if(fastaVm.count("help")){
                    std::cout << fastaDesc;
                } else {
                    po::notify(fastaVm);

                    std::string fileName = fastaVm["output-file"].as< std::string >();

                    bool aligned = false;
                    size_t num_cores = 0;
                    if(fastaVm.count("aligned")){
                        aligned = true;
                    }
                    if(fastaVm.count("num_cores")){
                        num_cores = fastaVm["num_cores"].as< size_t >();
                    }
                    
                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto fastaStart = std::chrono::high_resolution_clock::now();
                    
                    T->printFASTA(fout, aligned, num_cores);

                    auto fastaEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

                    fout.close();

                }
            } else if(strcmp(splitCommandArray[0], "subtree") == 0){
                // If command was subtree
                po::variables_map subtreeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(subtreeDesc).positional(subtreePositionArgumentDesc).run(), subtreeVm);

                if(subtreeVm.count("help")){
                    std::cout << subtreeDesc;
                } else {
                    po::notify(subtreeVm);
                    std::string outputFileName = subtreeVm["output-file"].as< std::string >();
                    std::vector< std::string > nodeIds;
                    std::string nodeId;

                    if(subtreeVm.count("input-file")){
                        std::string inputFileName = subtreeVm["input-file"].as< std::string >();
                        std::ifstream fin(inputFileName);

                        while(fin >> nodeId){
                            nodeIds.push_back(nodeId);
                        }

                        fin.close();
                    } else if(subtreeVm.count("node-ids")){
                        nodeIds = subtreeVm["node-ids"].as< std::vector< std::string > >();
                    } else {
                        printError("No source of node ids provided");
                        std::cout << subtreeDesc;
                    }

                    if(subtreeVm["newick"].as< bool >()){
                        std::filesystem::create_directory("./pmat");
                        std::ofstream fout("./newick/" + outputFileName + ".newick");

                        auto subtreeStart = std::chrono::high_resolution_clock::now();

                        fout << T->getNewickString(T->subtreeExtractParallel(nodeIds));

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";
                        fout.close();
                    } else {
                        std::filesystem::create_directory("./pmat");
                        std::ofstream fout("./pmat/" + outputFileName + ".pmat");

                        auto subtreeStart = std::chrono::high_resolution_clock::now();

                        T->writeToFile(fout, T->subtreeExtractParallel(nodeIds));

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";
                        fout.close();
                    }

                }

            } else if(strcmp(splitCommandArray[0], "vcf") == 0){
                // If command was vcf
                po::variables_map vcfVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(vcfDesc).positional(vcfPositionArgumentDesc).run(), vcfVm);

                if(vcfVm.count("help")){
                    std::cout << vcfDesc;
                } else {
                    po::notify(vcfVm);

                    std::string reference = vcfVm["reference"].as< std::string >();
                    std::string fileName = vcfVm["output-file"].as< std::string >();

                    std::filesystem::create_directory("./vcf");
                    std::ofstream fout("./vcf/" + fileName + ".vc");

                    auto vcfStart = std::chrono::high_resolution_clock::now();

                    T->printVCFParallel(reference, fout);

                    auto vcfEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;

                    std::cout << "\nVCF execution time: " << vcfTime.count() << " nanoseconds\n";

                    fout.close();

                    if(vcfVm.count("fasta-file")){
                        std::cout << "Generating FASTA File" << std::endl;
                        std::ifstream fin("./vcf/" + fileName + ".vc");
                        std::string fastaFileName = vcfVm["fasta-file"].as< std::string >();
                        std::filesystem::create_directory("./fasta");
                        fout.open("./fasta/" + fastaFileName + ".fasta");
                        T->vcfToFASTA(fin, fout);
                        fout.close();
                        fin.close();
                    }

                    // std::cout << "Verifying VCF File" << std::endl;

                    // std::ifstream fin("./vcf/" + fileName + ".vc");

                    // if(T->verifyVCFFile(fin)){
                    //     std::cout << "Success: VCF file generated correctly" << std::endl;
                    // } else {
                    //     std::cout << "Error: VCF extract generated incorrect sequences!" << std::endl;
                    // }

                    // fin.close();
                }

            } else if(strcmp(splitCommandArray[0], "write") == 0){
                // If command was write
                po::variables_map writeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(writeDesc).positional(writePositionArgumentDesc).run(), writeVm);

                if(writeVm.count("help")){
                    std::cout << writeDesc;
                } else {
                    po::notify(writeVm);
                    std::string fileName = writeVm["output-file"].as< std::string >();
                    std::filesystem::create_directory("./pmat");
                    std::ofstream fout("./pmat/" + fileName + ".pmat");

                    auto writeStart = std::chrono::high_resolution_clock::now();
                    
                    T->writeToFile(fout);

                    auto writeEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds writeTime = writeEnd - writeStart;

                    std::cout << "\nTree Write execution time: " << writeTime.count() << " nanoseconds\n";

                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "annotate") == 0){
                // If command was annotate
                po::variables_map annotateVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(annotateDesc).positional(annotatePositionArgumentDesc).run(), annotateVm);

                if(annotateVm.count("help")){
                    std::cout << annotateDesc;
                } else {
                    po::notify(annotateVm);
                    std::string fileName = annotateVm["input-file"].as< std::string >();
                    std::ifstream fin(fileName);

                    auto annotateStart = std::chrono::high_resolution_clock::now();

                    T->annotate(fin);

                    auto annotateEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds annotateTime = annotateEnd - annotateStart;

                    std::cout << "Annotate time: " << annotateTime.count() << " nanoseconds\n";
                }

            } else if(strcmp(splitCommandArray[0], "search") == 0){
                po::variables_map searchVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(searchDesc).positional(searchPositionArgumentDesc).run(), searchVm);

                std::vector< std::string > annotationVector = searchVm["keywords"].as< std::vector< std::string > >();

                for(auto word: annotationVector){
                    stripString(word);
                    std::cout << word << ": ";
                    auto result = T->searchByAnnotation(word);
                    for(auto r: result){
                        std::cout << r << ";";
                    }
                    std::cout << std::endl;
                }

            } else if(strcmp(splitCommandArray[0], "newick") == 0){
                std::cout << T->getNewickString(T->root) << std::endl;
            } else if(strcmp(splitCommandArray[0], "exit") == 0){
                return;
            }
        } catch (std::exception& e){
            std::cout << e.what() << std::endl;
        }
    }

#else

    PangenomeMAT::Tree *T;

    if(globalVm.count("help")){
        std::cout << globalDesc;
        return;
    } else if(globalVm.count("input-file")){
        std::string fileName = globalVm["input-file"].as< std::string >();
        std::ifstream inputStream(fileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new PangenomeMAT::Tree(inputStream);

        auto treeBuiltEnd= std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        inputStream.close();

    } else {
        printError("Incorrect Format");
        // std::cout << "\033[1;31m" << "Error: " << "\033[0m" << "Incorrect Format\n";
        std::cout << globalDesc;
        return;
    }

    char** splitCommandArray;

    while(true){
        std::cout << "> ";

        std::string command;
        std::getline (std::cin, command);
        stripString(command);
        std::vector< std::string > splitCommand;
        PangenomeMAT::stringSplit(command, ' ', splitCommand);

        splitCommandArray = new char*[splitCommand.size()];
        for(size_t i = 0; i < splitCommand.size(); i++){
            splitCommandArray[i] = new char[splitCommand[i].length() + 1];
            strcpy(splitCommandArray[i], splitCommand[i].c_str());
        }

        try{
            if(strcmp(splitCommandArray[0], "summary") == 0){
                // If command was summary
                auto summaryStart = std::chrono::high_resolution_clock::now();
                T->printSummary();
                auto summaryEnd = std::chrono::high_resolution_clock::now();

                std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

                std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
            } else if(strcmp(splitCommandArray[0], "fasta") == 0){
                // If command was fasta
                po::variables_map fastaVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(fastaDesc).positional(fastaPositionArgumentDesc).run(), fastaVm);

                if(fastaVm.count("help")){
                    std::cout << fastaDesc;
                } else {
                    po::notify(fastaVm);

                    std::string fileName = fastaVm["output-file"].as< std::string >();

                    bool aligned = false;
                    size_t num_cores = 0;
                    if(fastaVm.count("aligned")){
                        aligned = true;
                    }
                    if(fastaVm.count("num_cores")){
                        num_cores = fastaVm["num_cores"].as< size_t >();
                    }
                    
                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto fastaStart = std::chrono::high_resolution_clock::now();
                    
                    T->printFASTA(fout, aligned, num_cores);

                    auto fastaEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

                    fout.close();

                }
            } else if(strcmp(splitCommandArray[0], "write") == 0){
                // If command was write
                po::variables_map writeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(writeDesc).positional(writePositionArgumentDesc).run(), writeVm);

                if(writeVm.count("help")){
                    std::cout << writeDesc;
                } else {
                    po::notify(writeVm);
                    std::string fileName = writeVm["output-file"].as< std::string >();
                    std::filesystem::create_directory("./pmat");
                    std::ofstream fout("./pmat/" + fileName + ".pmat");

                    auto writeStart = std::chrono::high_resolution_clock::now();
                    
                    T->writeToFile(fout);

                    auto writeEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds writeTime = writeEnd - writeStart;

                    std::cout << "\nTree Write execution time: " << writeTime.count() << " nanoseconds\n";

                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "subtree") == 0){
                // If command was subtree
                po::variables_map subtreeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(subtreeDesc).positional(subtreePositionArgumentDesc).run(), subtreeVm);

                if(subtreeVm.count("help")){
                    std::cout << subtreeDesc;
                } else {
                    po::notify(subtreeVm);
                    std::string outputFileName = subtreeVm["output-file"].as< std::string >();
                    std::vector< std::string > nodeIds;
                    std::string nodeId;

                    if(subtreeVm.count("input-file")){
                        std::string inputFileName = subtreeVm["input-file"].as< std::string >();
                        std::ifstream fin(inputFileName);

                        while(fin >> nodeId){
                            nodeIds.push_back(nodeId);
                        }

                        fin.close();
                    } else if(subtreeVm.count("node-ids")){
                        nodeIds = subtreeVm["node-ids"].as< std::vector< std::string > >();
                    } else {
                        printError("No source of node ids provided");
                        std::cout << subtreeDesc;
                    }

                    if(subtreeVm["newick"].as< bool >()){
                        std::filesystem::create_directory("./pmat");
                        std::ofstream fout("./newick/" + outputFileName + ".newick");

                        auto subtreeStart = std::chrono::high_resolution_clock::now();

                        fout << T->getNewickString(T->subtreeExtractParallel(nodeIds));

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";
                        fout.close();
                    } else {
                        std::filesystem::create_directory("./pmat");
                        std::ofstream fout("./pmat/" + outputFileName + ".pmat");

                        auto subtreeStart = std::chrono::high_resolution_clock::now();

                        T->writeToFile(fout, T->subtreeExtractParallel(nodeIds));

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";
                        fout.close();
                    }

                }

            } else if(strcmp(splitCommandArray[0], "vcf") == 0){
                // If command was vcf
                po::variables_map vcfVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(vcfDesc).positional(vcfPositionArgumentDesc).run(), vcfVm);

                if(vcfVm.count("help")){
                    std::cout << vcfDesc;
                } else {
                    po::notify(vcfVm);

                    std::string reference = vcfVm["reference"].as< std::string >();
                    std::string fileName = vcfVm["output-file"].as< std::string >();

                    std::filesystem::create_directory("./vcf");
                    std::ofstream fout("./vcf/" + fileName + ".vc");

                    auto vcfStart = std::chrono::high_resolution_clock::now();

                    T->printVCFParallel(reference, fout);

                    auto vcfEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;

                    std::cout << "\nVCF execution time: " << vcfTime.count() << " nanoseconds\n";

                    fout.close();
                }

            } else if(strcmp(splitCommandArray[0], "newick") == 0){
                std::cout << T->getNewickString(T->root) << std::endl;
            } else if(strcmp(splitCommandArray[0], "annotate") == 0){
                // If command was annotate
                po::variables_map annotateVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(annotateDesc).positional(annotatePositionArgumentDesc).run(), annotateVm);

                if(annotateVm.count("help")){
                    std::cout << annotateDesc;
                } else {
                    po::notify(annotateVm);
                    std::string fileName = annotateVm["input-file"].as< std::string >();
                    std::ifstream fin(fileName);

                    auto annotateStart = std::chrono::high_resolution_clock::now();

                    T->annotate(fin);

                    auto annotateEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds annotateTime = annotateEnd - annotateStart;

                    std::cout << "Annotate time: " << annotateTime.count() << " nanoseconds\n";
                }

            } else if(strcmp(splitCommandArray[0], "search") == 0){
                po::variables_map searchVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(searchDesc).positional(searchPositionArgumentDesc).run(), searchVm);

                std::vector< std::string > annotationVector = searchVm["keywords"].as< std::vector< std::string > >();

                for(auto word: annotationVector){
                    stripString(word);
                    std::cout << word << ": ";
                    auto result = T->searchByAnnotation(word);
                    for(auto r: result){
                        std::cout << r << ";";
                    }
                    std::cout << std::endl;
                }

            } else if(strcmp(splitCommandArray[0], "exit") == 0){
                return;
            }

        } catch (std::exception& e){
            std::cout << e.what() << std::endl;
        }

        for(size_t i = 0; i < splitCommand.size(); i++){
            delete [] splitCommandArray[i];
        }
        delete [] splitCommandArray;

        std::cout << std::endl;
    }

    delete T;

#endif

}

int main(int argc, char* argv[]){

#ifdef NEW_PARSER
    updatedParser(argc, argv);
#else

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

            std::vector< std::string > splitCommand;

            PangenomeMAT::stringSplit(command, ' ', splitCommand);

            if(splitCommand.size() == 1 && splitCommand[0] == "summary"){
                auto summaryStart = std::chrono::high_resolution_clock::now();
                T.printSummary();
                auto summaryEnd = std::chrono::high_resolution_clock::now();

                std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

                std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
            } else if(splitCommand.size() >= 2 && splitCommand[0] == "fasta"){
                if(splitCommand.size() > 2 && splitCommand[1].substr(0,11) == "--parallel="){
                    int parallelism = std::stoi(splitCommand[1].substr(11));
                    std::string fileName = splitCommand[2];
                    bool aligned = false;
                    if(splitCommand.size() == 4 && splitCommand[2] == "--aligned"){
                        aligned = true;
                        fileName = splitCommand[3];
                    }

                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto fastaStart = std::chrono::high_resolution_clock::now();
                    
                    T.printFASTA(fout, aligned, parallelism);

                    auto fastaEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

                    fout.close();

                } else {
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

                    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

                    fout.close();
                }
            } else if(splitCommand.size() == 2 && splitCommand[0] == "write") {
                std::string fileName = splitCommand[1];
                std::filesystem::create_directory("./pmat");
                std::ofstream fout("./pmat/" + fileName + ".pmat");

                auto writeStart = std::chrono::high_resolution_clock::now();
                
                T.writeToFile(fout);

                auto writeEnd = std::chrono::high_resolution_clock::now();
                
                std::chrono::nanoseconds writeTime = writeEnd - writeStart;

                std::cout << "\nTree Write execution time: " << writeTime.count() << " nanoseconds\n";

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

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";

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

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";

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

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";

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

                    std::cout << "\nParallel Subtree Extract execution time: " << subtreeTime.count() << " nanoseconds\n";

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

                T.printVCFParallel(reference, fout);

                auto vcfEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;

                std::cout << "\nVCF execution time: " << vcfTime.count() << '\n';

                fout.close();

                // Debugging

                // std::ifstream fin("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_012892.2", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_007779.1", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_000913.3", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_002695.2", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_004431.1", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_013654.1", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_013364.1", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_011415.1", fin);
                // fin.close();

                // fin.open("./vcf/" + fileName + ".vc");
                // T.getSequenceFromVCF("NC_013353.1", fin);
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

                // std::vector< std::string > annotationVector = splitString(restOfCommand, ',');
                std::vector< std::string > annotationVector;
                
                PangenomeMAT::stringSplit(restOfCommand, ',', annotationVector);

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
#endif

}