#define MAT_V2

#include <iostream>
#include <chrono>
#include <filesystem>
#include <boost/program_options.hpp>
#include <json/json.h>

#include <fstream>

#include "PangenomeMAT.hpp"

#include "spoa/spoa.hpp"

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
po::options_description mafDesc("MAF Writer Command Line Arguments");
po::positional_options_description mafPositionArgumentDesc;
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
po::options_description generateGFADesc("Generate GFA Command Line Arguments");
po::positional_options_description generateGFAArgumentDesc;
po::options_description GFAToFASTADesc("GFA to Fasta writer Command Line Arguments");
po::positional_options_description GFAToFASTAArgumentDesc;
po::options_description groupWriteDesc("Group MAT Writer Command Line Arguments");
po::positional_options_description groupWritePositionArgumentDesc;
po::options_description sequenceExtractDesc("Sequence Extract Command Line Arguments");
po::positional_options_description sequenceExtractPositionArgumentDesc;
po::options_description groupFastaDesc("Tree Group FASTA writer Command Line Arguments");
po::positional_options_description groupFastaPositionArgumentDesc;

void setupOptionDescriptions(){
    // Global option descriptions
    globalDesc.add_options()
        ("help", "produce help message")
        ("input-file,I", po::value< std::string >(), "PanMAT input file path")
        ("gfa-in", po::value< std::string >(), "create PanMAT from GFA file")
        ("pangraph-in", po::value< std::string >(), "create PanMAT from Pangraph file")
        ("msa-in", po::value< std::string >(), "create PanMAT from MSA file")
        ("optimize", "currently UNSUPPORTED: whether given msa file should be optimized or not")
        ("newick-in", po::value< std::string >(), "Input file path for file containing newick string")
        ("tree-group", po::value< std::vector< std::string > >()->multitoken(), "File paths of PMATs to generate tree group")
        ("mutation-file", po::value< std::string >(), "File path of complex mutation file for tree group")
        ("panman-in", po::value< std::string >(), "Input file path for PanMAT Group")
    ;

    // Adding input file as positional argument
    globalPositionArgumentDesc.add("input-file", -1);

    // FASTA option descriptions
    fastaDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
        ("aligned", "print in aligned format")
        ("parallel", "Whether we should execute in parallel or not")
    ;

    // Adding output file as positional argument
    fastaPositionArgumentDesc.add("output-file", -1);

    // MAF option descriptions
    mafDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    mafPositionArgumentDesc.add("output-file", -1);

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
        ("node-ids", po::value< std::vector< std::string > >()->multitoken()->required(), "Node IDs to extract")
    ;

    // Adding output file as positional argument
    subtreePositionArgumentDesc.add("output-file", -1);

    // Sequence Extract option descriptions
    sequenceExtractDesc.add_options()
        ("help", "produce help message")
        ("sequence", po::value< std::string >()->required(), "Sequence name")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    sequenceExtractPositionArgumentDesc.add("output-file", -1);

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
        ("input-file", po::value< std::string >()->required(), "Name of the file containing annotation info")
    ;

    // Adding input file as positional argument
    annotatePositionArgumentDesc.add("input-file", -1);

    // Search by annotation option descriptions
    searchDesc.add_options()
        ("help", "produce help message")
        ("keywords", po::value< std::vector< std::string > >()->multitoken(), "list of keywords to search for")
    ;

    searchPositionArgumentDesc.add("keywords", -1);
    
    // Generate GFA option descriptions
    generateGFADesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    generateGFAArgumentDesc.add("output-file", -1);

    // GFA to FASTA option descriptions
    GFAToFASTADesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
        ("input-file", po::value< std::string >()->required(), "Input file name")
    ;

    GFAToFASTAArgumentDesc.add("output-file", -1);


    // Tree Group FASTA option descriptions
    groupFastaDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    groupFastaPositionArgumentDesc.add("output-file", -1);

    // Group MAT Writer option descriptions
    groupWriteDesc.add_options()
        ("help", "produce help message")
        ("output-file", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    groupWritePositionArgumentDesc.add("output-file", -1);

}

void printError(std::string e){
    std::cout << "\033[1;31m" << "Error: " << "\033[0m" << e << "\n";
}

void updatedParser(int argc, char* argv[]){

    // Setup boost::program_options
    setupOptionDescriptions();

    // Initial command line arguments consisting of input file types
    po::variables_map globalVm;
    po::store(po::command_line_parser(argc, argv).options(globalDesc).positional(globalPositionArgumentDesc).run(), globalVm);
    po::notify(globalVm);

    PangenomeMAT::Tree *T = nullptr;
    PangenomeMAT::TreeGroup *TG = nullptr;

    if(globalVm.count("help")){
        std::cout << globalDesc;
        return;
    } else if(globalVm.count("input-file")){
        std::string fileName = globalVm["input-file"].as< std::string >();
        std::ifstream inputStream(fileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new PangenomeMAT::Tree(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        inputStream.close();

        // std::ifstream fin("maf/ecoli_10_test.maf");
        // std::ofstream fout("maf/mafTest.fasta");
        // T->generateSequencesFromMAF(fin, fout);
        // fout.close();
        // fin.close();
        
    } else if(globalVm.count("panman-in")){
        std::string fileName = globalVm["panman-in"].as< std::string >();
        std::ifstream inputStream(fileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        TG = new PangenomeMAT::TreeGroup(inputStream);

        // for(auto& tree: TG->trees){
        //     std::cout << tree.root->identifier << std::endl;
        // }

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        inputStream.close();  
    } else if(globalVm.count("gfa-in")){
        std::string fileName = globalVm["gfa-in"].as< std::string >();
        if(!globalVm.count("newick-in")){
            printError("File containing newick string not provided!");
            return;
        }
        std::string newickFileName = globalVm["newick-in"].as< std::string >();

        std::cout << "Creating PanMAT from GFA" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new PangenomeMAT::Tree(inputStream, newickInputStream, PangenomeMAT::FILE_TYPE::GFA);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

    } else if(globalVm.count("pangraph-in")){
        std::string fileName = globalVm["pangraph-in"].as< std::string >();
        if(!globalVm.count("newick-in")){
            printError("File containing newick string not provided!");
            return;
        }
        std::string newickFileName = globalVm["newick-in"].as< std::string >();

        std::cout << "Creating PanMAT from Pangraph" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new PangenomeMAT::Tree(inputStream, newickInputStream, PangenomeMAT::FILE_TYPE::PANGRAPH);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

    } else if(globalVm.count("msa-in")){
        std::string fileName = globalVm["msa-in"].as< std::string >();
        if(!globalVm.count("newick-in")){
            printError("File containing newick string not provided!");
            return;
        }
        bool optimize = false;
        if(globalVm.count("optimize")){
            optimize = true;
        }

        std::string newickFileName = globalVm["newick-in"].as< std::string >();

        std::cout << "Creating PanMAT from MSA" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        if(!optimize){
            T = new PangenomeMAT::Tree(inputStream, newickInputStream, PangenomeMAT::FILE_TYPE::MSA);
        } else {
            T = new PangenomeMAT::Tree(inputStream, newickInputStream, PangenomeMAT::FILE_TYPE::MSA_OPTIMIZE);
        }

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

    } else if(globalVm.count("tree-group")){
        std::vector< std::string > fileNames;

        std::string mutationFileName;
        if(!globalVm.count("mutation-file")){
            printError("File containing complex mutations not provided!");
            return;
        }

        fileNames = globalVm["tree-group"].as< std::vector< std::string > >();
        mutationFileName = globalVm["mutation-file"].as< std::string >();
        
        std::ifstream mutationFile(mutationFileName);

        std::vector< std::ifstream > files;
        for(auto u: fileNames){
            files.emplace_back(u);
        }

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        TG = new PangenomeMAT::TreeGroup(files, mutationFile);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        mutationFile.close();
        for(auto& u: files){
            u.close();
        }

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
                    bool parallel = 0;
                    if(fastaVm.count("aligned")){
                        aligned = true;
                    }
                    if(fastaVm.count("parallel")){
                        parallel = true;
                    }
                    
                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto fastaStart = std::chrono::high_resolution_clock::now();
                    
                    if(parallel){
                        T->printFASTAParallel(fout, aligned);
                    } else {
                        T->printFASTA(fout, aligned);
                    }

                    auto fastaEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

                    fout.close();

                }
            } else if(strcmp(splitCommandArray[0], "maf") == 0){
                // If command was maf
                po::variables_map mafVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(mafDesc).positional(mafPositionArgumentDesc).run(), mafVm);

                if(mafVm.count("help")){
                    std::cout << mafDesc;
                } else {
                    po::notify(mafVm);

                    std::string fileName = mafVm["output-file"].as< std::string >();

                    std::filesystem::create_directory("./maf");
                    std::ofstream fout("./maf/" + fileName + ".maf");

                    auto mafStart = std::chrono::high_resolution_clock::now();
                    
                    T->printMAF(fout);

                    auto mafEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds mafTime = mafEnd - mafStart;

                    std::cout << "\nMAF execution time: " << mafTime.count() << " nanoseconds\n";

                    fout.close();

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

            } else if(strcmp(splitCommandArray[0], "sequence") == 0){
                // If command was maf
                po::variables_map sequenceExtractVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(sequenceExtractDesc).positional(sequenceExtractPositionArgumentDesc).run(), sequenceExtractVm);

                if(sequenceExtractVm.count("help")){
                    std::cout << sequenceExtractDesc;
                } else {
                    po::notify(sequenceExtractVm);

                    std::string fileName = sequenceExtractVm["output-file"].as< std::string >();
                    std::string sequenceName = sequenceExtractVm["sequence"].as< std::string >();

                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto sequenceExtractStart = std::chrono::high_resolution_clock::now();
                    
                    std::string sequenceString = T->getStringFromReference(sequenceName, false);
                    fout << ">" << sequenceName << '\n';
                    for(size_t i = 0; i < sequenceString.length(); i+=70){
                        fout << sequenceString.substr(i,std::min(sequenceString.length()-i, (size_t)70)) << '\n';
                    }

                    auto sequenceExtractEnd = std::chrono::high_resolution_clock::now();
                    
                    std::chrono::nanoseconds sequenceExtractTime = sequenceExtractEnd - sequenceExtractStart;

                    std::cout << "\nSequence Extract execution time: " << sequenceExtractTime.count() << " nanoseconds\n";

                    fout.close();

                }
            } else if(strcmp(splitCommandArray[0], "newick") == 0){
                if(T){
                    std::cout << T->getNewickString(T->root) << std::endl;
                } else if(TG){
                    std::cout << "Printing newick string of each PanMAT in PanMAN" << std::endl;
                    int index = 0;
                    for(auto& t: TG->trees){
                        std::cout << index++ << ": " << t.getNewickString(t.root) << std::endl;
                    }
                }
            } else if(strcmp(splitCommandArray[0], "genGFA") == 0){
                // If command was genGFA
                po::variables_map generateGFAVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(generateGFADesc).positional(generateGFAArgumentDesc).run(), generateGFAVM);
                if(generateGFAVM.count("help")){
                    std::cout << generateGFADesc;
                } else {
                    po::notify(generateGFAVM);
                    std::string fileName = generateGFAVM["output-file"].as< std::string >();
                    std::filesystem::create_directory("./gfa");
                    std::ofstream fout("./gfa/"+fileName+".gfa");

                    auto generateVGStart = std::chrono::high_resolution_clock::now();

                    T->convertToGFA(fout);

                    auto generateVGEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds generateVGTime = generateVGEnd - generateVGStart;

                    std::cout << "GFA generation time: " << generateVGTime.count() << " nanoseconds\n";
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "gfa-fasta") == 0){

                // If command was gfa-fasta

                po::variables_map GFAToFASTAVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(GFAToFASTADesc).positional(GFAToFASTAArgumentDesc).run(), GFAToFASTAVM);
                if(GFAToFASTAVM.count("help")){
                    std::cout << GFAToFASTADesc;
                } else {
                    po::notify(GFAToFASTAVM);
                    std::string inputFileName = GFAToFASTAVM["input-file"].as< std::string >();
                    std::string outputFileName = GFAToFASTAVM["output-file"].as< std::string >();
                    std::filesystem::create_directory("./fasta");

                    std::ifstream fin(inputFileName);
                    std::ofstream fout("./fasta/"+outputFileName+".fasta");

                    auto generateVGStart = std::chrono::high_resolution_clock::now();

                    T->printFASTAFromGFA(fin, fout);

                    auto generateVGEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds generateVGTime = generateVGEnd - generateVGStart;

                    std::cout << "FASTA generation time: " << generateVGTime.count() << " nanoseconds\n";

                    fin.close();
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "groupFasta") == 0){
                // If FASTA for tree group is required

                po::variables_map groupFastaVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(groupFastaDesc).positional(groupFastaPositionArgumentDesc).run(), groupFastaVM);
                if(groupFastaVM.count("help")){
                    std::cout << groupFastaDesc;
                }
                std::filesystem::create_directory("./fasta");

                std::string outputFileName = groupFastaVM["output-file"].as< std::string >();
                std::ofstream fout("./fasta/" + outputFileName + ".fasta");

                auto groupFastaStart = std::chrono::high_resolution_clock::now();

                TG->printFASTA(fout);

                auto groupFastaEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds groupFASTATime = groupFastaEnd - groupFastaStart;

                std::cout << "Group FASTA write time: " << groupFASTATime.count() << " nanoseconds\n";

                fout.close();

            } else if(strcmp(splitCommandArray[0], "groupWrite") == 0){
                // If FASTA for tree group is required

                po::variables_map groupWriteVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray).options(groupWriteDesc).positional(groupWritePositionArgumentDesc).run(), groupWriteVM);
                if(groupWriteVM.count("help")){
                    std::cout << groupWriteDesc;
                }
                std::filesystem::create_directory("./pman");

                std::string outputFileName = groupWriteVM["output-file"].as< std::string >();
                std::ofstream fout("./pman/" + outputFileName + ".pman");

                auto groupWriteStart = std::chrono::high_resolution_clock::now();

                TG->writeToFile(fout);

                auto groupWriteEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds groupWriteTime = groupWriteEnd - groupWriteStart;

                std::cout << "Group PanMAT write time: " << groupWriteTime.count() << " nanoseconds\n";

                fout.close();

            } else if(strcmp(splitCommandArray[0], "cplx-mutations") == 0){
                // Print Complex Mutations
                std::cout << "Complex Mutations:" << std::endl;
                TG->printComplexMutations();

            } else if(strcmp(splitCommandArray[0], "exit") == 0){
                return;
            }
        } catch (std::exception& e){
            std::cout << e.what() << std::endl;
        }
    }

}

void debuggingCode(){
    std::ifstream fin("../../ecoli_50.json");
    Json::Value pangraphData;

    fin >> pangraphData;

    std::cout << "Blocks orig: " << pangraphData["blocks"].size() << std::endl;

    std::map< std::string, std::map< std::pair< std::string, int >, bool > > paths;

    // load blocks
    // int s = 0;
    // std::map< std::string, int > blockLens;

    for(size_t i = 0; i < pangraphData["blocks"].size(); i++){
        std::string blockId = pangraphData["blocks"][(int)i]["id"].asString();
        if(blockId == "MBEUQKOMBV"){
            std::cout << pangraphData["blocks"][(int)i]["sequence"].asString() << std::endl;
        } else {
            continue;
        }
        // blockLens[blockId] = pangraphData["blocks"][(int)i]["sequence"].asString().length();
        // if(blockId == "YTRMIRRHDS"){
        //     for(int j = 0; j < pangraphData["blocks"][(int)i]["positions"].size(); j++){
        //         if(pangraphData["blocks"][(int)i]["positions"][j][0]["name"].asString() == "KX894803.1"){
        //             std::cout << pangraphData["blocks"][(int)i]["positions"][j][1][0].asInt() << " " << pangraphData["blocks"][(int)i]["positions"][j][1][1].asInt() << std::endl;
        //         }
        //     }
        // }

        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["mutate"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["number"].asInt();
            bool strand = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["strand"].asBool();
            if(seqName == "NZ_CP006027.1" && blockId == "MBEUQKOMBV"){
                for(size_t k = 0; k < pangraphData["blocks"][(int)i]["mutate"][(int)j][1].size(); k++){
                    std::string mutationString = pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString();
                    std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                    std::cout << "S " << number << " " << strand << " " << pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt() << " " << mutationString << std::endl;
                }
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["insert"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["number"].asInt();
            bool strand = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["strand"].asBool();
            if(seqName == "NZ_CP006027.1" && blockId == "MBEUQKOMBV"){
                for(size_t k = 0; k < pangraphData["blocks"][(int)i]["insert"][(int)j][1].size(); k++){
                    std::string mutationString = pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][1].asString();
                    std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                    std::cout << "I " << number << " " << strand << " " << pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][0].asInt() << " " << pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][1].asInt() << " " << " " << mutationString << std::endl;
                    // s+=mutationString.length();
                }
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["delete"].size(); j++){
            std::string seqName = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["name"].asString();
            size_t number = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["number"].asInt();
            bool strand = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["strand"].asBool();
            if(seqName == "NZ_CP006027.1" && blockId == "MBEUQKOMBV"){
                for(size_t k = 0; k < pangraphData["blocks"][(int)i]["delete"][(int)j][1].size(); k++){
                    std::cout << "D " << number << " " << strand << " " << pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][0].asInt() << " " << pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][1].asInt() << std::endl;
                    // s-=pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][1].asInt();
                }
            }
        }
    }

    for(size_t i = 0; i < pangraphData["paths"].size(); i++){
        Json::Value path = pangraphData["paths"][(int)i];
        // std::map< std::pair< std::string, bool >, int > blocks;
        // int ctr = 0;

        if(path["name"].asString() == "NZ_CP006027.1"){
            std::cout << "OFFSET: " << path["offset"].asInt() << std::endl;
        }

        for(size_t j = 0; j < path["blocks"].size(); j++){
            if(path["name"].asString() == "NZ_CP006027.1" && path["blocks"][(int)j]["id"].asString() == "MBEUQKOMBV"){
                std::cout << "Strand: " << path["blocks"][(int)j]["strand"].asBool() << std::endl;
            }
            if(path["name"].asString() == "NZ_CP006027.1" && j>0 && path["position"][(int)j].asInt() < path["position"][(int)j-1].asInt()){
                std::cout << path["blocks"][(int)j-1]["id"].asString() << " " << path["blocks"][(int)j-1]["strand"].asBool() << " " << path["blocks"][(int)j]["id"].asString() << " " << path["blocks"][(int)j]["strand"].asBool() << " " << path["position"][(int)j-1].asInt() << " " << path["position"][(int)j].asInt() << std::endl;
            }
            // if(path["name"].asString() == "NZ_CP013985.1" && path["blocks"][(int)j]["id"].asString() == "STETJDHNZS"){
            //     std::cout << path["blocks"][(int)j]["strand"].asBool() << std::endl;
            // }
            // s+=blockLens[path["blocks"][(int)j]["id"].asString()];
            // paths[path["name"].asString()][std::make_pair(path["blocks"][(int)j]["id"].asString(), path["blocks"][(int)j]["number"].asInt())] = path["blocks"][(int)j]["strand"].asBool();
            // if(path["name"].asString() == "KX894803.1"){
            //     std::cout << ctr << " " << path["blocks"][(int)j]["id"].asString() << std::endl;
            //     ctr++;
            // }
            // if(blocks.find(std::make_pair(path["blocks"][(int)j]["id"].asString(), !path["blocks"][(int)j]["strand"].asBool())) != blocks.end()){
            //     std::cout << path["name"].asString() << " " << path["blocks"][(int)j]["id"].asString() << std::endl;
            //     // std::cout << blocks[std::make_pair(path["blocks"][(int)j]["id"].asString(), !path["blocks"][(int)j]["strand"].asBool())] << " " << path["blocks"][(int)j]["number"].asInt() << std::endl;
            // }
            // blocks[std::make_pair(path["blocks"][(int)j]["id"].asString(), path["blocks"][(int)j]["strand"].asBool())] = path["blocks"][(int)j]["number"].asInt();
        }
    }
    // std::cout << s << std::endl;

    // // load blocks
    // for(size_t i = 0; i < pangraphData["blocks"].size(); i++){
    //     std::string blockId = pangraphData["blocks"][(int)i]["id"].asString();
    //     std::string consensusSeq = pangraphData["blocks"][(int)i]["sequence"].asString();
    //     if(consensusSeq == "CAGGCCCCACGGGGATACCGGGTGCTTACCAAAGTTACGACCTTCACCACCACCATGTGGGTGGTCTACTGGGTTCATCGCAGTACCGCGAACGGTAGGACGAACACCACGCCAGCGTGCAGCACCTGCTTTACCCAGAACGCGCAGCATATGCTCAGCATTGCCAACTTCGCCCAGAGTCGCGCGGCAGTCTGCTTCGACTTTACGCATTTCACCAGAACGCAGACGCAGGGTGACATAAGCACCGTCACGAGCAACGATCTGAACGTAAGTACCAGCGGAACGAGCCAGCTGACCGCCTTTACCTGGTTTCATTTCTACGTTATGAACGGTAGAACCAACCGGGATATTGCGCATCGGCAGGGTGTTGCCTGCTTTGATTGCAGCATCAACGCCAGACTGAATCTGGTCGCCAGCTTTCAGGCCTTTCGGGGCCAGGATGTAACGGCGCTCGCCGTCTTTGTACAGAACCAGCGCGATGTTCGCGGAGCGGTTCGGATCGTACTCAAGACGTTCAACAACTGCCGGGATACCATCTTTGTTGCGTTTGAAGTCAACAATACGGTAAGCCTGCTTGTGGCCACCACCGATATGACGGGTGGTGATACGGCCATTGTTGTTACGACCACCGGATTTGCTGTTTTTTTCCAGCAACGGAGCAAAAGGTTTGCCCTTGTGCAGCTCTGGGTTAACCACTTTAACGACGTGGCGACGACCCGGAGATGTCGGTTTACATTTAACAACTGCCATTGTATTACTCCTCCGACTTACTCAGCGCCGCCAACGAAGTCCAGATTCTGGCCTTCTTTCAGGGTGACGTAAGCTTTTTTCCAGTCGCTACGACGACCGATACGCTGTCCGTGACGTTTAACTTTCCCTTTAACAACCAGGGTGTTAACGACTTCGACTTCGACTTCAAACAGTTTCTGCACAGCAGCTTTGATTTCTGCTTTGGTCGCGTCTTTAGCAACTTTGAGAACGATGGTGTTTGTTTTTTCCATCGCAGTAGACGCTTTTTCAGAAACGTGCGGTGCGCGCAGCACCTTCAGCAGACGTTCTTCACGAATCATGCCAGCATCTCCTCAACTTGCTTAACAGCATCAGCAGTCATTACGACTTTGTCGAAGGCGATCAGGCTAACCGGGTCGATACCGTTCGCATCACGTACGTCAACCTTATGCAGGTTACGTGCGGCCAGGAACAGGTTCTCGTCCAGCTCACCGGTGATGATCAGCACATCTTCCAGAGCCATGTCTTTCAGTTTCTGTGCCAGCAGCTTAGTTTTCGGCGCTTCAACAGAGAACTTCTCGACAACGATCAGACGATCCTGACGTACCAGTTCGGACAGAATGCTTTTCAGCGCGCCGCGGTACATCTTCTTGTTAACTTTTTGACTGTGGTCCTGCGGACGAGCAGCGAAGGTCACGCCACCGGAACGCCAGATCGGGCTCTTGATAGAACCAGAACGCGCACGGCCGGTACCTTTCTGGCGCCACGGTTTTTTACCGGAGCCAGTGATTTCAGCACGAGTCTTCTGAGCACGAGTGCCCTGACGCGCACCGGCTGCATAAGCAACAACAACCTGGTGAACCAGCGCTTCGTTGAAATCACGACCGAAGGTAGTTTCGGAAACAGTCAGCGCGCTCTGCGCGTCTTTCAATACTAATTCCATTGCTATCTCCTCACGCCTTCACAGCTGGTTTAACGATCAGGTCGCTACCGGTCGCACCCGGGACAGCACCTTTAACCAGCAGCAGGTTGCGCTCAGCGTCAACACGTACTACGTCCAGGCTCTGAACGGTTACACGCTCATTACCCAGCTGACCTGCCATTTTCTTGCCTTTGAACACTTTGCCCGGAGTCTGGTTCTGACCGATAGAACCCGGAACGCGGTGGGACAAGGAGTTACCGTGAGTAGCGTCCTGGGTACGGAAGTTCCAGCGCTTAACGGTACCAGCAAAACCTTTACCTTTAGAGGTACCGGTTACGTCAACTTTTTTAACGTCAGCAAACAGTTCAACGCTAATGTTCTGACCTACGGTGAATTCTTCGCCGTCAGCCAGACGGAATTCCCACAGACCACGGCCAGCTTCAACGCCAGCTTTAGCGAAGTGACCCGCTTCCGGCTTGGTTACACGGTTAGCTTTTTTAGCACCGGTGGTAACCTGGATTGCGCGGTAGCCATCGTTAGCCAGGTCTTTAACCTGAGTAACACGGTTTGCTTCAACTTCGATTACGGTTACTGGGATAGAAACGCCATCTTCAGTGAAGATGCGGGTCATGCCCACTTTTTTACCGACTAAACCAATCATTGTTTCAACCTCTCAATCGCTCGATGACCTGATTAACCCAGGCTGATCTGCACGTCTACACCGGCAGCCAGATCCAGACGCATCAGAGCATCAACGGTTTTTTCAGTTGGCTCAACGATGTCAACCAGACGCTTGTGAGTGCGAATTTCGTACTGGTCACGCGCGTCTTTGTTGACGTGCGGGGAGATCAGAACGGTAAAGCGCTCTTTGCGGGTCGGCAGCGGGATCGGACCACGGACTTGCGCACCAGTGCGCTTAGCAGTCTCGACGATTTCCGCGGTTGATTGATCGATCAGACGATGATCAAACGCTTTCAGGCGGATACGGATTCTTTGGTTCTGCATGAGACCAGAGCTCCAATTATTTTATAGACGAAATGATTACTCCTCATACCCATTACGATTGATGGGAGAGTGTAACCGTTCTTACATAGTCCCCCGATTGGGGACATTGTTTGATAACCCAAATGCGATTACCAGGGTTCATATTGAACCGGCCGTCAACTTTGACAGGCCCGCGCATTATACGCAAATTCTGGCGCGACGCAAGTGCTGTTTAGAAAATCATCACCAGGATTTGCGC"){
    //         std::cout << blockId << std::endl;
    //     }
    //     // if(i == 3957){
    //     //     std::cout << blockId << " " << consensusSeq.substr(2774-5,5) << std::endl;
    //     // }
    //     // if(consensusSeq == "TGCTTTTT"){
    //     //     std::cout << "FOUND: " << std::endl;
    //     // }
    //     // if(consensusSeq.length() == 8){
    //     //     std::cout << blockId << " " << consensusSeq << std::endl;
    //     // }

    //     // for(size_t j = 0; j < pangraphData["blocks"][(int)i]["mutate"].size(); j++){
    //     //     std::string seqName = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["name"].asString();
    //     //     size_t number = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["number"].asInt();
    //     //     bool strand = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["strand"].asBool();
    //     //     if(paths[seqName][std::make_pair(blockId, number)] != strand){
    //     //         std::cout << "Problem!!!" << std::endl;
    //     //     }
    //     // }
    //     // for(size_t j = 0; j < pangraphData["blocks"][(int)i]["insert"].size(); j++){
    //     //     std::string seqName = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["name"].asString();
    //     //     size_t number = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["number"].asInt();
    //     //     bool strand = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["strand"].asBool();
    //     //     if(paths[seqName][std::make_pair(blockId, number)] != strand){
    //     //         std::cout << "Problem!!!" << std::endl;
    //     //     }
    //     // }
    //     // for(size_t j = 0; j < pangraphData["blocks"][(int)i]["delete"].size(); j++){
    //     //     std::string seqName = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["name"].asString();
    //     //     size_t number = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["number"].asInt();
    //     //     bool strand = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["strand"].asBool();
    //     //     if(paths[seqName][std::make_pair(blockId, number)] != strand){
    //     //         std::cout << "Problem!!!" << std::endl;
    //     //     }
    //     // }
    // }
}

int main(int argc, char* argv[]){

    // debuggingCode();
    tbb::task_scheduler_init init(32);

    updatedParser(argc, argv);

}