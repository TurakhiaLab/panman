#include <iostream>
#include <chrono>
#include <filesystem>
#include <tbb/parallel_for_each.h>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
// #include <boost/iostreams/filter/xz.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <json/json.h>

#include <fstream>

#include "panmanUtils.hpp"

namespace po = boost::program_options;

// Remove spaces from beginning and end of given string
void stripStringInPlace(std::string& s) {
    while(s.length() && s[s.length() - 1] == ' ') {
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++) {
        if(s[i] != ' ') {
            s = s.substr(i);
            return;
        }
    }
}

std::string printNucMut(int32_t mutInfo){
    std::string s = "Type ";
    s += std::to_string(mutInfo&0xf);
    s += " Length ";
    s += std::to_string(mutInfo>>4);
    return s;
}

std::string printNucs(int32_t mutInfo, int32_t nucs){
    std::string s = " Chars: ";
    int len = mutInfo >> 4;
    for (int i=0; i<len; i++) {
        s += std::to_string(nucs&0x4);
        s += " ";
        nucs = nucs >> 4;
    }
    return s;
}

void checkFunction(panmanUtils::Tree *T) {
    std::cout << T->root->identifier << std::endl;
    std::ofstream o("new.fa");
    T->printFASTA(o);

    // get node id
    panmanUtils::Node * node = T->allNodes["England/MILK-338D3D9/2022|OV784995.1|2022-01-21"];
    std::cout << "Found Node: " << node->identifier << std::endl;

    while (node != nullptr) {
        std::cout << "node: " << node->identifier << std::endl;
        // Pring block mutations
        std::cout << "Block mutations" << std::endl;
        for(auto &u: node->blockMutation) {
            std::cout << "\t" << u.blockMutInfo << " " <<
                                 u.inversion << " " <<
                                 u.primaryBlockId << " " <<
                                 u.secondaryBlockId << " " << std::endl;
        }

        // Pring Nuc mutations
        std::cout << "Nuc mutations" << std::endl;
        for(auto &u: node->nucMutation) {
            std::cout << "\t Position " << u.nucPosition << " Gap-position " <<
                                 u.nucGapPosition << " " <<
                                 printNucMut(u.mutInfo) << " " <<
                                 printNucs(u.mutInfo, u.nucs) << " " << std::endl;
        }
        node = node->parent;
    }

    return;
}



// program option description for building/loading a PanMAT into memory
po::options_description globalDesc("panmanUtils Command Line Arguments");
po::positional_options_description globalPositionArgumentDesc;

// program option descriptions of various command line functions
po::options_description summaryDesc("Summary Command Line Arguments");
po::options_description useDesc("Use Command Line Arguments");
po::options_description imputeDesc("Imputation Command Line Arguments");
po::options_description fastaDesc("FASTA Command Line Arguments");
po::positional_options_description fastaPositionArgumentDesc;
po::options_description fastaAlignDesc("FASTA Command Line Arguments");
po::positional_options_description fastaAlignPositionArgumentDesc;
po::options_description subnetDesc("Subnetwork Command Line Arguments");
po::positional_options_description subnetPositionArgumentDesc;
po::options_description vcfDesc("VCF writer Command Line Arguments");
po::positional_options_description vcfPositionArgumentDesc;
po::options_description gfaDesc("GFA writer Command Line Arguments");
po::positional_options_description gfaPositionArgumentDesc;
po::options_description mafDesc("MAF Writer Command Line Arguments");
po::positional_options_description mafPositionArgumentDesc;
po::options_description newickDesc("Newick Writer Command Line Arguments");
po::positional_options_description newickPositionArgumentDesc;
po::options_description extendNewickDesc("Extended Newick Writer Command Line Arguments");
po::positional_options_description extedNewickDescPositionArgumentDesc;
po::options_description annotateDesc("Annotate Command Line Arguments");
po::positional_options_description annotatePositionArgumentDesc;
po::options_description rerootDesc("Reroot Command Line Arguments");
po::positional_options_description rerootArgumentDesc;
po::options_description aaDesc("Amino Acid Translation Command Line Arguments");
po::positional_options_description aaTranslationArgumentDesc;
po::options_description createNetDesc("Create Network Command Line Arguments");
po::positional_options_description createNetDescTranslationArgumentDesc;
po::options_description printMutDesc("Print Mutations Command Line Arguments");
po::positional_options_description printMutPositionArgumentDesc;
po::options_description printPathDesc("Print Paths Command Line Arguments");
po::positional_options_description printPathsArgumentDesc;
po::options_description indexDesc("Indexing Command Line Arguments");
po::positional_options_description indexArgumentDesc;
po::options_description printRootDesc("Root Printer Command Line Arguments");
po::positional_options_description printRootPositionArgumentDesc;


void setupOptionDescriptions() {
    // Global option descriptions
    globalDesc.add_options()
    ("help,h", "Print help messages")
    ("input-panman,I", po::value< std::string >(), "Input PanMAN file path")
    // ("input-panmat,T", po::value< std::string >(), "Input PanMAT file path")
    ("input-pangraph,P", po::value< std::string >(), "Input PanGraph JSON file to build a PanMAN")
    ("input-gfa,G", po::value< std::string >(), "Input GFA file to build a PanMAN")
    ("input-msa,M", po::value< std::string >(), "Input MSA file (FASTA format) to build a PanMAN")
    ("input-newick,N", po::value< std::string >(), "Input tree topology as Newick string")
    ("impute", "Create new PanMAN with N sequences inputed")
    ("create-network,K",po::value< std::vector<std::string>>(), "Create PanMAN with network of trees from single or multiple PanMAN files")

    // ("optimize", "currently UNSUPPORTED: whether given msa file should be optimized or not")

    ("printTips", po::value< std::string >(),"Print PanMAN summary")
    ("summary,s", "Print PanMAN summary")
    ("newick,t", "Print newick string of all trees in a PanMAN")
    ("fasta,f", "Print tip sequences (FASTA format)")
    // ("fasta-fast", "Print tip/internal sequences (FASTA format)")
    ("fasta-aligned,m", "Print MSA of sequences for each PanMAT in a PanMAN (FASTA format)")
    ("subnet,b", "Extract subnet of given PanMAN to a new PanMAN file based on the list of nodes provided in the input-file")
    ("vcf,v", "Print variations of all sequences from any PanMAT in a PanMAN (VCF format)")
    ("gfa,g", "Convert any PanMAT in a PanMAN to a GFA file")
    ("maf,w", "Print m-WGA for each PanMAT in a PanMAN (MAF format)")
    ("annotate,a", "Annotate nodes of the input PanMAN based on the list provided in the input-file (TSV)")
    ("reroot,r", "Reroot a PanMAT in a PanMAN based on the input sequence id (--reference)")
    ("aa-translation,v", "Extract amino acid translations in TSV file")
    ("extended-newick,e", "Print PanMAN's network in extended-newick format")
    ("printMutations,p", "Print mutations from root to each node")
    ("acr,q", "ACR method [fitch(default), mppa]")
    ("index",po::value< bool >(0), "Generating indexes and print sequence (passed as reference) between x:y")
    ("refFile",po::value< std::string >() ,"reference sequence file")
    // ("printRoot", "Print root sequence")
    // ("printNodePaths", "Print mutations from root to each node")
    ("toUsher", "Convert a PanMAT in PanMAN to Usher-MAT")
    // ("protobuf2capnp", "Converts a Google Protobuf PanMAN to Capn' Proto PanMAN")
  
    ("low-mem-mode", "Perform Fitch Algrorithm in batch to save memory consumption")
    ("reference,n", po::value< std::string >(), "Identifier of reference sequence for PanMAN construction (optional), VCF extract (required), or reroot (required)")
    ("start,x", po::value< int64_t >(), "Start coordinate of protein translation/Start coordinate for indexing")
    ("end,y", po::value< int64_t >(), "End coordinate of protein translation/End coordinate for indexing")
    ("treeID,d", po::value< std::string >(), "Tree ID, required for --vcf")
    // ("tree-group", po::value< std::vector< std::string > >()->multitoken(), "File paths of PMATs to generate tree group")
    ("input-file,i", po::value< std::string >(), "Path to the input file, required for --subnet, --annotate, and --create-network")
    ("output-file,o", po::value< std::string >(), "Prefix of the output file name")
    ("threads", po::value< std::int32_t >(), "Number of threads")
    // ("complexmutation-file", po::value< std::string >(), "File path of complex mutation file for tree group")

    // ("panman-in", po::value< std::string >(), "Input file path for PanMAT Group")

    ;

    // Adding input file as positional argument (doesn't require the --input-file tag)
    globalPositionArgumentDesc.add("input-panman", -1);

    // Use option descriptions
    useDesc.add_options()
    ("help", "produce help message")
    ("index", po::value< size_t >()->required(), "PanMAT index")
    ;

    imputeDesc.add_options()
        ("input-file", po::value< std::string >(), "Input file name")
        ("output-file,o", po::value< std::string >(), "Output file name");

    summaryDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");

    // FASTA option descriptions
    fastaDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");
    
    fastaAlignDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");

    subnetDesc.add_options()
        ("input-file", po::value< std::string >(), "Input file name")
        ("output-file,o", po::value< std::string >(), "Output file name");

    vcfDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("reference", po::value< std::string >(), "Reference name")
        ("output-file,o", po::value< std::string >(), "Output file name");

    gfaDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("output-file,o", po::value< std::string >(), "Output file name");

    // MAF option descriptions
    mafDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("output-file,o", po::value< std::string >(), "Output file name");

    newickDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");

    extendNewickDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");

    annotateDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("input-file", po::value< std::string >(), "Input file name")
        ("output-file,o", po::value< std::string >(), "Output file name");

    rerootDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("reference", po::value< std::string >(), "Reference name")
        ("output-file,o", po::value< std::string >(), "Output file name");

    aaDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("start,s", po::value< int64_t >(), "Start coordinate of protein translation/Start coordinate for indexing")
        ("end,e", po::value< int64_t >(), "End coordinate of protein translation/End coordinate for indexing")
        ("output-file,o", po::value< std::string >(), "Output file name");

    createNetDesc.add_options()
        ("input-file", po::value< std::string >(), "File containing complex mutations")
        ("tree-group", po::value< int64_t >(), "List of PanMATs")
        ("output-file,o", po::value< std::string >(), "Output file name");
    
    printMutDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("output-file,o", po::value< std::string >(), "Output file name");
        
    printPathDesc.add_options()
        ("treeID", po::value< std::int64_t >(), "Tree ID [default 0]")
        ("output-file,o", po::value< std::string >(), "Output file name");
    
    indexDesc.add_options()
        ("reference", po::value< std::string >(), "Reference name")
        ("start,s", po::value< int64_t >(), "Start coordinate of protein translation/Start coordinate for indexing")
        ("end,e", po::value< int64_t >(), "End coordinate of protein translation/End coordinate for indexing")
        ("index",po::value< bool >(0), "Generating indexes and print sequence (passed as reference) between x:y")
        ("output-file,o", po::value< std::string >(), "Output file name");

    printRootDesc.add_options()
        ("output-file,o", po::value< std::string >(), "Output file name");
}

void writePanMAN(po::variables_map &globalVm, panmanUtils::TreeGroup *TG) {
    std::cout << "Writing PanMAN" << std::endl;
    std::string fileName = globalVm["output-file"].as< std::string >();
    std::filesystem::create_directory("./panman");

    std::ofstream outputFile("./panman/" + fileName + ".panman");
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

    auto writeStart = std::chrono::high_resolution_clock::now();

    // outPMATBuffer.push(boost::iostreams::gzip_compressor());
    boost::iostreams::lzma_params params;
    params.level = 9; // Highest compression level
    outPMATBuffer.push(boost::iostreams::lzma_compressor(params));
    outPMATBuffer.push(outputFile);
    std::ostream outstream(&outPMATBuffer);

    kj::std::StdOutputStream outputStream(outstream);

    TG->writeToFile(outputStream);
    boost::iostreams::close(outPMATBuffer);
    outputFile.close();

    auto writeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds writeTime = writeEnd - writeStart;
    std::cout << "\nNetwork Write execution time: " << writeTime.count()
              << " nanoseconds\n";

}

void writePanMAN(po::variables_map &globalVm, panmanUtils::Tree *T) {
    std::string fileName = globalVm["output-file"].as< std::string >();
    std::filesystem::create_directory("./panman");

    std::ofstream outputFile("./panman/" + fileName + ".panman");
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

    auto writeStart = std::chrono::high_resolution_clock::now();

    // outPMATBuffer.push(boost::iostreams::gzip_compressor());
    boost::iostreams::lzma_params params;
    params.level = 9; // Highest compression level
    outPMATBuffer.push(boost::iostreams::lzma_compressor(params));
    outPMATBuffer.push(outputFile);
    std::ostream outstream(&outPMATBuffer);
    kj::std::StdOutputStream outputStream(outstream);
    T->writeToFile(outputStream);
    boost::iostreams::close(outPMATBuffer);
    outputFile.close();

    auto writeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds writeTime = writeEnd - writeStart;
    std::cout << "\nNetwork Write execution time: " << writeTime.count()
              << " nanoseconds\n";

}

void impute(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // If command was impute, create a new PanMAN with imputed sequences
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    std::string fileName = globalVm["output-file"].as< std::string >();

    panmanUtils::TreeGroup tg = *TG;

    auto imputeStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        tg.trees[i].imputeNs();
    }

    auto imputeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds imputeTime = imputeEnd - imputeStart;
    std::cout << "\nImputation time: " << imputeTime.count() << " nanoseconds\n";

    writePanMAN(globalVm, &tg);
}

void summary(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // If command was summary, print the summary of the PanMAT
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto summaryStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree *T = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".summary");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);
        T->printSummary(fout);

        if(globalVm.count("output-file")) outputFile.close();
    }

    auto summaryEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;
    std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
}

void fasta(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print raw sequences to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree *T  = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".fasta");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        T->printFASTAUltraFast(fout, false, false);

        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";
}

void printTipsHelper(panmanUtils::Node* node) {
    if (node == nullptr) return;
    if (node->children.size() == 0) {
        std::cout << node->identifier << std::endl;
    } else {
        for (auto child: node->children) {
            printTipsHelper(child);
        }
    }
}

void printTips(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print raw sequences to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree *T  = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".fasta");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        std::string tip = globalVm["printTips"].as< std::string >();
        auto node = T->allNodes[tip];
        printTipsHelper(node);
        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";
}


void fastaAligned(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print multiple sequence alignment to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree *T  = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".msa");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);


        T->printFASTAUltraFast(fout, true);


        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";
}

void fastaFast(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print multiple sequence alignment to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree *T  = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".fasta");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);


        T->printFASTAUltraFast(fout, false, false);


        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";
}

void subnetwork(panmanUtils::Tree *T, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Extract subnet of PanMAN to new file

    if(T == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }
    if(!globalVm.count("output-file")) {
        panmanUtils::printError("Output file not provided!");
        std::cout << globalDesc;
        return;
    }

    // List of node identifiers that need to be extracted from the tree
    std::vector< std::string > nodeIds;
    std::string nodeId;

    if(globalVm.count("input-file")) {
        std::string inputFileName = globalVm["input-file"].as< std::string >();
        std::ifstream fin(inputFileName);
        while(fin >> nodeId) {
            nodeIds.push_back(nodeId);
        }
        fin.close();
    } else {
        panmanUtils::printError("No source of node ids provided");
        exit(0);
    }

    if(nodeIds.size() == 0) {
        std::cout << "No node identifiers provided!" << std::endl;
    }

    std::string outputFileName = globalVm["output-file"].as< std::string >();
    std::filesystem::create_directory("./panman");
    std::ofstream outputFiles("./panman/" + outputFileName + ".panman");
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

    auto subtreeStart = std::chrono::high_resolution_clock::now();

    // outPMATBuffer.push(boost::iostreams::gzip_compressor());
    boost::iostreams::lzma_params params;
    params.level = 9; // Highest compression level
    outPMATBuffer.push(boost::iostreams::lzma_compressor(params));
    outPMATBuffer.push(outputFiles);
    std::ostream outstream(&outPMATBuffer);
    kj::std::StdOutputStream outputStream(outstream);
    T->writeToFile(outputStream, T->subtreeExtractParallel(nodeIds));
    boost::iostreams::close(outPMATBuffer);
    outputFiles.close();

    auto subtreeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

    std::cout << "\nParallel Subtree Extract execution time: "
                << subtreeTime.count() << " nanoseconds\n";
}

void subnet(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Extract the subnetwork consisting of given node IDs from PanMAN

    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    std::string outputFileName;
    if(!globalVm.count("output-file")) {
        panmanUtils::printError("Output file not provided!");
        std::cout << globalDesc;
        return;
    } else outputFileName = globalVm["output-file"].as< std::string >();

    // List of node identifiers that need to be extracted from the tree
    std::unordered_map< int, std::vector< std::string > > nodeIds;
    std::string nodeId;

    if(globalVm.count("input-file")) {
        std::string inputFileName = globalVm["input-file"].as< std::string >();
        std::ifstream fin(inputFileName);
        std::string line;
        int treeId;
        while(std::getline(fin, line)) {
            std::stringstream ss(line);
            ss >> treeId;
            while(ss >> nodeId) {
                nodeIds[treeId].push_back(nodeId);
            }
        }
        fin.close();
    } else {
        panmanUtils::printError("Input file not provided!");
        std::cout << globalDesc;
        return;
    }

    if(nodeIds.size() == 0) {
        std::cout << "No node identifiers selected!" << std::endl;
    }

    std::filesystem::create_directory("./panman");
    std::ofstream outputFiles("./panman/" + outputFileName + ".panman");
    boost::iostreams::filtering_streambuf< boost::iostreams::output>
    outPMATBuffer;

    auto subtreeStart = std::chrono::high_resolution_clock::now();

    // outPMATBuffer.push(boost::iostreams::gzip_compressor());
    boost::iostreams::lzma_params params;
    params.level = 9; // Highest compression level
    outPMATBuffer.push(boost::iostreams::lzma_compressor(params));
    outPMATBuffer.push(outputFiles);
    std::ostream outstream(&outPMATBuffer);
    kj::std::StdOutputStream outputStream(outstream);
    panmanUtils::TreeGroup* subnetwork = TG->subnetworkExtract(nodeIds);
    subnetwork->writeToFile(outputStream);

    boost::iostreams::close(outPMATBuffer);
    outputFiles.close();

    auto subtreeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

    std::cout << "\nParallel Subnetwork Extract execution time: "
                << subtreeTime.count() << " nanoseconds\n";
}

void vcf(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }
    int treeID = 0;
    if(globalVm.count("treeID")) treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    std::string reference;
    if(!globalVm.count("reference")) {
        for (auto &n: T->allNodes) {
            reference = n.first;
            break;
        }
    } else reference = globalVm["reference"].as< std::string >();

    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".vcf");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);

    auto vcfStart = std::chrono::high_resolution_clock::now();

    panmanUtils::Node* refNode;
    for (auto &n: T->allNodes) {
        if (n.first == reference) {
            refNode = n.second;
            break;
        }
    }
    T->printVCFParallel(refNode, fout);

    auto vcfEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;
    std::cout << "\nVCF execution time: " << vcfTime.count() << " nanoseconds\n";
    if(globalVm.count("output-file")) outputFile.close();
}

void gfa(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // If GFA is to be extracted from PanMAN

    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    int treeID = 0;
    if(globalVm.count("treeID")) treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".gfa");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);

    auto generateVGStart = std::chrono::high_resolution_clock::now();

    T->convertToGFA(fout);

    auto generateVGEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds generateVGTime = generateVGEnd - generateVGStart;

    std::cout << "GFA generation time: " << generateVGTime.count()
                << " nanoseconds\n";
    if(globalVm.count("output-file")) outputFile.close();
}

void maf(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    if(TG == nullptr) {
        std::cout << "No PanMAN selected. Try groupFasta for FASTA of the whole"
                    " PanMAN" << std::endl;
        return;
    }

    int treeID = 0;
    if(globalVm.count("treeID")) treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".maf");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);

    auto mafStart = std::chrono::high_resolution_clock::now();

    T->printMAF(fout);

    auto mafEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mafTime = mafEnd - mafStart;
    std::cout << "\nMAF execution time: " << mafTime.count() << " nanoseconds\n";
    if(globalVm.count("output-file")) outputFile.close();
}

void newick (panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print newick string of the PanMAT or PanMAN loaded into memory
    if(TG) {
        int index = 0;
        for(auto& t: TG->trees) {
            if(globalVm.count("output-file")) {
                std::string fileName = globalVm["output-file"].as< std::string >();
                outputFile.open("./info/" + fileName + "_" + std::to_string(index) + ".newick");
                buf = outputFile.rdbuf();
            } else {
                buf = std::cout.rdbuf();
            }
            std::ostream fout (buf);
            fout << t.getNewickString(t.root) << std::endl;
            if(globalVm.count("output-file")) outputFile.close();
        }
    } else {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }
}

void extendNewick(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print Extended Newick String
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".extended-newick");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);


    auto writeStart = std::chrono::high_resolution_clock::now();

    for (auto& tree: TG->trees) {
        fout << tree.getNewickString(tree.root) << std::endl;
    }

    TG->printComplexMutations(fout);

    if(globalVm.count("output-file")) outputFile.close();

    auto writeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds writeTime = writeEnd - writeStart;
    std::cout << "\nExtended Newick execution time: " << writeTime.count()
                << " nanoseconds\n";
}

void annotate(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Annotate nodes of PanMAT
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    int treeID = 0;
    if(globalVm.count("treeID")) {
        treeID = std::stoi(globalVm["treeID"].as< std::string >());
    }

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    if(!globalVm.count("input-file")) {
        panmanUtils::printError("Input file not provided!");
        std::cout << globalDesc;
        return;
    }

    std::string fileName = globalVm["input-file"].as< std::string >();
    std::ifstream fin(fileName);
    auto annotateStart = std::chrono::high_resolution_clock::now();

    T->annotate(fin);

    auto annotateEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds annotateTime = annotateEnd - annotateStart;
    std::cout << "Annotate time: " << annotateTime.count() << " nanoseconds\n";

    writePanMAN(globalVm,TG);
}

void reroot(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Reroot the PanMAT to given sequence
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    int treeID;
    if(!globalVm.count("treeID")) {
        panmanUtils::printError("TreeID not provided!");
        std::cout << globalDesc;
        return;
    } else treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    if(!globalVm.count("reference")) {
        panmanUtils::printError("Refence ID not provided!");
        std::cout << globalDesc;
        return;
    }

    std::string sequenceName = globalVm["reference"].as< std::string >();

    auto rerootStart = std::chrono::high_resolution_clock::now();

    T->reroot(sequenceName);

    auto rerootEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds rerootTime = rerootEnd - rerootStart;
    std::cout << "\nReroot execution time: " << rerootTime.count()
                << " nanoseconds\n";

    TG->trees[treeID] = *T;


    writePanMAN(globalVm, TG);
}

void aa(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Extract amino acid translations in tsv file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    int treeID;
    if(!globalVm.count("treeID")) {
        panmanUtils::printError("TreeID not provided!");
        std::cout << globalDesc;
        return;
    } else treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &tg.trees[treeID];

    if(!globalVm.count("start") || !globalVm.count("end")) {
        std::cout << "Start/End Coordinate not provided" << std::endl;
        return;
    }

    int64_t startCoordinate = globalVm["start"].as< int64_t >();
    int64_t endCoordinate = globalVm["end"].as< int64_t >();

    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".tsv");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);

    auto aaStart = std::chrono::high_resolution_clock::now();

    T->extractAminoAcidTranslations(fout, startCoordinate, endCoordinate);

    auto aaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds aaTime = aaEnd - aaStart;
    std::cout << "\nAmino Acid translate execution time: " << aaTime.count()
                << " nanoseconds\n";
    if(globalVm.count("output-file")) outputFile.close();
}

void protobuf2capnp(panmanUtils::TreeGroup *TG, po::variables_map &globalVm) {
    std::string fileName = globalVm["input-panman"].as< std::string >();
    std::ifstream inputFile(fileName);
    boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::lzma_decompressor());
    inPMATBuffer.push(inputFile);
    std::istream inputStream(&inPMATBuffer);

    std::cout << "starting reading panman" << std::endl;
    TG = new panmanUtils::TreeGroup(inputStream, true);
    inputFile.close();

    writePanMAN(globalVm, TG);

}

void createNet(po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Create PanMAN from list of PanMAT files and a complex mutation file listing the complex
    // mutations relating these PanMATs

    std::vector< std::string > fileNames;

    std::string mutationFileName;
    if(!globalVm.count("input-file")) {
        panmanUtils::printError("Input File containing complex mutations not provided!");
        return;
    }

    fileNames = globalVm["create-network"].as< std::vector< std::string > >();
    mutationFileName = globalVm["input-file"].as< std::string >();

    std::ifstream mutationFile(mutationFileName);

    std::vector< std::ifstream > files;
    for(auto u: fileNames) {
        files.emplace_back(u);
    }

    auto treeBuiltStart = std::chrono::high_resolution_clock::now();

    // Currently handle only one file
    boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::lzma_decompressor());
    inPMATBuffer.push(files[0]);
    std::istream inputStream(&inPMATBuffer);
    panmanUtils::TreeGroup* TG = new panmanUtils::TreeGroup(inputStream);
    

    std::vector< panmanUtils::Tree* > tg;
    for (int i=0; i<TG->trees.size(); i++) {
        tg.push_back(&TG->trees[i]);
    }


    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

    panmanUtils::TreeGroup* TG_new = new panmanUtils::TreeGroup(tg, mutationFile);
    mutationFile.close();
    for(auto& u: files) {
       u.close();
    }

    writePanMAN(globalVm,TG_new);

}

void printMut(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    int treeID = 0;
    if(globalVm.count("treeID")) treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &TG->trees[treeID];
    // T = &tg.trees[treeID];


    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".mutations");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);


    auto substitutionsStart = std::chrono::high_resolution_clock::now();

    // std::cout << T->root->identifier << std::endl;

    // for (auto &n:T->allNodes){
    //     if (n.second->children.size() == 0) {
    //         std::cout << n.second->identifier << std::endl;
    //     }
    // }

    if (globalVm.count("input-file")) {
        std::string inputFileName = globalVm["input-file"].as< std::string >();
        std::ifstream fin(inputFileName);
        std::vector<std::string> nodeIds ={};
        // = {"node_841140", "node_225333", "node_1049429", "node_1691143", "node_99313"};
        // std::vector<std::string> nodeIds = {"node_3", "node_14","node_66", "node_140","node_261","node_341","node_431","node_1633","node_1698","node_1708187","node_1708254","node_1709987","node_1710263","node_29","node_81","node_150","node_200","node_445","node_1877","node_1707795","node_1707968","node_1708228","node_1708240","node_1708265","node_1708265","node_448","node_564","node_1119","node_1483","node_1504","node_1565","node_1708017","node_1708235","node_1709934","node_501","node_565","node_1136","node_1158","node_1267","node_1440","node_1465","node_1896","node_2006","node_2006","node_912","node_1278","node_1611","node_2028","node_2549","node_2554","node_3581","node_3603","node_3665","node_3671","node_3741","node_4002","node_4022","node_4309","node_4313","node_4627","node_4737","node_5767","node_5785","node_6007","node_6026","node_6046","node_6053","node_9217","node_1649541","node_1707067","node_1707135","node_1707184","node_1707189","node_1707199","node_1707219","node_1707261","node_1709169","node_1709394","node_1709396","node_1709458","node_1324","node_1666","node_1912","node_1938","node_2008","node_2584","node_2636","node_2641","node_2664","node_2670","node_2690","node_2804","node_2823","node_2854","node_2945","node_2970","node_3216","node_3337","node_3392","node_3483","node_3538","node_3746","node_3772","node_4080","node_4116","node_4134","node_4333","node_4514","node_4584","node_4641","node_4663","node_5792","node_9251","node_9714","node_9737","node_9750","node_9884","node_10068","node_10084","node_10092","node_10124","node_10155","node_32223","node_39456","node_70433","node_71745","node_72072","node_72110","node_72147","node_72147","node_1649668","node_1694127","node_1694131","node_1694132","node_1695137","node_1695211","node_1695484","node_1695533","node_1706913","node_1706980","node_1706990","node_1707001","node_1707237","node_1709151","node_1709972","node_1923","node_3011","node_3468","node_3831","node_4239","node_4548","node_5733","node_5815","node_9267","node_9745","node_39470","node_39863","node_39871","node_39884","node_40398","node_40420","node_40481","node_40500","node_40833","node_40889","node_67452","node_67475","node_67481","node_67520","node_67536","node_67548","node_67601","node_67932","node_68014","node_68251","node_68284","node_70217","node_70288","node_72498","node_72547","node_72553","node_80669","node_80860","node_81071","node_81278","node_82181","node_82776","node_82822","node_82949","node_83806","node_84133","node_84201","node_84206","node_84229","node_84390","node_84547","node_88313","node_88352","node_88367","node_88391","node_88434","node_88683","node_89457","node_89469","node_91185","node_91263","node_91298","node_91414","node_91437","node_91888","node_92676","node_92733","node_92932","node_95442","node_96451","node_96472","node_97523","node_97646","node_98033","node_98071","node_98158","node_98180","node_98240","node_98243","node_98279","node_98351","node_98608","node_1047501","node_1049252","node_1058569","node_1058961","node_1059029","node_1059072","node_1059240","node_1059861","node_1059938","node_1060087","node_1060100","node_1060208","node_1060282","node_1060304","node_1060414","node_1060605","node_1060756","node_1060891","node_1060904","node_1060979","node_1061013","node_1061016","node_1061051","node_1061108","node_1061132","node_1061134","node_1061137","node_1061155","node_1061212","node_1649561","node_1663677","node_1663715","node_1663773","node_1663804","node_1663883","node_1663923","node_1663974","node_1666325","node_1666524","node_1666588","node_1666654","node_1667081","node_1667355","node_1667360","node_1667483","node_1667540","node_1667808","node_1668591","node_1668647","node_1668649","node_1668663","node_1668911","node_1668953","node_1682190","node_1683647","node_1686835","node_1686996","node_1687206","node_1687467","node_1687569","node_1688811","node_1689701","node_1689707","node_1690784","node_1690906","node_1691020","node_1691143","node_1691185","node_1691218","node_1691267","node_1691344","node_1691371","node_1691378","node_1695162","node_1695323","node_1695398","node_1695968","node_1696108","node_1696208","node_1696335","node_1699125","node_1699399","node_1699832","node_1702928","node_1703154","node_1703242","node_1704193","node_1704206","node_1704293","node_1704398","node_1705201","node_1705641","node_1707357","node_1708288","node_1709022","node_1709081","node_1709174","node_1709182","node_1709189","node_1709385","node_2192","node_2590","node_10240","node_32254","node_32382","node_32396","node_32401","node_32436","node_32513","node_32562","node_37508","node_39394","node_39397","node_39466","node_39916","node_39973","node_40189","node_40851","node_67390","node_67578","node_67697","node_67753","node_67961","node_71808","node_71926","node_72102","node_72385","node_72528","node_72596","node_72679","node_72814","node_72846","node_73781","node_73790","node_78552","node_79634","node_79757","node_79834","node_80023","node_80043","node_80211","node_80249","node_80269","node_80486","node_80545","node_80576","node_80645","node_80680","node_80700","node_80715","node_80802","node_80895","node_81247","node_81273","node_81542","node_81558","node_82320","node_82400","node_82412","node_82434","node_82465","node_82759","node_82781","node_82868","node_82893","node_84140","node_84168","node_84257","node_84292","node_84537","node_87774","node_88411","node_88539","node_88584","node_88940","node_89364","node_89554","node_89574","node_90437","node_90451","node_91029","node_91273","node_91304","node_91472","node_91911","node_91925","node_91934","node_92131","node_92132","node_92150","node_92173","node_92747","node_92816","node_92996","node_94270","node_94408","node_94502","node_94516","node_94557","node_95359","node_95413","node_95427","node_95469","node_95581","node_95690","node_95729","node_95763","node_95782","node_95835","node_96015","node_96312","node_96445","node_96476","node_96507","node_96745","node_97529","node_97550","node_97562","node_97708","node_98019","node_98038","node_98140","node_98147","node_98181","node_98296","node_98584","node_98639","node_98666","node_98693","node_98736","node_224891","node_225095","node_225125","node_225252","node_225259","node_1047523","node_1047808","node_1047832","node_1047903","node_1048110","node_1048170","node_1048200","node_1048204","node_1048270","node_1048392","node_1048441","node_1048447","node_1048454","node_1048494","node_1048650","node_1049367","node_1058049","node_1058061","node_1058068","node_1058108","node_1058127","node_1058156","node_1058917","node_1059115","node_1059295","node_1059544","node_1059841","node_1059857","node_1059933","node_1060093","node_1060328","node_1060734","node_1060758","node_1060913","node_1060981","node_1061091","node_1662746","node_1663057","node_1663063","node_1663395","node_1663448","node_1663548","node_1663562","node_1663596","node_1663793","node_1663940","node_1663975","node_1666311","node_1666344","node_1666412","node_1666603","node_1666892","node_1667009","node_1667103","node_1667320","node_1667362","node_1667514","node_1667598","node_1668397","node_1668443","node_1668470","node_1668494","node_1683271","node_1683497","node_1683757","node_1683821","node_1683831","node_1683837","node_1683884","node_1685795","node_1685982","node_1686086","node_1686719","node_1686748","node_1686933","node_1686990","node_1687213","node_1687504","node_1687702","node_1687723","node_1687889","node_1689636","node_1689648","node_1689726","node_1689896","node_1690138","node_1690159","node_1690399","node_1690829","node_1691482","node_1695999","node_1696151","node_1698877","node_1702817","node_1704021","node_1704769","node_1707269","node_1708855","node_1709248","node_3501","node_31113","node_32549","node_32668","node_32729","node_37510","node_40261","node_40545","node_67358","node_67644","node_67754","node_67759","node_67763","node_67826","node_70311","node_72700","node_72724","node_72730","node_72848","node_78398","node_78416","node_78935","node_79529","node_79684","node_79719","node_80214","node_80561","node_80656","node_82478","node_82569","node_82649","node_82688","node_82712","node_82861","node_82918","(Theta)","node_83756","node_83954","node_84060","node_84254","node_84616","node_87253","node_87768","node_88178","node_88400","node_88416","node_88630","node_88823","node_88831","node_90487","node_91126","node_91388","node_91465","node_92110","node_92636","node_92654","node_92912","node_94549","node_95584","node_95751","node_95859","node_95945","node_96333","node_96553","node_97731","node_98171","node_98421","node_225035","node_225209","node_225236","node_1047584","node_1047595","node_1047684","node_1047732","node_1047931","node_1048120","node_1048182","node_1048716","node_1049300","node_1049326","node_1049342","node_1049349","node_1049403","node_1058077","node_1058316","node_1058735","node_1058746","node_1058833","node_1058875","node_1058894","node_1059003","node_1059495","node_1059512","node_1059669","node_1060339","node_1060684","node_1061186","node_1663013","node_1663370","node_1663437","node_1663538","node_1665054","node_1667054","node_1667525","node_1667697","node_1668502","node_1668722","node_1683447","node_1683451","node_1684442","node_1685040","node_1685216","node_1685485","node_1685733","node_1686766","node_1687472","node_1687631","node_1687661","node_1687691","node_1687912","node_1688064","node_1690790","node_1690791","node_1690795","node_1690803","node_1691195","node_1691273","node_1691464","node_1695884","node_1696867","node_1697100","node_1698601","node_1698772","node_1702836","node_1703118","node_1703946","node_1704407","node_1704962","node_1708321","node_1709093","node_1709207","node_2722","node_5621","node_32313","node_32329","node_32528","node_37511","node_37549","node_37601","node_38575","node_38828","node_39170","node_39191","node_39244","node_40536","node_40563","node_67733","node_67741","node_72151","node_72355","node_72399","node_79200","node_79554","node_79696","node_82955","node_83047","node_83912","node_84043","node_84617","node_87974","node_88606","node_89030","node_89105","node_89490","node_91091","node_91875","node_92600","node_92957","node_93175","node_93175","node_95637","node_95681","node_95694","node_97663","node_97696","node_98131","node_98641","node_225261","node_1048475","node_1049269","node_1049314","node_1058754","node_1058828","node_1058884","node_1059023","node_1059521","node_1059633","node_1060142","node_1060155","node_1060312","node_1060451","node_1060618","node_1060763","node_1649759","node_1649839","node_1662749","node_1663073","node_1663086","node_1663113","node_1663152","node_1663169","node_1663234","node_1663323","node_1666667","node_1667160","node_1667419","node_1667646","node_1669038","node_1683353","node_1683507","node_1683977","node_1684383","node_1684395","node_1686490","node_1686689","node_1688080","node_1689922","node_1690842","node_1694921","node_1704520","node_1706544","node_1706562","node_1708291","node_1708950","node_1709027","England/PHEC-31F487/2021|OX863249.1|2021-06-15","node_1709373","node_5875","node_10301","England/CAMC-11ABD0C/2021|OX605871.1|2021-01-27","node_37777","node_38529","node_38821","node_38904","node_38920","node_39044","node_39177","node_39278","node_70553","node_71919","node_80582","node_81179","node_83985","node_84658","node_84809","node_84821","node_84901","node_84988","node_84989","node_85646","node_85936","node_85940","node_86105","node_86117","node_86131","node_86967","node_87096","node_87097","node_87117","node_87163","node_87878","node_88467","node_88951","node_89131","node_92990","node_96758","node_97947","node_98376","node_98473","node_1058442","node_1058613","node_1058751","node_1059174","node_1059589","node_1662761","node_1663275","node_1663283","node_1663304","node_1666736","node_1667333","node_1667817","node_1683279","node_1683896","node_1685154","node_1685375","node_1686190","node_1688828","node_1690555","node_1693946","node_1703905","node_1708347","node_4756","node_10572","node_10830","node_13473","node_14049","(Eta)","node_31218","node_37922","node_38399","node_38621","node_39997","node_40938","node_67701","node_68325","node_78909","node_80226","node_80623","node_81602","node_83036","node_83961","node_84942","node_85020","node_85976","node_86995","node_94317","node_94570","node_1048293","node_1049184","node_1049429","node_1057971","node_1059123","node_1059261","node_1649353","node_1663267","node_1664579","node_1690175","node_1691727","node_1694321","node_1696519","node_5552","node_10542","node_11159","node_11311","node_11466","node_11761","node_12203","node_13489","node_16384","node_16531","node_17844","node_18083","node_29833","node_30767","node_32369","node_38419","node_38710","node_38964","node_41007","node_41481","node_41481","node_68326","(Mu)","node_81335","node_81370","node_86979","node_87101","node_95345","node_1048582","(Delta)","node_1061825","node_1691716","node_1692160","node_1702102","node_1708349","node_1708517","(Beta)","node_6102","node_12016","node_12192","node_13710","node_15231","node_16646","node_16740","node_16791","node_17806","node_17821","node_17847","node_18057","node_18118","node_29864","node_30750","node_30965","node_68352","node_68362","node_83079","node_85076","node_85790","node_85801","(Lambda)","node_87178","(Omicron)","(Epsilon)","node_1649859","node_1666370","node_1688360","node_1691721","node_1703761","node_1708481","node_4840","node_11585","node_13067","node_14139","node_15063","node_15204","node_15537","node_15923","node_17890","node_18031","node_18495","node_19402","node_19489","node_22387","node_23195","node_23260","node_23857","node_23942","node_29794","node_29909","node_32738","node_37781","node_70559","node_72901","node_85036","node_86139","(Gamma)","(Kappa)","node_1061426","node_1684105","node_1692162","node_1692199","node_1692576","node_1692618","node_1692653","node_13068","node_14435","node_14976","node_14989","node_15037","node_15041","node_15443","node_15566","node_16845","node_18163","node_18740","node_19563","node_20261","node_21782","node_22199","node_29918","node_70671","node_71171","node_71591","node_74202","node_79411","node_92346","node_99279","node_1678876","node_1684036","node_1692297","node_1693466","node_1693522","node_1693548","node_18280","node_18299","node_20781","node_21577","node_23898","node_24148","node_25295","node_38276","node_1649864","node_1669167","node_1693582","node_12055","node_14909","node_16679","node_22017","node_68820","node_85117","node_85843","node_96813","(Alpha)","node_841137","node_1665075","node_1682462","node_8472","node_8552","node_13303","node_13314","node_13336","node_13356","node_23870","node_225274","(BA.1)","node_841140","node_13090","node_13276","node_13294","node_16236","node_27562","node_86725","node_840821","node_1049616","node_1051876","node_1052953","node_1054203","node_1054370","node_1056244","node_1057086","(Iota)","node_37385","node_68376","node_96816","node_99313","node_225276","node_840589","node_840813","node_840931","Switzerland/VD-Risch-1202R15118/2022|OX388819.1|2022-12-02","node_1052136","node_1053136","node_1055028","node_1055443","node_1055734","node_1056622","node_7668","node_32956","node_33561","node_37146","node_37206","node_168662","node_840590","node_840810","node_840823","node_841333","node_843358","node_844060","node_844930","node_846005","node_846978","node_847183","node_855941","node_856232","node_857120","node_857388","node_861147","node_914059","node_915251","node_1046710","node_1049677","BRA/21891980IRT/2021|OQ521783.1|2021-07-02","(Delta)","node_1665406","node_1681623","node_33407","node_33488","node_33844","node_34244","node_34483","node_35059","node_35873","node_36167","node_37331","node_68652","node_96926","node_97213","node_97229","node_225281","node_840836","node_840928","node_840935","node_840971","node_841040","node_847386","node_848089","node_856600","node_856702","node_863725","node_864940","node_891513","node_915252","node_1052763","node_1612037","node_1612306","node_1612506","node_1613212","node_1628045","node_1630355","node_1648240","node_1648684","node_1648775","node_1648905","node_1648989","node_1649071","node_34265","node_34464","node_35788","node_35926","node_36351","node_37188","node_85167","node_99204","node_155626","node_170473","(BA.2)","node_225333","node_840594","node_840696","node_840710","node_840715","node_840718","node_840921","node_847817","node_857297","node_894358","node_913931","node_925303","node_928000","node_930371","node_931720","node_938224","node_938904","node_946212","node_948102","node_953101","node_962846","node_965259","node_974226","node_980728","USA/CO-CDC-MMB12218639/2021|OM074549.1|2021-12-16","node_1049509","node_1611528","node_1619015","node_1629119","node_1629893","node_1634594","node_1635680","node_1636254","node_33253","node_35659","node_35927","node_36169","node_96928","node_146018","node_840725","node_840729","node_840804","node_840839","node_840907","node_858534","node_864946","node_905761","node_906670","node_916700","node_919354","node_920040","node_940395","node_941651","node_992760","node_1025396","node_1049644","node_1603793","node_1611567","node_1611785","node_1612638","node_1628306","node_1628387","node_1629218","node_1636112","node_1657808","node_212818","node_213133","node_225421","node_230172","node_235672","node_248458","node_256371","node_256886","node_277886","node_393088","node_467852","node_470564","node_749812","node_750600","node_757664","node_839569","node_840596","USA/MN-MDH-37557/2023|OR872416.1|2023-10-16","node_840726","node_840730","node_840860","node_857998","node_894783","node_1049540","node_1612529","node_1614104","node_1628115","node_1628573","node_1630394","node_1634964","node_1647995","node_36261","node_226537","node_232246","node_234492","node_236306","node_238001","node_247696","node_253686","node_256085","node_256160","node_268026","node_272916","node_273969","node_274996","node_275500","node_279082","node_281160","node_403240","node_424118","node_427825","node_432816","node_433177","node_433201","node_434115","node_465328","node_465614","node_469082","node_469122","node_472924","node_474007","node_748181","node_748354","node_749022","node_756018","node_757819","node_757939","node_783125","node_786736","node_787451","node_792499","node_793758","node_794896","node_799727","node_800498","node_801187","node_803161","node_826972","node_834769","node_840664","node_840699","node_840797","node_840841","node_947605","node_1025398","node_1630380","node_1648816","node_35135","node_199113","Japan/sb_ncgm_sars_cov_2_03164/2022|BS008771.1|2022-04-10","node_230583","node_256161","node_256909","node_265677","node_274669","node_278243","node_280297","node_395722","node_395820","node_396776","node_397045","node_397579","node_397637","node_401113","node_401370","node_428362","node_428761","node_430131","node_432080","node_466587","node_468092","node_472544","node_473241","node_720779","node_758078","node_760911","node_793341","node_798766","node_803163","node_817267","node_821544","node_831614","node_834449","node_836871","node_840732","node_1025405","(Delta)","node_1626515","node_225709","node_227407","node_247072","node_248436","node_268102","node_281578","node_393666","node_395128","node_395405","node_396616","node_396745","node_398104","USA/FL-CDC-STM-UUWZN24D2/2022|OP464721.1|2022-08-31","node_399805","node_400704","node_401831","node_404599","node_423802","node_430145","(BA.2.12.1)","node_434119","node_464579","node_758033","node_758347","node_759301","node_775611","node_793765","node_809475","node_833658","node_1616030","node_227680","node_227694","node_397817","USA/MI-CDC-STM-7ND9EVZVA/2022|OP794825.1|2022-09-27","node_401905","USA/IL-CDC-LC0802864/2022|OP165785.1|2022-07-16","node_465638","node_465942","node_465969","node_474301","node_761944","node_786194","node_816974","node_831688","node_836930","node_840667","node_1062140","node_1088464","node_1137365","node_1139106","node_1149162","node_1215563","node_1261446","node_1285873","node_1287704","node_1308538","node_1309780","node_1319391","node_1329700","node_1334828","node_1335012","node_1337169","node_1344616","node_1348905","node_1363458","node_1401712","node_1415287","node_1423232","node_1583252","node_1592857","node_1601570","node_1601675","node_1601920","node_1602619","node_227409","node_227476","(BA.2.75)","node_404601","node_465932","node_465970","node_474302","node_474375","node_763527","node_764536","node_1085315","node_1087980","node_1088006","node_1137429","node_1137579","node_1137844","node_1138064","node_1138808","node_1208177","node_1218905","node_1218987","node_1219008","node_1219140","node_1315617","node_1319173","node_1323816","node_1326503","node_1328747","node_1329660","node_1332447","node_1335658","node_1337897","node_1338925","node_1340739","node_1343142","node_1343201","node_1343376","node_1361611","node_1361869","node_1415332","node_1415878","node_1417102","node_1422563","node_1422964","node_1572534","node_1599292","node_1601387","node_1612240","node_215198","node_227638","node_227720","node_272840","node_281585","node_400678","node_404725","node_404737","node_411172","node_421323","node_421347","node_421785","node_441814","(BA.4)","node_474376","node_763600","node_764616","node_1208161","node_1214825","node_1215035","node_1215085","node_1261282","node_1306169","node_1332140","node_1335125","node_1343689","node_1351287","node_1357612","node_1362090","node_1410466","node_1415238","node_1420345","node_1451723","node_1477459","node_1540665","node_1589045","node_1590627","node_1591484","node_1592502","node_1592518","node_1597775","node_1598956","node_227731","node_400685","node_404615","BHR/T103377662/2022|OP727490.1|2022-10-16","node_404738","node_404761","node_405023","node_411173","node_411211","node_411215","node_411234","node_421308","node_421320","node_421325","node_421350","node_421792","node_422999","node_445723","node_445984","node_474377","node_496704","node_763549","node_763615","node_764336","node_764394","node_764449","node_764479","node_764591","node_899014","node_1063780","node_1075601","node_1078122","node_1085867","node_1138113","node_1245805","node_1261014","node_1284023","node_1289413","node_1306192","node_1308559","node_1323818","node_1409423","node_1411486","node_1412587","England/PHEC-3W082W69/2021|OX785446.1|2021-12-15","node_1488122","node_1588829","node_1595823","Denmark/DCGC-598988/2022|OY835208.1|2022-10-16","node_404616","node_404826","node_405109","node_411237","node_421521","node_421528","node_421648","node_421655","node_421717","node_421730","node_421800","node_421821","node_421969","England/LSPA-3E3EBD4/2022|OX061351.1|2022-06-01","node_464523","node_474406","node_476028","node_477227","node_478068","node_480789","node_483932","node_486521","(BA.5)","node_496747","(BA.2.86)","node_721378","node_763530","node_803621","node_899005","node_899015","node_899031","node_1068734","node_1138501","node_1205673","node_1214322","node_1263158","node_1272595","node_1274698","node_1277429","node_1277862","node_1277925","node_1281011","node_1289020","node_1306335","node_1332177","node_1343832","node_1407596","node_1433514","node_1466637","node_1471437","node_1475096","node_1496084","node_1598997","node_404763","node_404841","node_405110","node_407789","node_409624","node_410540","node_410859","node_411082","node_411089","node_411176","node_411240","node_421822","node_421970","node_422248","node_422312","node_422507","node_423032","node_423187","node_434742","node_485790","node_487557","node_763532","node_764423","node_764501","node_764540","node_899019","node_899021","node_1081104","node_1086154","node_1086367","node_1217936","node_1241130","node_1244250","Switzerland/GE-HUG-35842111/2021|OU981925.1|2021-10-18","node_1329116","node_1332590","node_1335196","node_1343471","node_1412493","node_1419807","node_1496389","node_1497120","node_1526530","node_1561770","node_1568771","node_1595833","node_1636497","node_1637620","node_1638084","node_1638205","node_1638604","node_1638858","node_1638985","node_1641446","node_1643723","node_1643767","node_1643794","node_1645421","node_399814","node_408687","node_409569","node_410309","node_410935","node_419684","node_419697","node_421357","node_421405","node_421498","node_421541","node_422003","node_422388","node_422508","node_422526","node_423020","node_423028","node_423132","node_423214","node_475800","node_481234","node_484749","node_485212","node_500655","node_572678","node_613887","node_719621","node_719763","node_721380","node_810644","node_1142670","node_1217612","node_1249295","node_1278241","node_1324037","node_1325704","node_1329117","node_1331531","node_1353430","node_1373956","node_1435254","node_1435342","node_1508454","node_1638264","(XBB)","node_281590","node_399818","node_400187","node_400198","node_400537","node_400617","node_404846","node_408177","node_408414","node_409060","node_409101","node_409414","node_409504","node_410310","node_411188","node_411242","node_419643","node_419700","England/LSPA-3259EE4B/2022|OX360360.1|2022-10-03","node_422663","node_422956","node_422963","node_422988","node_423215","node_476783","node_481803","node_481861","node_487856","node_490254","node_490314","node_496125","node_500656","node_564428","node_564546","node_564868","node_572679","node_623457","node_720378","node_721492","node_721519","node_747936","node_747976","node_747988","node_747989","node_1333130","node_1342500","node_1435535","node_1444504","node_1507384","node_1511703","node_1644044","node_1647237","node_281591","node_384644","node_399832","node_399994","node_400055","node_400190","node_400204","node_400216","node_400242","node_400468","node_400581","node_405127","node_405480","node_406396","node_406538","node_407671","node_407928","node_408049","node_409317","node_409380","(CH.1.1)","node_411243","node_419644","node_421147","node_421213","node_422934","node_423026","Denmark/DCGC-552232/2022|OX289642.1|2022-07-18","node_481820","node_487883","node_497822","node_498102","node_499953","node_558829","node_565962","node_572140","node_572772","node_575024","node_577409","node_583313","node_587868","node_588504","node_592208","node_595643","node_599412","node_602919","node_606236","node_613324","node_613798","node_623458","node_720027","node_721384","node_721496","node_747889","node_764277","node_764302","node_1143179","node_1166894","node_1206604","node_1359222","node_1379395","node_1637668","node_1644381","node_1644961","node_384149","node_384162","node_384233","node_384262","node_384450","node_384582","node_384617","node_384686","node_384884","node_385093","node_385139","node_385165","node_391078","node_399969","node_400020","node_400086","node_400243","node_400469","node_400477","node_400492","node_400493","node_405128","node_405162","node_405521","node_405588","node_405957","node_406033","node_406103","node_406296","node_406477","node_407855","node_408690","node_409835","England/PHEC-YY8UB4H/2022|OX927065.1|2022-12-12","node_416334","node_419202","node_419318","node_419377","node_419707","node_419842","node_420840","node_422860","node_461172","node_495009","node_501499","node_501514","node_558382","node_563128","node_565520","node_572225","node_573987","node_575337","node_585361","node_589247","node_589783","node_595729","node_597108","node_597574","node_598221","node_604303","node_610304","node_610502","node_614095","node_619953","node_620845","node_623371","node_623666","node_623949","node_625276","node_625658","node_626023","node_626228","node_626964","node_645710","node_679593","node_679726","node_680517","node_714355","node_720039","node_721386","node_721466","Scotland/CLIMB-CM7YR4SZ/2024|OZ091091.1|2024-05-10","node_721627","node_721651","node_721702","node_721934","node_721937","node_722055","node_722205","node_722472","node_722557","node_722824","node_747295","node_747361","node_747774","England/CLIMB-CM7YGIMR/2024|OZ001435.1|2024-01-08","node_1288282","node_281635","node_281638","node_297766","node_297865","node_297928","node_297952","node_331500","node_331716","node_331767","node_331849","node_332363","node_332576","node_381924","node_383005","node_383467","node_383537","node_383708","node_384032","node_384403","node_384425","node_384545","node_384675","node_384746","node_385158","node_385163","(XBB.2.3)","node_385272","node_400021","node_400031","node_406781","node_409126","node_409160","node_411490","node_411539","node_412045","node_412988","node_419237","node_419591","node_419820","node_419843","node_420068","node_420098","node_421122","node_421195","node_432003","node_497031","node_500767","node_500805","node_501421","node_501591","node_501692","node_558483","node_558724","node_561234","node_562351","node_588776","node_589473","node_589877","node_592135","node_593035","node_596457","node_597110","node_598626","node_603884","node_604666","node_609375","node_609396","node_610483","node_612459","node_623845","node_626408","node_633449","node_635630","node_651526","node_656561","node_665902","node_679603","node_717568","node_718016","node_718794","node_720040","node_721468","node_722056","node_722526","node_722763","node_722882","node_747321","node_747436","node_747697","node_747727","node_763629","node_1072089","node_1365989","node_281644","node_282209","node_297585","BHR/320106594/2023|OQ439550.1|2023-02-05","node_297931","node_331328","node_331358","node_331750","node_331775","node_331834","node_331859","node_332581","node_333400","node_333437","node_333502","node_333550","node_333613","node_333629","node_333650","(XBB.1.5)","node_333709","node_381925","node_382156","node_382272","node_382697","node_382845","node_382958","node_382977","node_383506","node_383514","node_384033","node_384072","node_384172","node_384243","node_384267","node_384507","node_384552","node_384562","node_384588","node_384649","node_384775","node_384831","node_384874","node_400367","node_407605","node_411272","node_411304","node_411725","node_411820","node_412271","node_412989","node_413369","node_414086","node_419338","node_419709","node_419810","node_419941","node_420101","node_420901","node_421150","node_421153","node_432007","node_500663","node_562587","node_568043","node_569192","node_569856","node_572227","node_587075","node_590311","node_598630","node_598633","node_598673","node_626046","node_629215","node_635729","node_637261","node_638091","node_638670","node_640110","node_640496","node_640773","node_641081","node_641707","node_641998","node_642476","node_642989","node_643216","node_643493","node_643925","node_645720","node_648773","node_657491","node_662225","node_674297","node_674322","node_678913","node_713360","node_720142","node_722563","node_722829","node_722839","node_747824","node_763680","node_281667","node_282307","node_297886","(XBB.1.9)","node_298003","node_331320","node_331365","node_331839","node_331919","node_332365","node_332582","node_332586","USA/VA-VTVAS3-GSC42950/2023|OR223788.1|2023-05-31","node_332744","node_333634","node_333735","node_334558","node_336060","node_336493","node_339339","node_339483","node_340569","node_341331","node_342611","node_346521","node_346601","node_349658","node_349892","node_351252","node_352136","node_353481","node_356094","node_356557","node_356795","node_357086","node_357332","node_357552","node_358007","node_358449","node_360076","node_360631","node_361307","node_361570","node_361907","node_361975","node_362340","node_363958","node_364709","node_380999","node_381161","node_382039","USA/NY-PRL-230606_81A09/2023|OR144049.1|2023-06-03","node_383085","node_383710","node_384456","node_384602","node_384661","node_384832","node_384887","node_400023","node_400257","USA/CA-CDPH-A3000000314778/2023|OQ879503.1|2023-03-13","node_407457","node_409176","node_409185","node_409186","node_411434","node_412818","node_412880","node_413114","node_413287","node_413388","node_413643","node_413738","node_414322","node_414804","node_415339","node_415480","node_415603","node_415766","node_419822","node_420156","node_420257","node_420347","node_420409","node_420463","node_420587","node_420749","node_421151","node_484598","node_505825","node_508082","node_561340","node_561785","node_562357","node_577418","node_588659","node_589588","node_590348","node_590482","node_595869","node_597408","node_597472","node_629770","node_630456","node_635908","node_638924","node_639857","node_640819","node_642402","node_643953","node_644116","node_644210","node_651166","node_655115","node_661330","node_664436","node_664927","node_665172","node_668221","node_669103","node_670609","node_681445","node_682491","node_682943","node_687151","node_688738","node_690507","node_693096","node_693727","node_694653","node_695063","node_695870","node_696161","node_697436","node_697739","node_701890","node_702367","node_703462","node_705813","node_707247","node_708857","node_710004","node_710054","node_712058","node_712843","node_713042","node_713301","node_713495","node_713990","node_714062","node_720043","node_722567","node_722765","(JN.1)","node_723002","node_747416","node_764076","node_281663","node_281693","node_282006","node_282023","node_282318","node_297955","node_298021","node_307920","node_331208","node_331390","node_331402","node_331417","node_332170","node_332457","node_332745","node_333164","node_333618","node_333631","node_333726","node_334013","node_334067","node_334113","node_334837","node_335642","node_336545","node_337235","node_337384","node_341514","node_344979","node_345437","node_345890","node_348216","node_349081","node_350407","node_351354","node_352894","node_354031","node_354832","node_357380","node_357530","node_358452","node_358511","node_360104","node_361341","node_362249","node_362411","node_363917","node_364425","node_364502","node_364967","node_365154","node_365534","node_365674","node_366485","node_366699","node_367849","node_369018","node_369270","node_369946","node_370218","node_370776","node_371397","node_371716","node_373201","node_374607","node_375529","node_375940","node_376112","node_376869","node_377264","node_377600","node_378284","node_379471","node_380183","node_380785","node_381278","node_382212","node_382587","node_384007","node_400247","node_400321","node_412760","node_413527","node_413750","node_413924","node_414186","node_414347","node_415397","node_419085","node_420113","node_420490","node_420494","node_501246","node_557534","node_570150","node_572327","node_597308","node_626256","node_629533","node_637310","node_640316","node_642703","node_643929","node_644568","node_644590","node_659626","node_674421","node_674438","node_674498","node_686265","node_686403","node_688868","node_689383","node_691778","node_693181","node_693450","node_694658","node_696163","node_700261","node_701699","node_702383","node_702804","node_702835","node_703197","node_703561","node_705550","node_706125","node_706903","node_707102","node_708859","node_709423","node_709625","node_709911","node_710108","node_718507","node_722484","node_723263","node_723519","node_723749","node_723975","node_724287","node_724430","node_724479","node_725168","node_725879","node_726119","node_727372","node_727446","node_727732","node_727877","node_730777","node_730837","node_731182","node_731723","node_731983","node_732025","node_734134","node_734202","node_734566","node_734692","node_739685","node_739929","node_740835","node_741036","node_741144","node_741433","node_743810","node_743976","node_744046","node_744831","node_744870","node_744889","node_744913","node_746014","node_746041","node_746771","node_746913","England/CLIMB-CM7YKNRT/2024|OZ016822.1|2024-02-05","node_764080","node_282425","node_283473","node_297589","node_298063","node_331210","node_331220","node_331392","node_331753","node_331928","node_332107","node_332133","node_332149","node_333052","node_333176","node_333213","node_333354","node_333595","node_333737","node_334568","node_336547","node_337236","node_339341","node_347085","node_352897","node_352929","node_354057","node_355249","node_357405","node_357584","node_358686","node_360263","node_362280","node_366560","node_366764","node_367871","node_367882","node_368762","node_371274","node_373918","node_374733","node_375315","node_376026","node_376446","node_377601","node_378296","node_379108","node_380064","node_381184","node_381215","node_382055","node_382226","node_382243","node_383735","node_383745","node_383782","node_383804","node_383851","node_383901","node_383918","node_386023","node_386067","node_386110","node_386209","node_387741","node_387947","node_389040","node_389157","node_390558","node_390705","node_390750","node_390832","node_391013","node_400248","node_409172","node_413528","node_413568","Germany/IMS-10116-CVDP-E81A1A17-452E-4F01-8675-2FD09B76BA71/2023|OY262549.1|2023-03-15","node_501324","node_506297","(BQ.1)","node_508230","node_562425","node_562713","node_568047","node_591272","node_608110","node_640142","node_644503","node_644532","node_644552","node_644553","node_657185","node_663219","node_667120","node_669745","node_673678","node_683875","England/QEUH-325E9861/2022|OX379305.1|2022-11-05","node_687198","node_687373","node_687428","node_688239","node_689656","node_690899","node_697440","node_698389","node_698560","node_701017","node_702704","node_702881","node_704642","node_705620","node_707176","node_707216","node_708068","node_708973","node_709359","node_709559","node_709734","node_711081","node_713308","node_718215","node_723004","node_723597","node_723758","node_723962","node_724199","node_724306","node_724443","England/CLIMB-CM7Y8Z3Q/2024|OZ084055.1|2024-03-25","node_724480","node_724492","node_724502","node_724543","node_725340","node_725381","node_726121","node_726685","node_727477","node_727905","node_727942","node_730793","node_730926","node_731192","node_734938","node_735251","node_735407","node_735485","node_739844","node_740841","node_740929","node_741240","node_743104","node_743141","node_743329","node_743889","node_743924","node_743967","node_744062","node_744753","node_744866","node_744890","node_745112","node_745213","node_746357","node_746500","node_764081","node_281708","node_282426","node_282709","node_283112","node_283164","node_283213","node_283232","node_298183","node_299042","node_299063","node_299145","node_299223","node_299363","node_299396","node_299451","node_299601","node_299635","node_299756","node_300526","node_301022","node_301151","node_301751","node_301818","node_301825","node_301939","node_306376","node_306724","node_307716","node_307828","node_307858","node_307962","node_331961","node_331991","node_332949","node_332991","node_332998","node_333006","node_333037","node_333177","node_333801","node_333975","node_337387","node_337413","node_342586","node_354164","node_355097","USA/DC-CDC-2-6896453/2023|OQ608382.1|2023-02-07","node_355331","node_357592","node_359135","node_365056","node_365964","node_367905","node_368006","node_368763","node_368768","node_371824","node_372055","node_373459","USA/VA-GBW-H20-318-3972/2023|PQ018399.1|2023-04-19","node_379029","node_379167","node_380038","node_380068","node_381003","node_382056","node_382073","node_382231","node_383713","node_383736","node_383740","node_383853","node_383998","node_384005","node_385314","node_385325","node_385891","node_386148","node_387661","node_387768","node_387796","node_388867","node_389042","node_389078","node_389167","node_390574","node_390682","node_390706","node_413523","node_416400","node_416459","node_417848","node_418875","node_420261","node_499755","node_508231","node_508919","node_510940","node_510974","node_512574","node_512784","node_513107","node_514721","node_515657","node_516480","node_551368","node_556125","node_556219","node_556353","node_556521","node_557350","node_561498","node_592874","node_594796","node_598637","node_629847","node_636139","node_636586","node_675536","node_686457","node_686603","node_687681","node_689853","node_689973","node_697502","node_700663","node_702705","node_702706","node_702742","node_708287","node_709374","node_709468","node_711377","node_711421","node_722904","node_724200","node_724230","node_724448","node_724493","England/CLIMB-CM7YE7RY/2024|OZ081404.1|2024-06-18","node_725392","node_726122","node_726345","node_726686","node_730807","node_732090","node_732460","Switzerland/ZH-UZH-IMV-3ba6c99b/2024|OZ014401.1|2024-01-21","node_734802","node_735285","node_735350","node_735530","node_735679","node_737478","node_738193","node_738498","node_739687","node_739695","node_739949","node_740930","node_741040","England/CLIMB-CM7YGZFC/2024|OZ111159.1|2024-05-20","node_741714","node_743333","node_743896","node_744147","node_744932","node_745414","node_745802","node_746262","node_281709","node_281799","node_282321","node_282620","node_282850","node_282911","node_283214","node_283318","node_283341","node_298075","node_298153","node_298184","node_298388","node_299061","node_299169","node_299534","node_299614","node_299762","node_299994","node_300068","node_300131","node_300444","node_300822","node_301195","node_301352","node_301420","node_301567","node_302474","node_306823","node_306988","node_307728","node_307963","node_308093","node_308460","node_308681","node_309892","node_310170","node_330611","node_330626","node_330732","node_330748","node_331394","node_332608","node_332759","node_333581","node_333826","node_337395","node_339019","(XBB.1.5.70)","node_343075","node_346360","node_347115","node_347120","node_347138","node_350414","node_351420","node_355342","node_357669","node_357738","node_366218","node_366881","node_367057","node_367906","node_367999","node_368772","node_372361","node_377317","node_378410","node_379120","node_379168","node_382063","node_383730","node_383920","node_385821","node_386441","node_387608","node_388125","node_388202","node_388779","node_388829","node_388914","node_388992","node_390658","node_413412","node_417057","node_417816","node_417836","node_418784","node_419011","node_509108","node_511843","node_512446","node_513577","node_514480","node_514687","node_517545","node_518041","node_554225","node_555326","node_556066","node_568072","node_568285","node_629584","node_636308","node_637490","USA/VA-CAV_VAS3N_00018363_01/2023|OQ843743.1|2023-02-22","node_655141","node_676759","node_677544","node_686543","node_686547","node_693239","node_693377","node_698739","node_701443","node_701557","node_702647","node_709426","node_711083","node_711309","node_711418","node_722992","node_723305","node_723978","node_724449","USA/PR-CVL-025225/2024|PQ145336.1|2024-06-22","node_725882","(JN.1.11.1)","node_727959","node_730797","node_732044","node_732092","node_732461","node_733168","node_733206","node_733218","node_733243","node_733308","node_733936","node_735289","node_735291","node_739160","node_739690","node_741043","node_741046","node_741940","node_742301","node_742489","node_742728","node_742893","node_743113","node_743862","node_744755","node_744840","node_744892","node_745373","node_764267","node_281775","node_281786","node_282599","node_282711","node_282881","node_282933","node_283167","node_283225","node_298176","node_298286","node_298489","node_298677","node_299047","node_299624","node_299772","node_300107","node_300332","node_300502","node_300586","node_300697","node_301316","England/CLIMB-CM7Y8AHM/2023|OY983604.1|2023-12-07","node_301862","node_301944","node_302113","node_302330","node_302381","node_302535","node_302637","node_306677","node_307243","node_308015","node_308067","node_308068","node_308682","node_308691","node_308753","node_308819","node_308840","node_308904","node_331059","node_331195","node_333926","node_344105","node_344614","node_344630","node_344740","node_344743","Denmark/DCGC-663025/2023|OY784226.1|2023-11-06","node_347091","node_350415","node_352976","node_357632","node_366247","node_368765","node_371826","node_381026","node_381045","node_382084","node_383785","node_385354","node_385356","node_385772","node_385807","node_385953","node_386232","node_386298","node_387842","node_387854","node_387976","node_413484","node_413505","node_416461","node_418609","node_509109","node_509830","node_512916","node_512982","node_514172","node_518042","node_519622","node_520231","node_521378","node_521490","node_521557","node_521837","node_521954","node_522453","node_523363","node_523457","node_525527","node_526231","node_528216","node_529181","node_529856","node_529906","node_531696","node_531816","node_532037","node_532981","node_535183","node_535488","node_535866","node_536282","node_537022","node_540621","node_541009","node_541297","node_542795","node_543883","node_545398","node_547178","node_547754","node_548099","node_548202","node_548475","node_548591","node_548851","node_549570","node_550496","node_550519","USA/UT-UPHL-230127491739/2023|OQ371817.2|2023-01-11","node_555191","node_556221","node_556272","node_556586","node_568141","node_568189","node_674507","node_676837","node_702398","node_702504","node_702938","node_722934","USA/UT-UPHL-240502986258/2024|PP788897.1|2024-04-18","node_725397","node_725887","node_726521","node_728031","node_728046","node_728092","node_728147","node_729302","(KP.3)","node_729355","node_732463","node_732622","node_733309","node_733344","node_736551","node_744758","node_744760","node_745418","node_282333","node_282390","node_282440","node_282575","node_282935","node_283110","node_283185","node_283243","node_298067","node_298219","node_298316","node_299010","node_299773","node_299915","node_300186","node_300617","node_301222","node_301518","node_302016","node_302120","node_302246","node_302295","node_302684","node_306181","node_306218","node_306275","node_306732","USA/VA-GBW-H20-310-4001/2023|PQ019424.1|2023-05-08","node_308236","node_308607","node_309364","node_309465","node_309599","node_309701","node_309728","node_309894","node_310851","node_330750","node_331994","node_332615","node_333804","node_333944","node_337416","node_343029","node_343081","node_344489","node_383930","node_383937","node_383950","node_385668","node_385749","node_385815","node_386114","node_386132","node_386138","node_387551","node_387818","node_388210","node_388237","node_388804","node_389170","node_389281","node_389485","node_389963","node_390069","node_413394","node_512362","node_516609","node_518462","node_518650","node_518879","node_519316","node_519895","node_520750","node_521495","node_521566","node_522838","node_523458","node_524154","node_524226","node_525436","node_525528","node_526380","node_526576","node_526734","node_527038","node_528553","node_529077","node_529729","node_530314","node_530569","node_530768","node_531416","node_532300","node_533005","node_533067","node_536334","node_537496","node_538682","USA/NJ-CDC-LC0988560/2023|OQ295125.1|2023-01-09","node_539370","node_539480","node_539838","node_540623","node_540791","node_541181","node_541535","node_543446","node_543755","node_545318","node_545441","node_545617","node_545867","node_546241","node_546353","node_547788","node_547869","node_548838","node_549169","node_549191","node_549609","node_549775","node_550061","node_551470","node_702988","node_725893","node_726532","node_726566","node_726694","node_728104","node_728505","node_728614","node_728654","node_728928","node_728951","USA/CA-GBW-GKISBBBB46071/2024|PQ028885.1|2024-05-31","node_728976","node_729308","node_729537","node_730343","node_730349","node_733312","node_733940","node_737539","node_738501","node_741048","node_742738","node_743038","node_743090","USA/CA-GBW-GKISBBBC85779/2024|PQ214022.1|2024-08-05","USA/NJ-CDC-LC1102484/2024|PP600483.1|2024-03-27","node_745524","node_764101","node_281713","node_281812","node_282492","node_282631","node_282944","node_282986","node_283064","(XBB.1.16)","node_283494","node_298331","node_298920","node_298971","node_301519","node_301598","node_301864","node_301977","USA/TX-CDC-QDX49125560/2023|OQ892915.1|2023-04-11","node_302645","node_302673","node_307837","node_308004","node_308245","node_308955","node_309093","node_309327","node_310173","node_331011","node_332050","node_332763","node_337417","node_339030","node_339044","node_343104","node_343842","node_343994","node_344020","node_344052","node_344107","node_344119","node_344485","node_344503","node_344603","node_385289","node_385397","node_386133","node_386316","node_386444","node_388099","node_390179","node_413395","node_413396","node_508497","node_514176","node_514183","node_514526","node_516621","node_516951","node_517465","node_517565","node_526738","node_529869","node_530570","node_530643","node_531426","node_536369","node_538252","node_538578","node_540142","node_543539","node_544954","node_545178","node_546798","node_547212","node_551484","node_676830","node_722955","node_726536","OZ118874.1|2024-06-21","node_728226","node_728294","node_728346","node_728564","node_728582","node_728604","node_728649","node_728718","node_728745","node_728869","node_729339","node_729353","node_729357","node_729704","node_729844","node_733325","node_741062","node_745531","node_745637","England/CLIMB-CM7YF5KH/2024|OZ086013.1|2024-06-09","node_764157","node_282335","node_282770","node_282832","node_283889","node_284252","node_284312","node_286818","node_286948","node_286987","node_287276","node_287843","node_293022","node_294327","node_294547","node_296634","node_297057","node_297148","node_298492","node_298711","node_298729","node_299654","node_302701","node_309896","node_310667","node_310853","node_337502","node_337691","node_338515","node_343092","node_343293","node_343306","node_343858","node_343904","node_343967","node_344204","node_385959","node_386120","node_387250","node_387504","node_388281","node_510215","node_510679","node_514507","node_516625","node_520971","node_523581","node_523867","node_523933","node_530618","node_537680","node_539271","node_542333","node_545218","node_551594","node_728352","node_728411","USA/PA-CDC-QDX98316556/2024|PP973269.1|2024-06-12","node_728685","node_729095","node_729100","node_729103","node_729123","node_729157","node_729200","node_729324","node_729705","node_729738","node_729783","node_729845","node_729869","node_729933","England/CLIMB-CM7YRUYU/2024|OZ080971.1|2024-04-29","node_730273","node_733063","USA/WA-UW-22052946623/2022|ON740805.1|2022-05-29","node_764233","node_281814","node_282893","node_283779","node_284006","node_284150","node_285962","node_286149","node_286484","node_286591","node_287083","node_287352","node_289442","node_292886","node_292936","node_293173","node_294548","node_294962","node_295560","node_297285","node_298335","node_298730","node_300209","node_302702","node_309501","node_310508","node_311005","node_311016","(EG.5.1)","node_311180","node_331183","node_337447","node_337567","node_337666","node_337980","node_338950","node_343093","node_343432","node_343868","node_344073","node_344237","node_344329","node_344496","node_377330","node_385406","node_386672","node_386712","node_387251","node_390205","node_390208","node_413402","node_417063","node_417088","node_517574","node_523645","node_523877","node_525449","node_526751","node_536381","node_538424","node_541548","node_541960","node_542335","node_728125","node_728350","node_728412","node_728431","node_728433","node_728461","USA/MN-MDH-40801/2024|PP455087.1|2024-02-19","USA/CA-CDPH-500140864/2024|PQ210617.1|2024-07-12","node_728474","England/CLIMB-CM7YJ9HW/2024|OZ086542.1|2024-04-03","node_729238","node_729241","node_729248","node_729463","node_729479","node_729713","node_729796","USA/CO-CDPHE-42068884/2024|PQ158328.1|2024-06-28","node_729878","node_729982","node_281869","node_281888","node_282730","node_285605","node_287386","node_287554","node_287687","node_287712","node_287715","node_287754","node_287755","node_292907","node_294038","node_294594","node_295570","node_296031","node_296186","node_305882","node_311194","node_311238","node_311540","node_311740","node_312378","node_318659","node_318875","node_318956","node_321705","node_329710","node_330073","node_330346","node_330381","node_338383","node_343109","node_343893","node_343929","node_344123","node_344654","node_344657","node_374272","node_385682","node_386458","node_386556","node_387371","node_388577","node_390209","node_417377","node_517940","node_523658","node_525551","node_530649","node_541552","USA/VT-24-007204-WGS-01/2024|PQ187691.1|2024-07-24","node_728413","node_728425","node_728426","node_728443","node_728445","node_728446","node_729274","node_729379","USA/CA-CDC-LC1114866/2024|PQ237233.1|2024-08-07","USA/IL-CDC-LC1105937/2024|PP938830.1|2024-06-10","node_729726","node_729873","England/CLIMB-CM7YFOTR/2024|OZ091071.1|2024-04-23","node_729999","node_730028","node_730048","node_730091","node_730095","node_730111","node_730120","node_730122","England/CLIMB-CM7Y8J9S/2024|OZ124749.1|2024-07-08","USA/VA-CDC-LC1112451/2024|PQ207364.1|2024-07-29","node_287483","node_287516","node_293224","node_294304","node_297294","node_305499","node_305773","node_305904","node_311741","node_312239","node_312348","node_312583","node_319040","node_322570","node_338393","node_343897","node_344338","node_387253","node_390340","node_417090","node_417226","node_417369","node_417449","node_541664","node_541677","node_541694","node_284643","node_287719","node_294599","node_295330","node_298782","England/CLIMB-CM7YE1FM/2023|OY774288.1|2023-10-12","node_303122","node_303757","node_306018","node_310905","node_310988","node_312461","node_312678","node_312758","node_319233","node_319463","node_320201","node_321776","node_323376","node_343170","USA/FL-CDC-LC1089842/2024|PP231414.1|2024-01-16","USA/VA-GBW-H20-306-4687/2023|PQ020051.1|2023-05-15","node_387254","node_390215","node_417193","node_417241","USA/CA-LACPHL-AY02985/2023|OR736681.1|2023-10-05","node_530672","node_541585","node_290946","node_293635","node_294068","node_303151","node_303276","node_311821","node_312381","node_312413","node_312759","node_312955","node_313182","node_313627","node_313667","node_313777","node_313985","node_314233","node_314695","node_314697","node_315004","node_315751","node_315806","node_315837","node_315933","node_316026","(HK.3)","node_316174","node_318365","node_320409","node_325722","node_326237","node_326633","node_330024","node_417363","node_551614","THA/MTM_01_1170/2023|PP577694.1|2023-11-29","node_289969","node_290947","node_292157","node_293236","node_302797","node_303128","node_303183","node_310992","node_313183","node_313630","node_313986","node_314562","node_314569","node_314696","node_315037","node_315052","node_315320","node_315520","node_315982","node_316710","node_317218","node_317417","node_317431","node_317439","node_317592","node_317610","node_318040","node_320602","node_321062","node_321314","node_324260","node_326096","node_326709","node_328845","node_328921","node_374287","node_417286","node_284951","node_285066","node_290830","node_294095","node_294224","node_302802","node_303312","node_313200","node_313308","node_313744","node_314566","node_315041","node_315433","node_315549","node_315930","node_316259","node_316947","node_317031","node_317150","node_317563","node_318041","node_321406","node_323969","node_325269","node_329993","node_387265","node_417295","node_291323","node_314952","node_315327","node_315368","node_322950","node_323970","node_324835","node_302806","node_302881","node_323972","node_316269","node_316483","node_317069","node_323983"};


        std::string refFileName = globalVm["refFile"].as< std::string >();
        std::ifstream reffin(refFileName);
        std::string ref;
        std::getline(reffin, ref);
        T->printMutationsNew(fout, nodeIds, ref);
    } else if (globalVm.count("refFile")){
        std::string refFileName = globalVm["refFile"].as< std::string >();
        std::ifstream reffin(refFileName);
        std::string ref;
        std::getline(reffin, ref);
        T->printMutationsNew(fout, ref);
    } else {
        T->printMutationsNew(fout);
    }

    auto substitutionsEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds substitutionsTime = substitutionsEnd - substitutionsStart;
    std::cout << "\nMutation extract execution time: "
                << substitutionsTime.count() << " nanoseconds\n";

    if(globalVm.count("output-file")) outputFile.close();
}

void printPath(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }
    int treeID = 0;
    if(globalVm.count("treeID")) treeID = std::stoi(globalVm["treeID"].as< std::string >());

    panmanUtils::TreeGroup tg = *TG;
    panmanUtils::Tree * T = &TG->trees[treeID];
    // T = &tg.trees[treeID];


    if(globalVm.count("output-file")) {
        std::string fileName = globalVm["output-file"].as< std::string >();
        outputFile.open("./info/" + fileName + ".mutations");
        buf = outputFile.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream fout (buf);


    auto substitutionsStart = std::chrono::high_resolution_clock::now();

    std::cout << T->root->identifier << std::endl;

    // T->printMutations(fout);
    T->printNodePaths(fout);

    auto substitutionsEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds substitutionsTime = substitutionsEnd - substitutionsStart;
    std::cout << "\nMutation extract execution time: "
                << substitutionsTime.count() << " nanoseconds\n";

    if(globalVm.count("output-file")) outputFile.close();
}

void index(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // indexing
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    // Get start and end coordinate
    int64_t startCoordinate = 0;
    int64_t endCoordinate = -1;
    if(!globalVm.count("start")) {
        std::cout << "Start Coordinate not provided, setting it to 0" << std::endl;
    } else {
        startCoordinate = globalVm["start"].as< int64_t >();
    }

    if(!globalVm.count("end")) {
        std::cout << "End Coordinate not provided, setting it to length of seqeunce - 1" << std::endl;
    } else {
        endCoordinate = globalVm["end"].as< int64_t >();
    }

    // get sequence
    std::string reference="";
    if(!globalVm.count("reference")) {
        std::cout << "Error: Reference not provided" << std::endl;
        return;
    } else {
        reference = globalVm["reference"].as< std::string >();
    }

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree * T = &tg.trees[i];
        if (T->allNodes.find(reference) == T->allNodes.end()) {
            std::cout << "Error: reference " << reference << " does not exist in PanMAN\n";
            exit(0);
        }
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".index");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        bool allIndex = globalVm["index"].as< bool >();

        T->extractPanMATIndex(fout, startCoordinate,endCoordinate, reference, allIndex);

        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nIndexing execution time: " << fastaTime.count() << " nanoseconds\n";
}

void printRoot(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::ofstream &outputFile, std::streambuf * buf) {
    // Print raw sequences to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree * T = &tg.trees[i];
        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + "_" + std::to_string(i) + ".fasta");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);


        T->printFASTA(fout, true, true);

        if(globalVm.count("output-file")) outputFile.close();
    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nFASTA execution time: " << fastaTime.count() << " nanoseconds\n";

}

void toUsher(panmanUtils::TreeGroup *TG, po::variables_map &globalVm) {
    // Print raw sequences to output file
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    auto fastaStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        panmanUtils::Tree * T = &tg.trees[i];
        std::string fileName;
        if(globalVm.count("output-file")) {
            fileName = globalVm["output-file"].as< std::string >();
        } else {
            std::cout << "Output File not provided" << std::endl;
            return;
        }
        std::string refName;
        if(globalVm.count("reference")) {
            refName = globalVm["reference"].as< std::string >();
        } else {
            std::cout << "Reference not provided" << std::endl;
            return;
        }

        panmanUtils::panmanToUsher(T, refName, fileName);

    }

    auto fastaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;
    std::cout << "\nUsher Conversion time: " << fastaTime.count() << " nanoseconds\n";

}

void parseAndExecute(int argc, char* argv[]) {

    // Setup boost::program_options
    setupOptionDescriptions();

    // Initial command line arguments consisting of input file types
    po::variables_map globalVm;
    po::store(po::command_line_parser(argc, argv).options(globalDesc)
              .positional(globalPositionArgumentDesc).allow_unregistered().run(), globalVm);
    po::notify(globalVm);

    int threads = 64;
    if (globalVm.count("threads")) threads = globalVm["threads"].as<std::int32_t>();
    tbb::task_scheduler_init init(threads);


    // If the data structure loaded into memory is a PanMAT, it is pointed to by T
    panmanUtils::Tree *T = nullptr;

    // If the data structure loaded into memory is a PanMAN, it is pointed to by TG
    panmanUtils::TreeGroup *TG = nullptr;

    if(globalVm.count("help")) {
        std::cout << globalDesc;
        return;
    } else if (globalVm.count("protobuf2capnp")) {
        protobuf2capnp(TG, globalVm);
    } else if(globalVm.count("input-panmat")) {
        // Load PanMAT file directly into memory

        std::string fileName = globalVm["input-panmat"].as< std::string >();
        std::ifstream inputFile(fileName, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        inPMATBuffer.push(boost::iostreams::lzma_decompressor());
        inPMATBuffer.push(inputFile);
        std::istream inputStream(&inPMATBuffer);

        T = new panmanUtils::Tree(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        std::vector<panmanUtils::Tree*> tg;
        tg.push_back(T);

        TG = new panmanUtils::TreeGroup(tg);

        inputFile.close();


        writePanMAN(globalVm, TG);

        std::filesystem::create_directory("./info");


    } else if(globalVm.count("input-panman")) {
        // Load PanMAN file directly into memory

        std::string fileName = globalVm["input-panman"].as< std::string >();
        std::ifstream inputFile(fileName);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        inPMATBuffer.push(boost::iostreams::lzma_decompressor());
        // inPMATBuffer.push(boost::iostreams::gzip_decompressor());
        inPMATBuffer.push(inputFile);
        std::istream inputStream(&inPMATBuffer);

        std::cout << "starting reading panman" << std::endl;
        TG = new panmanUtils::TreeGroup(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";
        inputFile.close();

        std::filesystem::create_directory("./info");

    } else if(globalVm.count("input-gfa")) {
        // Create PanMAT from GFA and Newick files

        std::string fileName = globalVm["input-gfa"].as< std::string >();
        if(!globalVm.count("input-newick")) {
            panmanUtils::printError("File containing newick string not provided!");
            return;
        }
        if(!globalVm.count("output-file")) {
            panmanUtils::printError("Output file not provided!");
            std::cout << globalDesc;
            return;
        }
        std::string newickFileName = globalVm["input-newick"].as< std::string >();

        std::cout << "Creating PanMAN from GFA and Newick" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new panmanUtils::Tree(inputStream, newickInputStream, panmanUtils::FILE_TYPE::GFA);

        std::vector<panmanUtils::Tree*> tg;
        tg.push_back(T);

        TG = new panmanUtils::TreeGroup(tg);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, TG);

    } else if(globalVm.count("input-pangraph")) {
        // Create PanMAT from PanGraph and Newick files

        std::string fileName = globalVm["input-pangraph"].as< std::string >();
        if(!globalVm.count("input-newick")) {
            panmanUtils::printError("File containing newick string not provided!");
            std::cout << globalDesc;
            return;
        }
        if(!globalVm.count("output-file")) {
            panmanUtils::printError("Output file not provided!");
            std::cout << globalDesc;
            return;
        }

        std::string newickFileName = globalVm["input-newick"].as< std::string >();
        std::string referenceSequence;
        if(globalVm.count("reference")) {
            referenceSequence = globalVm["reference"].as< std::string >();
        }

        std::cout << "Creating PanMAN from PanGraph and Newick" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        T = new panmanUtils::Tree(inputStream, newickInputStream,
                                  panmanUtils::FILE_TYPE::PANGRAPH, referenceSequence);

        std::vector<panmanUtils::Tree*> tg;
        tg.push_back(T);

        TG = new panmanUtils::TreeGroup(tg);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, TG);

    } else if(globalVm.count("input-msa")) {
        // Create PanMAT from MSA and Newick files

        std::string fileName = globalVm["input-msa"].as< std::string >();
        if(!globalVm.count("input-newick")) {
            panmanUtils::printError("File containing newick string not provided!");
            return;
        }

        if(!globalVm.count("output-file")) {
            panmanUtils::printError("Output file not provided!");
            std::cout << globalDesc;
            return;
        }

        bool optimize = false;
        if(globalVm.count("low-mem-mode")) {
            optimize = true;
        }

        std::string reference = "";
        if (globalVm.count("reference")) {
            reference = globalVm["reference"].as<std::string>();
        }

        std::string newickFileName = globalVm["input-newick"].as< std::string >();

        std::cout << "Creating PanMAN from MSA and Newick" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        if(!optimize) {
            T = new panmanUtils::Tree(inputStream, newickInputStream,
                                  panmanUtils::FILE_TYPE::MSA, reference);
        } else {
            T = new panmanUtils::Tree(inputStream, newickInputStream,
                                    panmanUtils::FILE_TYPE::MSA_OPTIMIZE, reference);
        }

        // checkFunction(T);

        std::vector<panmanUtils::Tree*> tg;
        tg.push_back(T);

        TG = new panmanUtils::TreeGroup(tg);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, TG);

    } else if (globalVm.count("create-network")) {
        std::ofstream outputFile;
        std::streambuf * buf;
        createNet(globalVm, outputFile, buf);
        return;
    } else {
        panmanUtils::printError("Incorrect Format");
        std::cout << globalDesc;
        return;
    }

    // If only one function needs to be performed on the loaded PanMAT/PanMAN, do not start the
    // command line utility.
    std::ofstream outputFile;
    std::streambuf * buf;

    if(globalVm.count("summary")) {
        summary(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("printTips")) {
        printTips(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("impute")) {
        impute(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("fasta")) {
        fasta(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("fasta-aligned")) {
        fastaAligned(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("subnetwork")) { // for PanMAT -> Old
        subnetwork(T, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("subnet")) {
        subnet(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("vcf")) {
        vcf(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("gfa")) {
        gfa(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("maf")) {
        maf(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("newick")) {
        newick(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("extended-newick")) {
        extendNewick(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("annotate")) {
        annotate(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("reroot")) {
        reroot(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("aa-mutations")) {
        aa(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("printMutations")) {
        printMut(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("printNodePaths")) {
        printPath(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("index")) {
        index(TG, globalVm, outputFile, buf);
        return;
    } else if(globalVm.count("printRoot")) {
        printRoot(TG, globalVm, outputFile, buf);
        return;
    }  else if(globalVm.count("toUsher")) {
        toUsher(TG, globalVm);
        return;
    // } else if(globalVm.count("fasta-fast")){
    //     fastaFast(TG, globalVm, outputFile, buf);
    //     return;
    } else {
        char** splitCommandArray;

        while(true) {
            std::cout << "> ";

            std::string command;
            std::getline (std::cin, command);
            stripStringInPlace(command);

            // Split command by spaces
            std::vector< std::string > splitCommand;
            panmanUtils::stringSplit(command, ' ', splitCommand);
            splitCommandArray = new char*[splitCommand.size()];
            for(size_t i = 0; i < splitCommand.size(); i++) {
                splitCommandArray[i] = new char[splitCommand[i].length() + 1];
                strcpy(splitCommandArray[i], splitCommand[i].c_str());
            }

            try{
                if(strcmp(splitCommandArray[0], "help") == 0) {
                    std::cout << globalDesc;
                } 
                else if (strcmp(splitCommandArray[0], "root") == 0) {
                    buf = std::cout.rdbuf();
                    std::ostream fout (buf);
                    TG->trees[0].printFASTA(fout, true, true);
                } else if(strcmp(splitCommandArray[0], "impute") == 0) {
                    po::variables_map imputeVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(imputeDesc)
                        .run(), imputeVm);

                   impute(TG, imputeVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "summary") == 0) {
                    po::variables_map summaryVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(summaryDesc)
                        .run(), summaryVm);

                    summary(TG, summaryVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "fasta") == 0) {
                    po::variables_map fastaVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(fastaDesc)
                        .run(), fastaVm);

                    fasta(TG, fastaVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "fasta-aligned") == 0) {
                    po::variables_map fastaAlignVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(fastaAlignDesc)
                        .run(), fastaAlignVm);

                    fastaAligned(TG, fastaAlignVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "subnet") == 0) {
                    po::variables_map subnetVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(subnetDesc)
                        .run(), subnetVm);

                    subnet(TG, subnetVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "vcf") == 0) {
                    po::variables_map vcfVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(vcfDesc)
                        .run(), vcfVm);
                    vcf(TG, vcfVm, outputFile, buf);

                } else if(strcmp(splitCommandArray[0], "gfa") == 0) {
                    po::variables_map gfaVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(gfaDesc)
                        .run(), gfaVm);
                    gfa(TG, gfaVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "maf") == 0) {
                    po::variables_map mafVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(mafDesc)
                        .run(), mafVm);
                    
                    maf(TG, mafVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "newick") == 0) {
                    po::variables_map newickVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(newickDesc)
                        .run(), newickVm);
                    newick(TG, newickVm, outputFile, buf);

                } else if(strcmp(splitCommandArray[0], "extended-newick") == 0) {
                    po::variables_map extendNewickVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(extendNewickDesc)
                        .positional(fastaPositionArgumentDesc).run(), extendNewickVm);
                    extendNewick(TG, extendNewickVm, outputFile, buf);
                } else if (strcmp(splitCommandArray[0], "annotate") == 0) {
                    po::variables_map annotateVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(annotateDesc)
                        .run(), annotateVm);
                    annotate(TG, annotateVm, outputFile, buf);

                } else if (strcmp(splitCommandArray[0], "reroot") == 0) {
                    po::variables_map rerootVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(rerootDesc)
                        .run(), rerootVm);
                    reroot(TG, rerootVm, outputFile, buf);

                } else if (strcmp(splitCommandArray[0], "aa-mutations") == 0) {
                    po::variables_map aaVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(aaDesc)
                        .run(), aaVm);
                    aa(TG, aaVm, outputFile, buf);

                // } else if(strcmp(splitCommandArray[0], "create-network") == 0) {
                //     po::variables_map createNetVm;
                //     po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                //         .options(createNetDesc)
                //         .run(), createNetVm);
                //     createNet(TG, createNetVm, outputFile, buf);

                } else if(strcmp(splitCommandArray[0], "printMutations") == 0) {
                    po::variables_map printMutVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(printMutDesc)
                        .run(), printMutVm);
                    printMut(TG, printMutVm, outputFile, buf);

                } else if(strcmp(splitCommandArray[0], "printNodes") == 0) {
                    po::variables_map printNodeVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(printPathDesc)
                        .run(), printNodeVm);
                    printPath(TG, printNodeVm, outputFile, buf);
                } else if (strcmp(splitCommandArray[0], "index") == 0) {
                    po::variables_map indexVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(indexDesc)
                        .run(), indexVm);
                    index(TG, indexVm, outputFile, buf);
                } else if(strcmp(splitCommandArray[0], "printRoot") == 0) {
                    po::variables_map printRootVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(printRootDesc)
                        .run(), printRootVm);
                    printRoot(TG, printRootVm, outputFile, buf);
                } else if (strcmp(splitCommandArray[0], "exit") == 0 || strcmp(splitCommandArray[0], "q") == 0) {
                    return;
                } else {
                    std::cout << "type exit or q to exit" << std::endl;
                }
            } catch (std::exception& e) {
                std::cout << e.what() << std::endl;
            }
        }
    }
}

void debuggingCode() {
    std::string sequenceName = "NZ_AP019856.1";

    std::ifstream fin("/home/AD.UCSD.EDU/swalia/data/ecoli/pangraph/ecoli_1000.json");
    Json::Value pangraphData;
    fin >> pangraphData;

    std::cout << "LOADED" << std::endl;

    std::vector< std::string > blocks;
    std::vector< int > blockNumbers;
    std::vector< bool > strands;

    std::vector< std::vector< std::pair< char, std::vector< char > > > > sequence;
    std::unordered_map< std::string, std::string > stringIdToConsensusSeq;
    std::unordered_map< std::string, std::vector< std::pair< int, int > > > stringIdToGaps;
    std::unordered_map< std::string, std::unordered_map< int, std::vector< std::pair< int, std::string > > > > substitutions;
    std::unordered_map< std::string, std::unordered_map< int, std::vector< std::tuple< int, int, std::string > > > > insertions;
    std::unordered_map< std::string, std::unordered_map< int, std::vector< std::pair< int, int > > > > deletions;

    for(int i = 0; i < pangraphData["blocks"].size(); i++) {
        std::string blockId = pangraphData["blocks"][(int)i]["id"].asString();
        std::string stringSequence = pangraphData["blocks"][(int)i]["sequence"].asString();
        std::transform(stringSequence.begin(), stringSequence.end(),stringSequence.begin(), ::toupper);
        stringIdToConsensusSeq[blockId] = stringSequence;
        std::vector< std::string > gapMemberNames = pangraphData["blocks"][(int)i]["gaps"].getMemberNames();
        for(auto member: gapMemberNames) {
            stringIdToGaps[blockId].push_back( std::make_pair( std::stoi(member), pangraphData["blocks"][(int)i]["gaps"][member].asInt() ) );
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["mutate"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["name"].asString();
            if(seqName != sequenceName) {
                continue;
            }
            size_t number = pangraphData["blocks"][(int)i]["mutate"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["mutate"][(int)j][1].size(); k++) {
                std::string mutationString = pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                substitutions[blockId][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["mutate"][(int)j][1][(int)k][0].asInt(), mutationString) );
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["insert"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["name"].asString();
            if(seqName != sequenceName) {
                continue;
            }
            size_t number = pangraphData["blocks"][(int)i]["insert"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["insert"][(int)j][1].size(); k++) {
                std::string mutationString = pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][1].asString();
                std::transform(mutationString.begin(), mutationString.end(),mutationString.begin(), ::toupper);
                insertions[blockId][number].push_back( std::make_tuple( pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][0].asInt(), pangraphData["blocks"][(int)i]["insert"][(int)j][1][(int)k][0][1].asInt(), mutationString ) );
            }
        }
        for(size_t j = 0; j < pangraphData["blocks"][(int)i]["delete"].size(); j++) {
            std::string seqName = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["name"].asString();
            if(seqName != sequenceName) {
                continue;
            }
            size_t number = pangraphData["blocks"][(int)i]["delete"][(int)j][0]["number"].asInt();

            for(size_t k = 0; k < pangraphData["blocks"][(int)i]["delete"][(int)j][1].size(); k++) {
                deletions[blockId][number].push_back( std::make_pair( pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][0].asInt(), pangraphData["blocks"][(int)i]["delete"][(int)j][1][(int)k][1].asInt() ) );
            }
        }
    }

    std::cout << "blocks and mutations loaded" << std::endl;

    int totLength = 0;

    for(int i = 0; i < pangraphData["paths"].size(); i++) {
        if(pangraphData["paths"][i]["name"].asString() == sequenceName) {
            std::cout << "FOUND" << std::endl;
            std::unordered_map< std::string, int > numbers;
            for(int j = 0; j < pangraphData["paths"][i]["blocks"].size(); j++) {
                blocks.push_back(pangraphData["paths"][i]["blocks"][(int)j]["id"].asString());
                numbers[pangraphData["paths"][i]["blocks"][(int)j]["id"].asString()]++;
                blockNumbers.push_back(numbers[pangraphData["paths"][i]["blocks"][(int)j]["id"]
                                               .asString()]);
                strands.push_back(pangraphData["paths"][i]["blocks"][(int)j]["strand"].asBool());

                totLength += stringIdToConsensusSeq[blocks[j]].length();
            }
            sequence.resize(blocks.size());
            for(int j = 0; j < blocks.size(); j++) {
                // if(!strands[j]) {
                //     std::cout << "REVERSE STRAND FOUND" << std::endl;
                // }

                std::string sequenceString = stringIdToConsensusSeq[blocks[j]];
                sequence[j].resize(sequenceString.length()+1, {'-',{}});
                for(int k = 0; k < sequenceString.size(); k++) {
                    sequence[j][k].first = sequenceString[k];
                }
                for(size_t k = 0; k < stringIdToGaps[blocks[j]].size(); k++) {
                    sequence[j][stringIdToGaps[blocks[j]][k].first].second.resize(stringIdToGaps[blocks[j]][k].second, '-');
                }
                for(const auto& v: substitutions[blocks[j]][blockNumbers[j]]) {
                    sequence[j][v.first-1].first = v.second[0];
                }
                for(const auto& v: insertions[blocks[j]][blockNumbers[j]]) {
                    totLength += std::get<2>(v).length();
                    for(size_t k = 0; k < std::get<2>(v).length(); k++) {
                        sequence[j][std::get<0>(v)].second[std::get<1>(v)+k] = std::get<2>(v)[k];
                    }
                }
                for(const auto& v: deletions[blocks[j]][blockNumbers[j]]) {
                    totLength -= v.second;
                    for(size_t k = v.first; k < v.first + v.second; k++) {
                        sequence[j][k-1].first = '-';
                    }
                }
            }
            std::string sequenceString;
            for(int j = 0; j < sequence.size(); j++) {
                for(int k = 0; k < sequence[j].size(); k++) {
                    for(int w = 0; w < sequence[j][k].second.size(); w++) {
                        if(sequence[j][k].second[w] != '-') {
                            sequenceString += sequence[j][k].second[w];
                        }
                    }
                    if(sequence[j][k].first != '-') {
                        sequenceString += sequence[j][k].first;
                    }
                }
            }
            std::cout << sequenceString.substr(0, 10) << std::endl;
            std::cout << sequenceString.length() << std::endl;

            std::cout << "TOTAL LENGTH COMPUTED: " << totLength << std::endl;

            break;
        }
    }
    fin.close();
}

int main(int argc, char* argv[]) {
    parseAndExecute(argc, argv);
}
