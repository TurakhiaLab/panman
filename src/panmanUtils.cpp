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
po::options_description maskDesc("Masking Command Line Arguments");
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
    ("max-insertion-impute-distance,D", po::value< int64_t >(), "Maximum total branch length to move a node for insertion imputation [default 5]")
    ("mask", "Mask nucleotides and report imputation success")
    ("mask-p", po::value< std::int64_t >(), "p out of 100 chance of masking")
    ("create-network,K",po::value< std::vector<std::string>>(), "Create PanMAN with network of trees from single or multiple PanMAN files")

    // ("optimize", "currently UNSUPPORTED: whether given msa file should be optimized or not")

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
        ("max-insertion-impute-distance,D", po::value< int64_t >(), 
        "Maximum total branch length to move a node for insertion imputation [default 5]")
        ("output-file,o", po::value< std::string >(), "Output file name");
    
    maskDesc.add_options()
        ("input-file", po::value< std::string >(), "Input file name")
        ("mask-p", po::value< std::int64_t >(), "Mask nucleotides and report imputation success")
        ("max-insertion-impute-distance,D", po::value< int64_t >(), 
        "Maximum total branch length to move a node for insertion imputation [default 5]");

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

    int allowedDistance = 5;
    if(globalVm.count("max-insertion-impute-distance")) {
        allowedDistance = globalVm["max-insertion-impute-distance"].as< int64_t >();
    }

    auto imputeStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        tg.trees[i].imputeNs(allowedDistance);
    }

    auto imputeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds imputeTime = imputeEnd - imputeStart;
    std::cout << "\nImputation time: " << imputeTime.count() << " nanoseconds\n";

    writePanMAN(globalVm, &tg);
}

void mask(panmanUtils::TreeGroup *TG, po::variables_map &globalVm, std::streambuf * buf) {
    // If command was impute, create a new PanMAN with imputed sequences
    if(TG == nullptr) {
        std::cout << "No PanMAN selected" << std::endl;
        return;
    }

    panmanUtils::TreeGroup tg = *TG;

    int allowedDistance = 5;
    if(globalVm.count("max-insertion-impute-distance")) {
        allowedDistance = globalVm["max-insertion-impute-distance"].as< int64_t >();
    }
    int maskP = globalVm["mask-p"].as< int64_t >();

    auto maskStart = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < tg.trees.size(); i++) {
        tg.trees[i].testImputation(maskP, allowedDistance);
    }

    auto maskEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds maskTime = maskEnd - maskStart;
    std::cout << "\nMasking time: " << maskTime.count() << " nanoseconds\n";
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

    std::cout << T->root->identifier << std::endl;

    // T->printMutations(fout);
    T->printMutationsNew(fout);

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

    int threads = 32;
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
    } else if (globalVm.count("impute")) {
        impute(TG, globalVm, outputFile, buf);
        return;
    } else if (globalVm.count("mask")) {
        mask(TG, globalVm, buf);
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
                } else if(strcmp(splitCommandArray[0], "mask") == 0) {
                    po::variables_map maskVm;
                    po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                        .options(maskDesc)
                        .run(), maskVm);

                   mask(TG, maskVm, buf);
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
