#include <iostream>
#include <chrono>
#include <filesystem>
#include <tbb/parallel_for_each.h>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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

// program option description for building/loading a PanMAT into memory
po::options_description globalDesc("panmanUtils Command Line Arguments");
po::positional_options_description globalPositionArgumentDesc;

// program option descriptions of various command line functions
po::options_description summaryDesc("Summary Command Line Arguments");
po::options_description useDesc("Use Command Line Arguments");
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
po::options_description rerootDesc("Reroot Command Line Arguments");
po::positional_options_description rerootArgumentDesc;
po::options_description substitutionsDesc("Substitutions Command Line Arguments");
po::positional_options_description substitutionsArgumentDesc;
po::options_description aaTranslationDesc("Amino Acid Translation Command Line Arguments");
po::positional_options_description aaTranslationArgumentDesc;
po::options_description segmentExtractDesc("Segment Extract Command Line Arguments");
po::positional_options_description segmentExtractArgumentDesc;
po::options_description GFAToFASTADesc("GFA to Fasta writer Command Line Arguments");
po::positional_options_description GFAToFASTAArgumentDesc;
po::options_description groupWriteDesc("Group MAT Writer Command Line Arguments");
po::positional_options_description groupWritePositionArgumentDesc;
po::options_description sequenceExtractDesc("Sequence Extract Command Line Arguments");
po::positional_options_description sequenceExtractPositionArgumentDesc;
po::options_description groupFastaDesc("Tree Group FASTA writer Command Line Arguments");
po::positional_options_description groupFastaPositionArgumentDesc;

void setupOptionDescriptions() {
    // Global option descriptions
    globalDesc.add_options()
        ("help,h", "Print help messages")
        ("input-panman,I", po::value< std::string >(), "Input PanMAN file path")
        ("input-pangraph,P", po::value< std::string >(), "Input PanGraph JSON file to build a PanMAN")
        ("input-gfa,G", po::value< std::string >(), "Input GFA file to build a PanMAN")
        ("input-msa,M", po::value< std::string >(), "Input MSA file (FASTA format) to build a PanMAN")
        ("input-newick,N", po::value< std::string >(), "Input tree topology as Newick string")

        // ("optimize", "currently UNSUPPORTED: whether given msa file should be optimized or not")
        
        ("summary,s", "Print PanMAN summary")
        ("newick,t", "Print newick string of all trees in a PanMAN")
        ("fasta,f", "Print tip/internal sequences (FASTA format)")
        ("fasta-aligned,m", "Print MSA of sequences for each PanMAT in a PanMAN (FASTA format)")
        ("subnet,b", "Extract subnet of given PanMAN to a new PanMAN file based on the list of nodes provided in the input-file")
        ("vcf,v", "Print variations of all sequences from any PanMAT in a PanMAN (VCF format)")
        ("gfa,g", "Convert any PanMAT in a PanMAN to a GFA file")
        ("maf,w", "Print m-WGA for each PanMAT in a PanMAN (MAF format)")
        ("annotate,a", "Annotate nodes of the input PanMAN based on the list provided in the input-file")
        ("reroot,r", "Reroot a PanMAT in a PanMAN based on the input sequence id (--reference)")
        ("aa-translation,v", "Extract amino acid translations in tsv file")
        ("reference,n", po::value< std::string >(), "Identifier of reference sequence for PanMAN construction (optional), VCF extract (required), or reroot (required)")
        ("start,s", po::value< std::string >(), "Start coordinate of protein translation")
        ("end,e", po::value< std::string >(), "End coordinate of protein translation")

        ("input-file", po::value< std::string >(), "Path to the input file, required for --subnet and --annotate")
        ("output-file,o", po::value< std::string >(), "Prefix of the output file name")
    
        ("tree-group", po::value< std::vector< std::string > >()->multitoken(), "File paths of PMATs to generate tree group")
        ("mutation-file", po::value< std::string >(), "File path of complex mutation file for tree group")
        ("panman-in", po::value< std::string >(), "Input file path for PanMAT Group")
    
    ;

    // Adding input file as positional argument (doesn't require the --input-file tag)
    globalPositionArgumentDesc.add("input-panman", -1);

    // Use option descriptions
    useDesc.add_options()
        ("help", "produce help message")
        ("index", po::value< size_t >()->required(), "PanMAT index")
    ;

    // FASTA option descriptions
    fastaDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
        ("aligned", "print in aligned format (MSA)")
        ("parallel", "Whether we should execute in parallel or not")
    ;

    // Adding output file as positional argument (doesn't require the --output-file tag)
    fastaPositionArgumentDesc.add("output-file,o", -1);

    // MAF option descriptions
    mafDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument (doesn't require the --output-file tag)
    mafPositionArgumentDesc.add("output-file,o", -1);

    // MAT Writer option descriptions
    writeDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument (doesn't require the --output-file tag)
    writePositionArgumentDesc.add("output-file,o", -1);

    // Subtree Extract option descriptions
    subtreeDesc.add_options()
        ("help", "produce help message")
        ("newick", po::value< bool >()->default_value(false), "just print newick string")
        ("input-file", po::value< std::string >(), "Input file name if reading node IDs from file")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
        ("node-ids", po::value< std::vector< std::string > >()->multitoken(), "Node \
IDs to extract")
    ;

    // Adding output file as positional argument
    subtreePositionArgumentDesc.add("output-file,o", -1);

    // Sequence Extract option descriptions
    sequenceExtractDesc.add_options()
        ("help", "produce help message")
        ("list", po::value< std::vector< std::string > >()->multitoken()->required(), "Sequence names")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument (doesn't require the --output-file tag)
    sequenceExtractPositionArgumentDesc.add("output-file,o", -1);

    // VCF Writer option descriptions
    vcfDesc.add_options()
        ("help", "produce help message")
        ("reference", po::value< std::string >()->required(), "Sequence ID of the reference \
sequence")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
        ("fasta-file", po::value< std::string >(), "FASTA file name if it should also be created \
            from VCF File. Mainly used to verify the correctness of VCF file")
    ;

    // Adding output file as positional argument
    vcfPositionArgumentDesc.add("output-file,o", -1);

    // MAT Annotate option descriptions
    annotateDesc.add_options()
        ("help", "produce help message")
        ("input-file", po::value< std::string >()->required(), "Name of the file containing \
annotation info")
    ;

    // Adding input file as positional argument
    annotatePositionArgumentDesc.add("input-file", -1);

    // Search by annotation option descriptions
    searchDesc.add_options()
        ("help", "produce help message")
        ("keywords", po::value< std::vector< std::string > >()->multitoken(), "list of keywords to \
search for")
    ;

    searchPositionArgumentDesc.add("keywords", -1);
    
    // Generate GFA option descriptions
    generateGFADesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    generateGFAArgumentDesc.add("output-file,o", -1);

    // GFA to FASTA option descriptions
    GFAToFASTADesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
        ("input-file", po::value< std::string >()->required(), "Input file name")
    ;

    GFAToFASTAArgumentDesc.add("output-file,o", -1);

    rerootDesc.add_options()
        ("help", "produce help message")
        ("sequence-name", po::value< std::string >()->required(), "Name of sequence to reroot to")
    ;

    rerootArgumentDesc.add("sequence-name", -1);

    substitutionsDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Name of the output file")
    ;

    substitutionsArgumentDesc.add("output-file,o", -1);

    aaTranslationDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(),
            "Name of output file to store tsv file")
        ("start", po::value< int64_t >()->required(), "Root coordinate to start transcription")
        ("end", po::value< int64_t >()->required(), "Root coordinate to end transcription")
    ;

    aaTranslationArgumentDesc.add("output-file,o", -1);

    segmentExtractDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(),
            "Name of output file to store tsv file")
        ("start", po::value< int64_t >()->required(), "Root coordinate to start extraction")
        ("end", po::value< int64_t >()->required(), "Root coordinate to end extraction")
    ;

    segmentExtractArgumentDesc.add("output-file,o", -1);

    // Tree Group FASTA option descriptions
    groupFastaDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    groupFastaPositionArgumentDesc.add("output-file,o", -1);

    // Group MAT Writer option descriptions
    groupWriteDesc.add_options()
        ("help", "produce help message")
        ("output-file,o", po::value< std::string >()->required(), "Output file name")
    ;

    // Adding output file as positional argument
    groupWritePositionArgumentDesc.add("output-file,o", -1);
}

void writePanMAN(po::variables_map &globalVm, panmanUtils::Tree *T)
{
    std::string fileName = globalVm["output-file"].as< std::string >();
    std::filesystem::create_directory("./panman");

    std::ofstream outputFile("./panman/" + fileName + ".panman");
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

    auto writeStart = std::chrono::high_resolution_clock::now();

    outPMATBuffer.push(boost::iostreams::gzip_compressor());
    outPMATBuffer.push(outputFile);
    std::ostream outstream(&outPMATBuffer);
    T->writeToFile(outstream);
    boost::iostreams::close(outPMATBuffer);
    outputFile.close();

    auto writeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds writeTime = writeEnd - writeStart;
    std::cout << "\nTree Write execution time: " << writeTime.count()
        << " nanoseconds\n";

}

void parseAndExecute(int argc, char* argv[]) {

    // Setup boost::program_options
    setupOptionDescriptions();

    // Initial command line arguments consisting of input file types
    po::variables_map globalVm;
    po::store(po::command_line_parser(argc, argv).options(globalDesc)
        .positional(globalPositionArgumentDesc).allow_unregistered().run(), globalVm);
    po::notify(globalVm);

    // If the data structure loaded into memory is a PanMAT, it is pointed to by T
    panmanUtils::Tree *T = nullptr;

    // If the data structure loaded into memory is a PanMAN, it is pointed to by TG
    panmanUtils::TreeGroup *TG = nullptr;

    if(globalVm.count("help")) {
        std::cout << globalDesc;
        return;
    } else if(globalVm.count("input-panman")) {
        // Load PanMAT file directly into memory

        std::string fileName = globalVm["input-panman"].as< std::string >();
        std::ifstream inputFile(fileName);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        inPMATBuffer.push(boost::iostreams::gzip_decompressor());
        inPMATBuffer.push(inputFile);
        std::istream inputStream(&inPMATBuffer);

        T = new panmanUtils::Tree(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        inputFile.close();    

        std::filesystem::create_directory("./info");
        

    } else if(globalVm.count("panman-in")) {
        // Load PanMAN file directly into memory

        std::string fileName = globalVm["panman-in"].as< std::string >();
        std::ifstream inputFile(fileName);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        inPMATBuffer.push(boost::iostreams::gzip_decompressor());
        inPMATBuffer.push(inputFile);
        std::istream inputStream(&inPMATBuffer);

        TG = new panmanUtils::TreeGroup(inputStream);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;

        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";
        inputFile.close();
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

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, T);

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

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, T);

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
        if(globalVm.count("optimize")) {
            optimize = true;
        }

        std::string newickFileName = globalVm["input-newick"].as< std::string >();

        std::cout << "Creating PanMAN from MSA and Newick" << std::endl;

        std::ifstream inputStream(fileName);
        std::ifstream newickInputStream(newickFileName);

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        if(!optimize) {
            T = new panmanUtils::Tree(inputStream, newickInputStream,
                panmanUtils::FILE_TYPE::MSA);
        } else {
            T = new panmanUtils::Tree(inputStream, newickInputStream,
                panmanUtils::FILE_TYPE::MSA_OPTIMIZE);
        }

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        newickInputStream.close();
        inputStream.close();

        writePanMAN(globalVm, T);

    } else if(globalVm.count("tree-group")) {
        // Create PanMAN from list of PanMAT files and a complex mutation file listing the complex
        // mutations relating these PanMATs

        std::vector< std::string > fileNames;

        std::string mutationFileName;
        if(!globalVm.count("mutation-file")) {
            panmanUtils::printError("File containing complex mutations not provided!");
            return;
        }

        fileNames = globalVm["tree-group"].as< std::vector< std::string > >();
        mutationFileName = globalVm["mutation-file"].as< std::string >();

        std::ifstream mutationFile(mutationFileName);

        std::vector< std::ifstream > files;
        for(auto u: fileNames) {
            files.emplace_back(u);
        }

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();

        TG = new panmanUtils::TreeGroup(files, mutationFile);

        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Data load time: " << treeBuiltTime.count() << " nanoseconds \n";

        mutationFile.close();
        for(auto& u: files) {
            u.close();
        }
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
        // If command was summary, print the summary of the PanMAT
        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }

        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + ".summary");
            buf = outputFile.rdbuf();    
        } else {
            buf = std::cout.rdbuf();
        }
        auto summaryStart = std::chrono::high_resolution_clock::now();
        std::ostream fout (buf);
        T->printSummary(fout);
        auto summaryEnd = std::chrono::high_resolution_clock::now();

        if(globalVm.count("output-file")) outputFile.close();

        std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

        std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
        return;
    } else if(globalVm.count("fasta")) {
        // Print raw sequences to output file

        if(T == nullptr) {
            std::cout << "No PanMAT selected" << std::endl;
            return;
        }

        // std::string fileName = globalVm["output-file"].as< std::string >();

        // std::filesystem::create_directory("./fasta");
        // std::ofstream fout("./fasta/" + fileName + ".fasta");

        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + ".fasta");
            buf = outputFile.rdbuf();    
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        auto fastaStart = std::chrono::high_resolution_clock::now();
        
        T->printFASTA(fout, false);

        auto fastaEnd = std::chrono::high_resolution_clock::now();

        std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

        std::cout << "\nFASTA execution time: " << fastaTime.count()
            << " nanoseconds\n";
        if(globalVm.count("output-file")) outputFile.close();
        return;
    } else if(globalVm.count("fasta-aligned")) {
        // Print multiple sequence alignment to output file

        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }

        // std::string fileName = globalVm["output-file"].as< std::string >();

        // std::filesystem::create_directory("./fasta");
        // std::ofstream fout("./fasta/" + fileName + ".fasta");

        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + ".fasta");
            buf = outputFile.rdbuf();    
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        auto fastaStart = std::chrono::high_resolution_clock::now();
        
        T->printFASTA(fout, true);

        auto fastaEnd = std::chrono::high_resolution_clock::now();

        std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

        std::cout << "\nFASTA execution time: " << fastaTime.count()
            << " nanoseconds\n";

        if(globalVm.count("output-file")) outputFile.close();
        return;
    } else if(globalVm.count("subnet")) {
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
            std::cout << subtreeDesc;
        }

        if(nodeIds.size() == 0) {
            std::cout << "No node identifiers provided!" << std::endl;
        }

        std::string outputFileName = globalVm["output-file"].as< std::string >();
        std::filesystem::create_directory("./pmat");
        std::ofstream outputFile("./pmat/" + outputFileName + ".pmat");
        boost::iostreams::filtering_streambuf< boost::iostreams::output>
            outPMATBuffer;

        auto subtreeStart = std::chrono::high_resolution_clock::now();

        outPMATBuffer.push(boost::iostreams::gzip_compressor());
        outPMATBuffer.push(outputFile);
        std::ostream outstream(&outPMATBuffer);
        T->writeToFile(outstream, T->subtreeExtractParallel(nodeIds));
        boost::iostreams::close(outPMATBuffer);
        outputFile.close();

        auto subtreeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

        std::cout << "\nParallel Subtree Extract execution time: "
            << subtreeTime.count() << " nanoseconds\n";
        return;
    } else if(globalVm.count("vcf")) {
        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }

        if(!globalVm.count("reference")) {
            panmanUtils::printError("Reference ID not provided!");
            std::cout << globalDesc;
            return;
        }

        std::string reference = globalVm["reference"].as< std::string >();
        
        
        // std::string fileName = globalVm["output-file"].as< std::string >();

        // std::filesystem::create_directory("./vcf");
        // std::ofstream fout("./vcf/" + fileName + ".vcf");

        if(globalVm.count("output-file")) {
            std::string fileName = globalVm["output-file"].as< std::string >();
            outputFile.open("./info/" + fileName + ".vcf");
            buf = outputFile.rdbuf();    
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);

        auto vcfStart = std::chrono::high_resolution_clock::now();

        T->printVCFParallel(reference, fout);

        auto vcfEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;
        std::cout << "\nVCF execution time: " << vcfTime.count() << " nanoseconds\n";
        if(globalVm.count("output-file")) outputFile.close();
        return;
    } else if(globalVm.count("gfa")) {
        // If GFA is to be extracted from PanMAN

        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }

        // std::string fileName = globalVm["output-file"].as< std::string >();
        // std::filesystem::create_directory("./gfa");
        // std::ofstream fout("./gfa/"+fileName+".gfa");

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
        return;
    } else if(globalVm.count("maf")) {
        if(T == nullptr) {
            std::cout << "No PanMAN selected. Try groupFasta for FASTA of the whole"
                " PanMAN" << std::endl;
            return;
        }

        // std::string fileName = globalVm["output-file"].as< std::string >();

        // std::filesystem::create_directory("./maf");
        // std::ofstream fout("./maf/" + fileName + ".maf");

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
        return;
    } else if(globalVm.count("newick")) {
        // Print newick string of the PanMAT or PanMAN loaded into memory

        if(T) {
            std::cout << T->getNewickString(T->root) << std::endl;

            // std::cout << T->getNewickString(T->root) << std::endl;
        } else if(TG) {
            std::cout << "Printing newick string of each PanMAT in the PanMAN" << std::endl;
            int index = 0;
            for(auto& t: TG->trees) {
                std::cout << index++ << ": " << t.getNewickString(t.root) << std::endl;
            }
        } else {
            std::cout << "No PanMAN selected. Try groupFasta for FASTA of the whole PanMAN" << std::endl;
            return;
        }
    } else if (globalVm.count("anotate")) {
        // Annotate nodes of PanMAT

        
        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }

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
    } else if (globalVm.count("reroot")) {
        // Reroot the PanMAT to given sequence

        if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }
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
    } else if (globalVm.count("aa-mutations")) {
        // Extract amino acid translations in tsv file
        
       if(T == nullptr) {
            std::cout << "No PanMAN selected" << std::endl;
            return;
        }


        // std::string fileName = globalVm["output-file"].as< std::string >();
        // std::filesystem::create_directory("./aa_translations");
        // std::ofstream fout("./aa_translations/" + fileName + ".tsv");

        if(!globalVm.count("start") || !globalVm.count("start")) {
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
        return;
    }



    /*
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
            if(strcmp(splitCommandArray[0], "use") == 0) {
                // If command was use, select the PanMAT with the given index from the PanMAN
                po::variables_map useVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(useDesc)
                    .run(), useVm);

                if(useVm.count("help")) {
                    std::cout << useDesc;
                } else {
                    po::notify(useVm);
                    size_t treeIndex = useVm["index"].as< size_t >();
                    if(TG == nullptr) {
                        std::cout << "No PanMAN loaded" << std::endl;
                    } else {
                        if(TG->trees.size() > treeIndex) {
                            T = &TG->trees[treeIndex];
                            std::cout << "PanMAT loaded" << std::endl;
                        } else {
                            std::cout << "PanMAT with index " << treeIndex << " doesn't exist."
                                " There are only " << TG->trees.size() << " PanMATs." << std::endl;
                        }
                    }
                }
            } else if(strcmp(splitCommandArray[0], "summary") == 0) {
                // If command was summary, print the summary of the PanMAT
                if(T == nullptr) {
                    std::cout << "No PanMAT selected" << std::endl;
                    continue;
                }

                auto summaryStart = std::chrono::high_resolution_clock::now();
                T->printSummary();
                auto summaryEnd = std::chrono::high_resolution_clock::now();

                std::chrono::nanoseconds summaryTime = summaryEnd - summaryStart;

                std::cout << "\nSummary creation time: " << summaryTime.count() << " nanoseconds\n";
            } else if(strcmp(splitCommandArray[0], "fasta") == 0) {
                // Print the PanMAT sequences in FASTA format in output-file
                po::variables_map fastaVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(fastaDesc)
                    .positional(fastaPositionArgumentDesc).run(), fastaVm);

                if(fastaVm.count("help")) {
                    std::cout << fastaDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(fastaVm);
                    std::string fileName = fastaVm["output-file"].as< std::string >();

                    bool aligned = false;
                    bool parallel = 0;
                    if(fastaVm.count("aligned")) {
                        // Print MSA
                        aligned = true;
                    }
                    if(fastaVm.count("parallel")) {
                        // Compute sequences in parallel - often slower than default version since
                        // the algorithm is different
                        parallel = true;
                    }

                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto fastaStart = std::chrono::high_resolution_clock::now();
                    
                    if(parallel) {
                        T->printFASTAParallel(fout, aligned);
                    } else {
                        T->printFASTA(fout, aligned);
                    }

                    auto fastaEnd = std::chrono::high_resolution_clock::now();

                    std::chrono::nanoseconds fastaTime = fastaEnd - fastaStart;

                    std::cout << "\nFASTA execution time: " << fastaTime.count()
                        << " nanoseconds\n";

                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "maf") == 0) {
                // Print the MAF file for sequences in PanMAT

                po::variables_map mafVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(mafDesc).positional(mafPositionArgumentDesc).run(), mafVm);

                if(mafVm.count("help")) {
                    std::cout << mafDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected. Try groupFasta for FASTA of the whole"
                            " PanMAN" << std::endl;
                        continue;
                    }

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
            } else if(strcmp(splitCommandArray[0], "reroot") == 0) {
                // Reroot the PanMAT to given sequence
                po::variables_map rerootVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(rerootDesc).positional(rerootArgumentDesc).run(), rerootVm);

                if(rerootVm.count("help")) {
                    std::cout << rerootDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }


                    po::notify(rerootVm);

                    std::string sequenceName = rerootVm["sequence-name"].as< std::string >();

                    auto rerootStart = std::chrono::high_resolution_clock::now();

                    T->reroot(sequenceName);

                    auto rerootEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds rerootTime = rerootEnd - rerootStart;
                    std::cout << "\nReroot execution time: " << rerootTime.count()
                        << " nanoseconds\n";
                }
            } else if(strcmp(splitCommandArray[0], "substitutions") == 0) {
                // Print list of substitutions in sequence.
                po::variables_map substitutionsVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(substitutionsDesc).positional(substitutionsArgumentDesc).run(),
                    substitutionsVm);
                if(substitutionsVm.count("help")) {
                    std::cout << substitutionsDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(substitutionsVm);

                    std::string fileName = substitutionsVm["output-file"].as< std::string >();
                    
                    std::filesystem::create_directory("./mutations");
                    std::ofstream fout("./mutations/" + fileName);


                    auto substitutionsStart = std::chrono::high_resolution_clock::now();

                    T->printMutations(fout);

                    auto substitutionsEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds substitutionsTime = substitutionsEnd - substitutionsStart;
                    std::cout << "\nMutation extract execution time: "
                        << substitutionsTime.count() << " nanoseconds\n";

                    fout.close();
                }

            } else if(strcmp(splitCommandArray[0], "aaTranslation") == 0) {
                // Extract amino acid translations in tsv file
                po::variables_map aaTranslationVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(aaTranslationDesc)
                    .positional(aaTranslationArgumentDesc).run(), aaTranslationVm);

                if(aaTranslationVm.count("help")) {
                    std::cout << aaTranslationDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(aaTranslationVm);
                    std::string fileName = aaTranslationVm["output-file"].as< std::string >();
                    int64_t startCoordinate = aaTranslationVm["start"].as< int64_t >();
                    int64_t endCoordinate = aaTranslationVm["end"].as< int64_t >();

                    std::filesystem::create_directory("./aa_translations");
                    std::ofstream fout("./aa_translations/" + fileName + ".tsv");

                    auto aaStart = std::chrono::high_resolution_clock::now();

                    T->extractAminoAcidTranslations(fout, startCoordinate, endCoordinate);

                    auto aaEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds aaTime = aaEnd - aaStart;
                    std::cout << "\nAmino Acid translate execution time: " << aaTime.count()
                            << " nanoseconds\n";
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "extract-segment") == 0) {
                // Extract PanMAT representing a segment
                po::variables_map segmentExtractVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(segmentExtractDesc)
                    .positional(segmentExtractArgumentDesc).run(), segmentExtractVm);

                if(segmentExtractVm.count("help")) {
                    std::cout << segmentExtractDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(segmentExtractVm);
                    std::string fileName = segmentExtractVm["output-file"].as< std::string >();
                    int64_t startCoordinate = segmentExtractVm["start"].as< int64_t >();
                    int64_t endCoordinate = segmentExtractVm["end"].as< int64_t >();

                    std::filesystem::create_directory("./pmat");
                    std::ofstream fout("./pmat/" + fileName + ".pmat");
                    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

                    auto segmentExtractStart = std::chrono::high_resolution_clock::now();
                    outPMATBuffer.push(boost::iostreams::gzip_compressor());
                    outPMATBuffer.push(fout);
                    std::ostream outstream(&outPMATBuffer);

                    T->extractPanMATSegment(outstream, startCoordinate, endCoordinate);

                    auto segmentExtractEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds segmentExtractTime = segmentExtractEnd - segmentExtractStart;
                    std::cout << "\nSegment Extract execution time: " << segmentExtractTime.count()
                            << " nanoseconds\n";

                    boost::iostreams::close(outPMATBuffer);
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "vcf") == 0) {
                // Print the VCF file for sequences in PanMAT
                po::variables_map vcfVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(vcfDesc).positional(vcfPositionArgumentDesc).run(), vcfVm);

                if(vcfVm.count("help")) {
                    std::cout << vcfDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(vcfVm);

                    std::string reference = vcfVm["reference"].as< std::string >();
                    std::string fileName = vcfVm["output-file"].as< std::string >();

                    std::filesystem::create_directory("./vcf");
                    std::ofstream fout("./vcf/" + fileName + ".vcf");

                    auto vcfStart = std::chrono::high_resolution_clock::now();

                    T->printVCFParallel(reference, fout);

                    auto vcfEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds vcfTime = vcfEnd - vcfStart;
                    std::cout << "\nVCF execution time: " << vcfTime.count() << " nanoseconds\n";
                    fout.close();

                    // Generate FASTA file from VCF file to verify that the VCF file is correct
                    if(vcfVm.count("fasta-file")) {
                        std::cout << "Generating FASTA File" << std::endl;
                        std::ifstream fin("./vcf/" + fileName + ".vcf");
                        std::string fastaFileName = vcfVm["fasta-file"].as< std::string >();
                        std::filesystem::create_directory("./fasta");
                        fout.open("./fasta/" + fastaFileName + ".fasta");
                        T->vcfToFASTA(fin, fout);
                        fout.close();
                        fin.close();
                    }
                }
            } else if(strcmp(splitCommandArray[0], "subnet") == 0) {
                // Extract the subnet consisting of given node IDs from PanMAT
                po::variables_map subtreeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(subtreeDesc).positional(subtreePositionArgumentDesc).run(), subtreeVm);

                if(subtreeVm.count("help")) {
                    std::cout << subtreeDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(subtreeVm);
                    std::string outputFileName = subtreeVm["output-file"].as< std::string >();

                    // List of node identifiers that need to be extracted from the tree
                    std::vector< std::string > nodeIds;
                    std::string nodeId;

                    if(subtreeVm.count("input-file")) {
                        std::string inputFileName = subtreeVm["input-file"].as< std::string >();
                        std::ifstream fin(inputFileName);
                        while(fin >> nodeId) {
                            nodeIds.push_back(nodeId);
                        }
                        fin.close();
                    } else if(subtreeVm.count("node-ids")) {
                        nodeIds = subtreeVm["node-ids"].as< std::vector< std::string > >();
                    } else {
                        panmanUtils::printError("No source of node ids provided");
                        std::cout << subtreeDesc;
                    }

                    if(nodeIds.size() == 0) {
                        std::cout << "No node identifiers selected!" << std::endl;
                    }

                    if(subtreeVm["newick"].as< bool >()) {
                        // Only write newick string

                        std::filesystem::create_directory("./pmat");
                        std::ofstream fout("./newick/" + outputFileName + ".newick");

                        auto subtreeStart = std::chrono::high_resolution_clock::now();
                        fout << T->getNewickString(T->subtreeExtractParallel(nodeIds));

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: "
                            << subtreeTime.count() << " nanoseconds\n";
                        fout.close();
                    } else {
                        std::filesystem::create_directory("./pmat");
                        std::ofstream outputFile("./pmat/" + outputFileName + ".pmat");
                        boost::iostreams::filtering_streambuf< boost::iostreams::output>
                            outPMATBuffer;

                        auto subtreeStart = std::chrono::high_resolution_clock::now();

                        outPMATBuffer.push(boost::iostreams::gzip_compressor());
                        outPMATBuffer.push(outputFile);
                        std::ostream outstream(&outPMATBuffer);
                        T->writeToFile(outstream, T->subtreeExtractParallel(nodeIds));
                        boost::iostreams::close(outPMATBuffer);
                        outputFile.close();

                        auto subtreeEnd = std::chrono::high_resolution_clock::now();
                        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;

                        std::cout << "\nParallel Subtree Extract execution time: "
                            << subtreeTime.count() << " nanoseconds\n";
                    }
                }
            } else if(strcmp(splitCommandArray[0], "write") == 0) {
                // Write PanMAT to a file

                po::variables_map writeVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(writeDesc).positional(writePositionArgumentDesc).run(), writeVm);

                if(writeVm.count("help")) {
                    std::cout << writeDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(writeVm);
                    std::string fileName = writeVm["output-file"].as< std::string >();
                    std::filesystem::create_directory("./pmat");

                    std::ofstream outputFile("./pmat/" + fileName + ".pmat");
                    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

                    auto writeStart = std::chrono::high_resolution_clock::now();

                    outPMATBuffer.push(boost::iostreams::gzip_compressor());
                    outPMATBuffer.push(outputFile);
                    std::ostream outstream(&outPMATBuffer);
                    T->writeToFile(outstream);
                    boost::iostreams::close(outPMATBuffer);
                    outputFile.close();

                    auto writeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds writeTime = writeEnd - writeStart;
                    std::cout << "\nTree Write execution time: " << writeTime.count()
                        << " nanoseconds\n";
                }
            } else if(strcmp(splitCommandArray[0], "annotate") == 0) {
                // Annotate nodes of PanMAT

                po::variables_map annotateVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(annotateDesc).positional(annotatePositionArgumentDesc)
                    .run(), annotateVm);

                if(annotateVm.count("help")) {
                    std::cout << annotateDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(annotateVm);
                    std::string fileName = annotateVm["input-file"].as< std::string >();
                    std::ifstream fin(fileName);
                    auto annotateStart = std::chrono::high_resolution_clock::now();

                    T->annotate(fin);

                    auto annotateEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds annotateTime = annotateEnd - annotateStart;
                    std::cout << "Annotate time: " << annotateTime.count() << " nanoseconds\n";
                }
            } else if(strcmp(splitCommandArray[0], "search") == 0) {
                // Search nodes of PanMAT by annotations

                po::variables_map searchVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(searchDesc).positional(searchPositionArgumentDesc).run(), searchVm);

                if(searchVm.count("help")) {
                    std::cout << searchDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    std::vector< std::string > annotationVector = searchVm["keywords"]
                        .as< std::vector< std::string > >();

                    for(auto word: annotationVector) {
                        stripStringInPlace(word);
                        std::cout << word << ": ";
                        auto result = T->searchByAnnotation(word);
                        for(auto r: result) {
                            std::cout << r << "; ";
                        }
                        std::cout << std::endl << std::endl;
                    }
                }
            } else if(strcmp(splitCommandArray[0], "sequences") == 0) {
                // Print FASTA file consisting of given sequences from PanMAT

                po::variables_map sequenceExtractVm;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(sequenceExtractDesc).positional(sequenceExtractPositionArgumentDesc)
                    .run(), sequenceExtractVm);

                if(sequenceExtractVm.count("help")) {
                    std::cout << sequenceExtractDesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(sequenceExtractVm);

                    std::string fileName = sequenceExtractVm["output-file"].as< std::string >();
                    std::vector< std::string > sequenceNames = sequenceExtractVm["list"]
                        .as< std::vector< std::string > >();

                    std::filesystem::create_directory("./fasta");
                    std::ofstream fout("./fasta/" + fileName + ".fasta");

                    auto sequenceExtractStart = std::chrono::high_resolution_clock::now();
                    std::mutex fastaMutex;

                    tbb::parallel_for_each(sequenceNames, [&](auto sequenceName) {
                        std::string sequenceString = T->getStringFromReference(sequenceName, false);
                        if(sequenceString.find("Error:") == std::string::npos) {
                            fastaMutex.lock();
                            fout << ">" << sequenceName << '\n';
                            for(size_t i = 0; i < sequenceString.length(); i+=70) {
                                fout << sequenceString.substr(i,std::min(sequenceString.length()-i,
                                    (size_t)70)) << '\n';
                            }
                            fastaMutex.unlock();
                        } else {
                            std::cout << sequenceString << std::endl;
                        }
                    });

                    auto sequenceExtractEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds sequenceExtractTime = sequenceExtractEnd
                        - sequenceExtractStart;
                    std::cout << "\nSequence Extract execution time: "
                        << sequenceExtractTime.count() << " nanoseconds\n";
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "newick") == 0) {
                // Print newick string of the PanMAT or PanMAN loaded into memory

                if(T) {
                    std::cout << T->getNewickString(T->root) << std::endl;

                    // std::cout << T->getNewickString(T->root) << std::endl;
                } else if(TG) {
                    std::cout << "Printing newick string of each PanMAT in the PanMAN" << std::endl;
                    int index = 0;
                    for(auto& t: TG->trees) {
                        std::cout << index++ << ": " << t.getNewickString(t.root) << std::endl;
                    }
                }
            } else if(strcmp(splitCommandArray[0], "gfa") == 0) {
                // If command was genGFA

                po::variables_map generateGFAVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(generateGFADesc).positional(generateGFAArgumentDesc).run(),
                    generateGFAVM);
                if(generateGFAVM.count("help")) {
                    std::cout << generateGFADesc;
                } else {
                    if(T == nullptr) {
                        std::cout << "No PanMAT selected" << std::endl;
                        continue;
                    }

                    po::notify(generateGFAVM);
                    std::string fileName = generateGFAVM["output-file"].as< std::string >();
                    std::filesystem::create_directory("./gfa");
                    std::ofstream fout("./gfa/"+fileName+".gfa");

                    auto generateVGStart = std::chrono::high_resolution_clock::now();

                    T->convertToGFA(fout);

                    auto generateVGEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds generateVGTime = generateVGEnd - generateVGStart;

                    std::cout << "GFA generation time: " << generateVGTime.count()
                        << " nanoseconds\n";
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "gfa-fasta") == 0) {
                // If command was gfa-fasta

                po::variables_map GFAToFASTAVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(GFAToFASTADesc).positional(GFAToFASTAArgumentDesc).run(),
                    GFAToFASTAVM);
                if(GFAToFASTAVM.count("help")) {
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

                    std::cout << "FASTA generation time: " << generateVGTime.count()
                        << " nanoseconds\n";

                    fin.close();
                    fout.close();
                }
            } else if(strcmp(splitCommandArray[0], "groupFasta") == 0) {
                // Print FASTA of PanMAN

                po::variables_map groupFastaVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(groupFastaDesc).positional(groupFastaPositionArgumentDesc).run(),
                    groupFastaVM);
                if(groupFastaVM.count("help")) {
                    std::cout << groupFastaDesc;
                }
                std::filesystem::create_directory("./fasta");

                std::string outputFileName = groupFastaVM["output-file"].as< std::string >();
                std::ofstream fout("./fasta/" + outputFileName + ".fasta");

                auto groupFastaStart = std::chrono::high_resolution_clock::now();

                TG->printFASTA(fout);

                auto groupFastaEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds groupFASTATime = groupFastaEnd - groupFastaStart;
                std::cout << "Group FASTA write time: " << groupFASTATime.count()
                    << " nanoseconds\n";
                fout.close();

            } else if(strcmp(splitCommandArray[0], "groupWrite") == 0) {
                // Write PanMAN to file

                po::variables_map groupWriteVM;
                po::store(po::command_line_parser((int)splitCommand.size(), splitCommandArray)
                    .options(groupWriteDesc).positional(groupWritePositionArgumentDesc).run(),
                    groupWriteVM);
                if(groupWriteVM.count("help")) {
                    std::cout << groupWriteDesc;
                } else {
                    std::filesystem::create_directory("./pman");
                    std::string outputFileName = groupWriteVM["output-file"].as< std::string >();
                    std::ofstream outputFile("./pman/" + outputFileName + ".pman");
                    boost::iostreams::filtering_streambuf< boost::iostreams::output> outPMATBuffer;

                    auto groupWriteStart = std::chrono::high_resolution_clock::now();
                    outPMATBuffer.push(boost::iostreams::gzip_compressor());
                    outPMATBuffer.push(outputFile);
                    std::ostream outstream(&outPMATBuffer);

                    TG->writeToFile(outstream);
                    
                    auto groupWriteEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds groupWriteTime = groupWriteEnd - groupWriteStart;

                    std::cout << "Group PanMAT write time: " << groupWriteTime.count()
                        << " nanoseconds\n";

                    boost::iostreams::close(outPMATBuffer);
                    outputFile.close();
                }
            } else if(strcmp(splitCommandArray[0], "cplx-mutations") == 0) {
                // Print Complex Mutations of PanMAN

                std::cout << "Complex Mutations:" << std::endl;
                TG->printComplexMutations();

            } else if(strcmp(splitCommandArray[0], "exit") == 0) {
                return;
            }
        } catch (std::exception& e) {
            std::cout << e.what() << std::endl;
        }
    }

    */
}

void debuggingCode(){
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
            if(seqName != sequenceName){
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
            if(seqName != sequenceName){
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
            if(seqName != sequenceName){
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
            for(int j = 0; j < blocks.size(); j++){
                // if(!strands[j]) {
                //     std::cout << "REVERSE STRAND FOUND" << std::endl;
                // }

                std::string sequenceString = stringIdToConsensusSeq[blocks[j]];
                sequence[j].resize(sequenceString.length()+1, {'-',{}});
                for(int k = 0; k < sequenceString.size(); k++){
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
    // std::ifstream fin("./gfa/ecoli_100_parallel.gfa");
    // std::ofstream fout("./fasta/ecoli_100_parallel.fasta");
    // std::unordered_map< size_t, std::string > nodeIdToSequence;
    // std::string line;
    // while(getline(fin, line)) {
    //     if(line[0] == 'S') {
    //         std::vector< std::string > words;
    //         std::string str;
    //         for(size_t i = 0; i < line.length(); i++) {
    //             if(line[i] == '\t') {
    //                 if(str.length())words.push_back(str);
    //                 str = "";
    //             } else {
    //                 str += line[i];
    //             }
    //         }
    //         if(str.length()) {
    //             words.push_back(str);
    //             str = "";
    //         }
    //         nodeIdToSequence[std::stoll(words[1])] = words[2];
    //     }
    // }

    // std::cout << nodeIdToSequence.size() << std::endl;
    // fin.clear();
    // fin.seekg(0);
    // while(getline(fin, line)) {
    //     if(line[0] == 'P'){
    //         std::vector< std::string > words;
    //         std::string str;
    //         for(size_t i = 0; i < line.length(); i++) {
    //             if(line[i] == '\t' || line[i] == ',' || line[i] == '*') {
    //                 if(str.length())words.push_back(str);
    //                 str = "";
    //             } else {
    //                 str += line[i];
    //             }
    //         }
    //         if(str.length()) {
    //             words.push_back(str);
    //             str = "";
    //         }
    //         fout << ">" << words[1] << "\n";
    //         for(size_t i = 2; i < words.size(); i++) {
    //             char strand = words[i][words[i].length()-1];
    //             words[i].pop_back();
    //             size_t blockId = std::stoll(words[i]);
    //             std::string seq = nodeIdToSequence[blockId];
    //             if(strand == '+') {
    //                 fout << seq;
    //             } else {
    //                 // if(words[1] == "NZ_CP006027.1") {
    //                 //     std::cout << seq << std::endl;
    //                 //     return;
    //                 // }
    //                 for(int j = seq.length()-1; j>=0; j--) {
    //                     fout << panmanUtils::getComplementCharacter(seq[j]);
    //                 }
    //                 // for(int j = seq.length()-1; j>=0; j--) {
    //                 //     fout << panmanUtils::getComplementCharacter(seq[j]);
    //                 // }
    //             }
    //         }
    //         fout << "\n";
    //     }
    // }
    // fout.close();
    // fin.close();
}

int main(int argc, char* argv[]) {
    tbb::task_scheduler_init init(32);
    // debuggingCode();
    parseAndExecute(argc, argv);
}
