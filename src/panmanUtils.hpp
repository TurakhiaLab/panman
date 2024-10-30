#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <json/json.h>
#include "panman.capnp.h"

#include "panman.hpp"

namespace panmanUtils {

static inline void printError(std::string e) {
    std::cout << "\033[1;31m" << "Error: " << "\033[0m" << e << std::endl;
}

// Get nucleotide character from 4-bit PanMAT code
char getNucleotideFromCode(int code);

// Get 4-bit PanMAT code from nucleotide character
char getCodeFromNucleotide(char nuc);

// Get complement Nucleotide character of given nucleotide character. Used to compute reverse
// complement
char getComplementCharacter(char nuc);

// Given a sequence and block presence/strand information, print the sequence in FASTA format
// where each line has length lineSize
void printSequenceLines(const sequence_t& sequence,
                        const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize,
                        bool aligned, std::ostream& fout, int offset = 0, bool debug = false);

void printSubsequenceLines(const sequence_t& sequence,\
                                     const blockExists_t& blockExists, blockStrand_t& blockStrand, size_t lineSize, 
                                     const std::tuple<int, int, int, int>& panMATStart, 
                                     const std::tuple<int, int, int, int>& panMATEnd, 
                                     bool aligned, std::ostream& fout, int offset=0, bool debug=false);

// Remove '-' character from sequence string
std::string stripGaps(std::string sequenceString);
std::string getDate();
std::string stripString(std::string s);

void stringSplit (std::string const& s, char delim, std::vector<std::string>& words);




// Represents input PanGraph information for PanMAT generation
class Pangraph {
  public:
    // Graph adjacency list
    size_t numNodes;
    std::vector< std::vector< size_t > > adj;
    std::unordered_map< std::string, std::vector< size_t > > intSequences;
    std::vector< size_t > topoSortedIntSequences;
    std::unordered_map< std::string, std::vector< std::string > > paths;

    // Added as a patch to incorporate strands
    std::unordered_map< std::string, std::vector< int > > strandPaths;

    // Represents the "number" parameter of each block
    std::unordered_map< std::string, std::vector< size_t > > blockNumbers;

    // Circular offset of given sequence. Zero if not circular
    std::unordered_map< std::string, int > circularSequences;

    // If a sequence was rotated, which block was it rotated by
    std::unordered_map< std::string, int > rotationIndexes;

    // Specifies whether sequence is inverted by the rotation algorithm or not
    std::unordered_map< std::string, bool > sequenceInverted;

    std::unordered_map< std::string, std::string > stringIdToConsensusSeq;

    // block identifier to list of gaps - < position, size > pairs
    std::unordered_map< std::string,
        std::vector< std::pair< size_t, size_t > > > stringIdToGaps;

    // Mapping from new integer ID of block to old string ID. Duplicated blocks have same string
    // ID
    std::unordered_map< size_t, std::string > intIdToStringId;

    // List of substitutions for each block
    tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
        tbb::concurrent_unordered_map< size_t,
        std::vector< std::pair< size_t, std::string > > > > > substitutions;
    // List of insertions for each block
    tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
        tbb::concurrent_unordered_map< size_t, std::vector< std::tuple< size_t, size_t,
        std::string > > > > > insertions;
    // List of deletions for each block
    tbb::concurrent_unordered_map< std::string, tbb::concurrent_unordered_map< std::string,
        tbb::concurrent_unordered_map< size_t,
        std::vector< std::pair< size_t, size_t > > > > > deletions;

    Pangraph(Json::Value& pangraphData);
    std::vector< size_t > getTopologicalSort();
    std::unordered_map< std::string,std::vector< int > >
    getAlignedSequences(const std::vector< size_t >& topoArray);
    // Patch to incorporate strands
    std::unordered_map< std::string,std::vector< int > >
    getAlignedStrandSequences(const std::vector< size_t >& topoArray);
};

// Represents input GFA information for PanMAT construction
class GfaGraph {
  public:
    size_t numNodes;
    // Graph adjacency list
    std::vector< std::vector< size_t > > adj;

    // Names of the sequences
    std::vector< std::string > pathIds;

    // The sequences themselves where each node has an integer ID
    std::vector< std::vector< size_t > > intSequences;

    std::vector< size_t > topoSortedIntSequences;

    // The strands of each block in the sequences
    std::vector< std::vector< int > > strandPaths;

    // Raw sequence corresponding to each node
    std::vector< std::string > intNodeToSequence;
    GfaGraph(const std::vector< std::string >& pathNames,
             const std::vector< std::vector< std::pair< std::string, bool > > >& sequences,
             std::map< std::string, std::string >& nodes);

    std::vector< size_t > getTopologicalSort();
    std::vector< std::vector< int64_t > > getAlignedSequences(const
            std::vector< size_t >& topoArray);
    std::vector< std::vector< int > > getAlignedStrandSequences(const
            std::vector< size_t >& topoArray);

    
};


};
