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

#define PMAT_VERSION "2.0-beta"
#define VCF_VERSION "4.2"


static const int SANKOFF_INF = 100000001;

typedef std::vector< // block id
            std::pair< 
                std::vector< std::pair< char, std::vector< char > > >, // vector - nuc id, char - char at nuc id, vector <char> - nuc at gap id
                std::vector< 
                    std::vector< std::pair< char, std::vector< char > > > 
                > 
            > 
        > sequence_t;
// Individual block
typedef std::vector< std::pair< char, std::vector< char > > > block_t;

typedef  std::vector< std::pair< bool, std::vector< bool > > > blockExists_t;
// Forward or reverse strand
typedef  std::vector< std::pair< bool, std::vector< bool > > > blockStrand_t;

namespace panmanUtils {

enum FILE_TYPE {
    PANMAT = 0,
    GFA = 1,
    PANGRAPH=2,
    MSA = 3,
    MSA_OPTIMIZE = 4,
    // FASTA = 5
};


};