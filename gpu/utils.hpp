#ifndef UTILS
#define UTILS
#include "kseq.h"
#include "zlib.h"
#include <unordered_map>
#include <string>
#include <algorithm>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>

namespace utility {
    struct util {
        int global_batch_size;
        int local_batch_size;
        int global_batch_number;
        int local_batch_number;
        int num_nodes;
        int num_tips;
        int max_order;
        std::string ref_seq;
        std::string ref_name;
        int msa_len;
        std::string consensus;
        std::string seq_file_name;
        int start_coordinate;
        util(int gbs, int lbs, int gbn, int lbn);
    };
};


void read_seqs(std::string fname, std::unordered_map<std::string, std::string>& seqs, int start_idx=-1, int length=-1, bool print_details=0);
void make_seqs_equal(std::string fname, std::unordered_map<std::string, std::string>& seqs, bool start, int length, bool print_details=0);
#endif