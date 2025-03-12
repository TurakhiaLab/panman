#ifndef UTILS
#include "utils.hpp"
#endif
KSEQ_INIT2(, gzFile, gzread);

utility::util::util(int gbs, int lbs, int gbn, int lbn) {
    global_batch_size = gbs;
    local_batch_size = lbs;
    global_batch_number = gbn;
    local_batch_number = lbn;
    return;
}

void read_seqs(std::string fname, std::unordered_map<std::string, std::string>& seqs, int start_idx, int length, bool print_details) {
    gzFile f_rd = gzopen(fname.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fname.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);
    size_t align_len = 0;
    while (kseq_read(kseq_rd) >= 0) {
        std::string seq = (kseq_rd->seq.s);
        if (start_idx ==-1 && length==-1)
            seqs[kseq_rd->name.s] = seq;
        else
            seqs[kseq_rd->name.s] = seq.substr(start_idx, length);
        if (align_len == 0) align_len = kseq_rd->seq.l;
        else if (align_len != kseq_rd->seq.l) {
            fprintf(stderr, "ERROR: alignment length mismatch\n");
            exit(1);
        }
    }

    kseq_destroy(kseq_rd);
    gzclose(f_rd);

if (print_details) {
    printf("=== Sequence information ===\n");
    printf("Number : %d\n", seqs.size());
    printf("Alignment Length: %d\n",align_len);
    printf("============================\n");
}

    return;
}
