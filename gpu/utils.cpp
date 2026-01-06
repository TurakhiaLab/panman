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

void make_seqs_equal(std::string fname, std::unordered_map<std::string, std::string>& seqs, bool start, int length, bool print_details) {
    gzFile f_rd = gzopen(fname.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fname.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);
    size_t align_len = 0;
    while (kseq_read(kseq_rd) >= 0) {
        std::string seq = (kseq_rd->seq.s);
        
        seqs[kseq_rd->name.s] = seq;
        if (kseq_rd->seq.l < length) {
            fprintf(stderr, "ERROR: alignment length mismatch %s %d\n", kseq_rd->name.s, kseq_rd->seq.l);
            if (start)
                seqs[kseq_rd->name.s].insert(seqs[kseq_rd->name.s].begin(), '-');  
            else
                seqs[kseq_rd->name.s].push_back('-');
        }
        
        if (align_len == 0) align_len = seqs[kseq_rd->name.s].size();
        else if (align_len != seqs[kseq_rd->name.s].size()) {
            fprintf(stderr, "%s %s %d %d\n", kseq_rd->name.s, kseq_rd->seq.s, align_len, seqs[kseq_rd->name.s].size());
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
    FILE *file = fopen("seqs_new.aln", "w"); 
    for (auto &s: seqs)
        fprintf(file, ">%s\n%s\n",s.first.c_str(), s.second.c_str());
    fclose(file);

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
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        if (start_idx ==-1 && length==-1)
            seqs[kseq_rd->name.s] = seq;
        else {
            if (kseq_rd->seq.l < start_idx+length) {
                fprintf(stderr, "ERROR: alignment length mismatch %s %d\n", kseq_rd->name.s, kseq_rd->seq.l);
                seqs[kseq_rd->name.s] = seq.substr(start_idx, length-1);
                seqs[kseq_rd->name.s].push_back('-');   
            }
            else {
                seqs[kseq_rd->name.s] = seq.substr(start_idx, length);
            }

        }
        if (align_len == 0) align_len = kseq_rd->seq.l;
        else if (align_len != kseq_rd->seq.l) {
            fprintf(stderr, "%s %s\n", kseq_rd->name.s, kseq_rd->seq.s);
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
