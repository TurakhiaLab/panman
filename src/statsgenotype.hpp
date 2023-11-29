#ifndef STATS_GENOTYPE
#define STATS_GENOTYPE

#include <iostream>
#include <vector>
// #include "PangenomeMATV2.hpp"

using namespace std;

namespace statsgenotype {

static int getIndexFromNucleotide(char nucleotide) {
    switch(nucleotide) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '*':
            return 4;
        default:
            return -1;
    }
}

struct mutationMatrices {
    vector< vector<double> > submat;
    vector<double> insmat = {0};
    vector<double> delmat = {0};
    
    vector<double> total_submuts;
    double total_insmut = 0;
    double total_delmut = 0;
    
    mutationMatrices() {
        total_submuts.resize(4);
        submat.resize(4);
        for (size_t i = 0; i < 4; ++i) {
            submat[i].resize(4);
        }
    }
};

enum variationType {
    SNP = 1,
    INS = 2,
    DEL = 4
};

struct variationSite {
    // substitution only
    variationSite(
        size_t sid, char ref, size_t position, int variation_types, const string& nucs,
        const vector<string>& insertion_seqs, const vector<size_t>& deletion_lens, const string& errors
    ) {
        site_id = sid;
        ref_position = position;
        size_t offset = 0;
        site_info = (getIndexFromNucleotide(ref) << 3) + variation_types;
        read_errs.resize(5);

        if (variation_types && variationType::SNP) {
            for (auto i = 0; i < nucs.size(); ++i) {
                read_errs[getIndexFromNucleotide(nucs[i])].push_back(double(errors[i]) - 33.0);
            }
            offset += nucs.size();
        }

        if (variation_types && variationType::INS) {
            for (auto i = 0; i < insertion_seqs.size(); ++i) {
                insertions[insertion_seqs[i]].push_back(double(errors[i + offset] - 33.0));
            }
            offset += insertion_seqs.size();
        }

        if (variation_types && variationType::DEL) {
            for (auto i = 0; i < deletion_lens.size(); i++) {
                deletions[deletion_lens[i]].push_back(double(errors[i + offset] - 33.0));
            }
        }
    }

    size_t site_id;
    size_t ref_position;
    int8_t site_info; // 2 bit -> reference nuc, 3 bit -> varaition types
    
    // substitution
    vector< vector<double> > read_errs;

    // deletion
    // map<size_t, size_t> deletions;
    map<size_t, vector<double> > deletions;
    
    // insertion
    // map<string, size_t> insertions;
    map<string, vector<double> > insertions;
    
};



}


#endif //STATS_GENOTYPE
