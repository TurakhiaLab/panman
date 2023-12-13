#ifndef STATS_GENOTYPE
#define STATS_GENOTYPE

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

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

static char getNucleotideFromIndex(int index) {
    switch(index) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return '*';
        default:
            return 'N';
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

static double phred_complement(double q) {
    double p = pow(10, (-q / 10));
    return -10 * log10(1 - p);
}

static double likelihood(
    int genotype_idx,
    const vector< vector<double> >& read_errs,
    const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions,
    int variation_type
) {
    vector<double> genotype_probs;
    vector<double> variants_probs;

    for (int i = 0; i < read_errs.size(); i++) {
        const auto& row = read_errs[i];
        if ((variation_type & variationType::SNP) && (genotype_idx == i)) {
            for (const auto& prob : row) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), row.begin(), row.end());
        }
    }

    int ins_i = 0;
    for (const auto& insertion : insertions) {
        if ((variation_type & variationType::INS) && (genotype_idx == ins_i)) {
            for (const auto& prob : insertion.second) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), insertion.second.begin(), insertion.second.end());
        }
        ins_i++;
    }
    
    int del_i = 0;
    for (const auto& deletion : deletions) {
        if ((variation_type & variationType::DEL) && (genotype_idx == del_i)) {
            for (const auto& prob : deletion.second) {
                genotype_probs.push_back(phred_complement(prob));
            }
        } else {
            variants_probs.insert(variants_probs.end(), deletion.second.begin(), deletion.second.end());
        }
        del_i++;
    }
    
    double genotype_prob = accumulate(genotype_probs.begin(), genotype_probs.end(), 0.0);
    double variant_prob = accumulate(variants_probs.begin(), variants_probs.end(), 0.0);
    
    return genotype_prob + variant_prob;

}


static vector<double> genotype_likelihoods(
    const vector< vector<double> >& read_errs, const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions, const int8_t& site_info
) {
    vector<double> likelihoods;
    likelihoods.resize(5);
    for (auto i = 0; i < 5; ++i) {
        likelihoods[i] = numeric_limits<double>::max();
    }

    auto variation_types = site_info & 7;
    auto ref_nuc = site_info >> 3;

    for (int i = 0; i < read_errs.size(); ++i) {
        if (read_errs[i].empty() && i != ref_nuc) { continue; }
        likelihoods[i] = likelihood(i, read_errs, deletions, insertions, variationType::SNP);
    }

    if (variation_types & variationType::INS) {
        for (auto i = 0; i < insertions.size(); i++) {
            likelihoods.push_back(likelihood(i, read_errs, deletions, insertions, variationType::INS));
        }
    }

    if (variation_types & variationType::DEL) {
        for (auto i = 0; i < deletions.size(); i++) {
            likelihoods.push_back(likelihood(i, read_errs, deletions, insertions, variationType::DEL));
        }
    }

    return likelihoods;
}

static vector<double> genotype_posteriors(
    const vector<double>& likelihoods, const map<string, vector<double> >& deletions,
    const map<string, vector<double> >& insertions, const int8_t& site_info, const mutationMatrices& mutmat
) {
    vector<double> posteriors;
    posteriors.resize(4);
    for (auto i = 0; i < 4; i++) {
        posteriors[i] = numeric_limits<double>::max();
    }

    auto ref_nuc = site_info >> 3;
    for (auto i = 0; i < 4; ++i) {
        if (likelihoods[i] != numeric_limits<double>::max()) {
            posteriors[i] = likelihoods[i] + mutmat.submat[ref_nuc][i];
        }
    }

    size_t insertion_idx = 0;
    for (const auto& insertion : insertions) {
        posteriors.push_back(likelihoods[5 + insertion_idx] + mutmat.insmat[insertion.first.size()]);
        insertion_idx++; 
    }

    size_t deletion_idx = 0;
    for (const auto& deletion : deletions) {
        posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] + mutmat.insmat[deletion.first.size()]);
        insertion_idx++;
    }

    double min_score = *min_element(posteriors.begin(), posteriors.end());
    for (int i = 0; i < posteriors.size(); ++i) {
        posteriors[i] -= min_score;
    }
    return posteriors;
}

struct variationSite {
    variationSite(
        size_t sid, char ref, size_t position, int variation_types, const string& nucs,
        const vector<string>& insertion_seqs, const vector<string>& deletion_seqs, const string& errors,
        const mutationMatrices& mutmat, const string& tmp_string
    ) {
        site_id = sid;
        ref_position = position;
        size_t offset = 0;
        site_info = (getIndexFromNucleotide(ref) << 3) + variation_types;
        read_errs.resize(5);
        assert(errors.size() == nucs.size() + insertion_seqs.size() + deletion_seqs.size());


        for (auto i = 0; i < nucs.size(); ++i) {
            read_errs[getIndexFromNucleotide(nucs[i])].push_back(double(errors[i]) - 33.0);
        }
        offset += nucs.size();


        if (variation_types & variationType::INS) {
            for (auto i = 0; i < insertion_seqs.size(); ++i) {
                insertions[insertion_seqs[i]].push_back(double(errors[i + offset] - 33.0));
            }
            offset += insertion_seqs.size();
        }

        if (variation_types & variationType::DEL) {
            for (auto i = 0; i < deletion_seqs.size(); i++) {
                deletions[deletion_seqs[i]].push_back(double(errors[i + offset] - 33.0));
            }
        }

        likelihoods = genotype_likelihoods(read_errs, deletions, insertions, site_info);
        posteriors = genotype_posteriors(likelihoods, deletions, insertions, site_info, mutmat);
        
        for (size_t i = 0; i < 4; i++) {
            read_depth.push_back(read_errs[i].size());
        }
        for (const auto& ins : insertions) {
            read_depth.push_back(ins.second.size());
        }
        for (const auto& del : deletions) {
            read_depth.push_back(del.second.size());
        }

        for (size_t i = 0; i < posteriors.size(); i++) {
            if (posteriors[i] == 0.0) {
                most_probable_idx = i;
                break;
            }
        }

        tmp_readbase_string = tmp_string;
    }

    size_t site_id;
    size_t ref_position;
    int8_t site_info; // 2 bit -> reference nuc, 3 bit -> varaition types
    
    // substitution
    vector< vector<double> > read_errs;

    // deletion
    // map<size_t, size_t> deletions;
    map<string, vector<double> > deletions;
    
    // insertion
    // map<string, size_t> insertions;
    map<string, vector<double> > insertions;

    size_t most_probable_idx;
    vector<double> likelihoods;
    vector<double> posteriors;
    vector<size_t> read_depth;

    string tmp_readbase_string;
};



static void printSiteGenotypeLikelihoods(const variationSite& site) {
    cout << site.ref_position << '\t' << getNucleotideFromIndex(site.site_info >> 3) << "\t";
    for (int i = 0; i < site.read_errs.size(); i++) {
        if (site.likelihoods[i] != numeric_limits<double>::max()) {
            cout << getNucleotideFromIndex(i) << ":" << site.likelihoods[i] << "\t";
        }
    }

    size_t insertion_idx = 0;
    for (const auto& insertion : site.insertions) {
        cout << "+" << insertion.first.size() << insertion.first << ":" << site.likelihoods[5 + insertion_idx] << "\t";
        insertion_idx++;
    }

    size_t deletion_idx = 0;
    for (const auto& deletion : site.deletions) {
        cout << "-" << deletion.first << ":" << site.likelihoods[5 + site.insertions.size() + deletion_idx] << "\t";
        deletion_idx++;
    }
    cout << endl;
}

static void printSiteGenotypePosteriors(const variationSite& site) {
    auto ref_nuc_idx = site.site_info >> 3;
    
    bool no_print = true;
    for (int i = 0; i < site.posteriors.size(); i++) {
        if (i != ref_nuc_idx && site.posteriors[i] != numeric_limits<double>::max()) {
            no_print = false;
        }
    }
    if (no_print) {
        return;
    }
    
    cout << site.ref_position << '\t' << getNucleotideFromIndex(ref_nuc_idx) << "\t";
    for (int i = 0; i < 4; i++) {
        if (site.posteriors[i] != numeric_limits<double>::max()) {
            cout << getNucleotideFromIndex(i) << ":" << site.posteriors[i] << ";";
        }
    }

    size_t insertion_idx = 0;
    for (const auto& insertion : site.insertions) {
        cout << "+" << insertion.first.size() << insertion.first << ":" << site.posteriors[4 + insertion_idx] << ";";
        insertion_idx++;
    }

    size_t deletion_idx = 0;
    for (const auto& deletion : site.deletions) {
        cout << "-" << deletion.first.size() << deletion.first << ":" << site.posteriors[4 + site.insertions.size() + deletion_idx] << ";";
        deletion_idx++;
    }
    
    cout << endl;
    // cout << "\t" << site.tmp_readbase_string << endl;
}

static void printVCFLine(const statsgenotype::variationSite& site) {
    size_t position  = site.ref_position + 1;
    int ref_nuc_idx  = site.site_info >> 3;
    size_t readDepth = 0;
    vector<string> altAlleles;
    vector<size_t> ad;
    vector<int> pl;
    string refAllele;
    int gt;
    
    int quality = site.likelihoods[ref_nuc_idx];
    for (const auto& depth : site.read_depth) {
        readDepth += depth;
    }
   
    refAllele += getNucleotideFromIndex(ref_nuc_idx);

    // find longest deletion
    size_t ldl = 0; // longest deletion length
    string lds;     // longest deletion string
    for (const auto& del : site.deletions) {
        if (del.first.size() > ldl) {
            ldl = del.first.size();
            lds = del.first;
        }
    }

    if (!lds.empty()) {
        refAllele += lds;
    }

    // depth and likelihood for reference
    ad.push_back(site.read_depth[ref_nuc_idx]);
    pl.push_back(site.posteriors[ref_nuc_idx]);
    if (pl.back() == 0.0) {
        gt = 0;
    }

    // substitutions
    for (int i = 0; i < 4; i++) {
        if (i == ref_nuc_idx) {
            continue;
        }
        if (site.posteriors[i] != numeric_limits<double>::max()) {
            altAlleles.push_back(getNucleotideFromIndex(i) + lds);
            ad.push_back(site.read_depth[i]);
            pl.push_back(site.posteriors[i]);
            if (pl.back() == 0.0) {
                gt = altAlleles.size();
            }
        }
    }

    // insertions
    size_t indelIdx = 0;
    for (const auto& ins : site.insertions) {
        altAlleles.push_back(getNucleotideFromIndex(ref_nuc_idx) + ins.first + lds);
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        if (pl.back() == 0.0) {
            gt = altAlleles.size();
        }
        indelIdx += 1;
    }

    // deletions
    for (const auto& del : site.deletions) {
        size_t delSize = del.first.size();
        if (delSize == ldl) {
            altAlleles.push_back(refAllele.substr(0, 1));
        } else {
            altAlleles.push_back(getNucleotideFromIndex(ref_nuc_idx) + lds.substr(delSize, ldl - delSize));
        }
        ad.push_back(site.read_depth[4 + indelIdx]);
        pl.push_back(site.posteriors[4 + indelIdx]);
        if (pl.back() == 0.0) {
            gt = altAlleles.size();
        }
        indelIdx += 1;
    }



    cout << "ref" << "\t"                              // #CHROM
         << position << "\t"                           // POS
         << "." << "\t"                                // ID
         << refAllele << "\t";                         // REF
    for (size_t i = 0; i < altAlleles.size() - 1; i++) {
        cout << altAlleles[i] << ",";
    }                                    
    cout << altAlleles[altAlleles.size() - 1] << "\t"; // ALT
    cout << quality << "\t"                            // QUAL
         << "." << "\t"                                // FILTER
         << "DP=" << readDepth << "\t"                 // INFO
         << "GT:AD:PL" << "\t"                         // FORMAT
         << gt << ":";                                 // SAMPLE
    for (size_t i = 0; i < ad.size() - 1; i++) {
        cout << ad[i] << ",";
    }
    cout << ad[ad.size() - 1] << ":";
    for (size_t i = 0; i < pl.size() - 1; i++) {
        cout << pl[i] << ",";
    }
    cout << pl[pl.size() - 1] << "\t\t" << site.tmp_readbase_string << endl;
}

}


#endif //STATS_GENOTYPE
