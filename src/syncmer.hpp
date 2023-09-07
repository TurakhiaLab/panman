#include <iostream>
#include <set>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>


struct kmer_t {
    std::string seq;
    size_t pos;
    size_t pos2;
    bool reversed;
    bool operator<(const kmer_t &rhs) const {
        if (seq == "") {
            return false;
        } else if (rhs.seq == "") {
            return true;
        } else {
            return seq < rhs.seq;
        }
    };
};


struct read_t {
    std::string seq;
    std::set<kmer_t> kmers;
    //std::vector<int> read_coord;
    //std::vector<int> ref_coord;
    //std::vector<bool> reversed;
};


struct seedIndex {
    std::set<kmer_t> rootSeeds;
    std::unordered_map<std::string, std::set<size_t>> deletions;
    std::unordered_map<std::string, std::set<kmer_t>> insertions;
};

struct dynamicJaccard {
    size_t intersectionSize;
    size_t unionSize;
    float jaccardIndex;
};

inline bool is_syncmer(std::string seq, int s, bool open) {
    std::set<std::string> submers;
    for (size_t i = 0; i < seq.size() -  s + 1; i++) {
        std::string submer = seq.substr(i, s);
        submers.insert(submer);
    }
    if (open) {
        if ((*(submers.begin())) == seq.substr(0, s)) {
            return true;
        }
    } else {
        if ((*(submers.begin())) == seq.substr(0, s) || (*(submers.begin())) == seq.substr(seq.length()-s, s)) {
            return true;
        }
    }
    
    return false;
}

inline std::set<kmer_t> syncmerize(std::string seq, int k, int s, bool open, bool aligned) {
    std::set<kmer_t> ret;
    if (aligned) {
        std::unordered_map<size_t, size_t> degap;
        size_t pos = 0;
        std::string ungapped = "";
        for (size_t i = 0; i < seq.size(); i++) {
            char c = seq[i];
            if (c != '-') {
                pos++;
                ungapped += c;
            }
            degap[i] = pos;
        }
    
        for (size_t i = 0; i < ungapped.size() - k + 1; i++) {
            std::string kmer = ungapped.substr(i, k);
            if (is_syncmer(kmer, s, open)) {
                ret.insert(kmer_t{kmer, degap[i]});
            }
        }
    } else {
        for (size_t i = 0; i < seq.size() - k + 1; i++) {
            std::string kmer = seq.substr(i, k);
            if (is_syncmer(kmer, s, open)) {
                ret.insert(kmer_t{kmer, i});
            }
        }
    }
    return ret;
}
