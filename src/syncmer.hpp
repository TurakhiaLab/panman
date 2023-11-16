#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>


struct kmer_t {
    
    std::string seq;
    size_t pos;
    int32_t idx; // index of kmer in vector used in indexing DFS
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
    bool operator==(const kmer_t &rhs) const {
        return seq == rhs.seq;
    };
};
class KHash {
public:
    size_t operator()(const kmer_t& t) const
    {
        const std::hash<std::string> h;
        return h(t.seq);
    }
    
};

struct read_t {
    std::string seq;
    std::vector<kmer_t> kmers;
    //std::vector<int> read_coord;
    //std::vector<int> ref_coord;
    //std::vector<bool> reversed;
};

struct seedIndex {
    std::vector<kmer_t> rootSeeds;
    std::unordered_map<std::string, std::vector<kmer_t>> deletions;
    std::unordered_map<std::string, std::vector<kmer_t>> insertions;
}; 

struct dynamicJaccard {
    size_t intersectionSize;
    size_t unionSize;
    float jaccardIndex;
};

inline bool is_syncmer(std::string &seq, int s, bool open) {
    if (seq.size() < s) {
        return false;
    }
    std::string min(s, 'Z');
    for (size_t i = 0; i < seq.size() -  s + 1; i++) {
        std::string submer = seq.substr(i, s);
        if (submer < min) {
            min = submer;
        }
    }
    if (open) {
        if (min == seq.substr(0, s)) {
            return true;
        }
    } else {
        if (min == seq.substr(0, s) || min == seq.substr(seq.length()-s, s)) {
            return true;
        }
    }
    
    return false;
}

inline std::vector<kmer_t> syncmerize(std::string seq, int k, int s, bool open, bool aligned, int pad) {
    std::vector<kmer_t> ret;
    if (aligned) {
        std::unordered_map<int32_t, int32_t> degap;
        int32_t pos = 0;
        std::string ungapped = "";
        for (int32_t i = 0; i < seq.size(); i++) {
            char c = seq[i];
            degap[pos] = i;
            if (c != '-') {
                ungapped += c;
                pos++;
            }
        }
    // 012345 6789
        if (ungapped.size() < k + 1) {
            return ret;
        }
        for (int32_t i = 0; i < ungapped.size() - k + 1; i++) {
            std::string kmer = ungapped.substr(i, k);
            if (is_syncmer(kmer, s, open)) {
                ret.push_back(kmer_t{kmer, degap[i]+pad, -1});
            }
        }
    } else {
        for (int32_t i = 0; i < seq.size() - k + 1; i++) {
            std::string kmer = seq.substr(i, k);
            if (is_syncmer(kmer, s, open)) {
                ret.push_back(kmer_t{kmer, i+pad, -1});
            }
        }
    }
    return ret;
}
