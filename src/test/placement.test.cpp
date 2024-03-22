 #define BOOST_TEST_MODULE Index construction
#include <boost/test/included/unit_test.hpp>
#include "../PangenomeMAT.hpp"
#include "../statsgenotype.hpp"
#include <stack>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <json/json.h>
#include <filesystem>
#include <unordered_set>
#include <regex>
#include <random>
#include <algorithm>

namespace fs = boost::filesystem;


using namespace PangenomeMAT;

// BOOST_AUTO_TEST_CASE(_removeIndices) {
// /* should erase elements at the given positions */
//     std::vector<kmer_t> a = {
//         kmer_t{"zero"},
//         kmer_t{"one"},
//         kmer_t{"two"},
//         kmer_t{"three"},
//         kmer_t{"four"},
//         kmer_t{"five"}      
//     };
    
//     // smaller values at the top
//     std::stack<int32_t> indices;

//     indices.push(4);
//     indices.push(2);
//     indices.push(0);

//     std::vector<kmer_t> expected = {
//         kmer_t{"one"},
//         kmer_t{"three"},
//         kmer_t{"five"}
//     };

//     removeIndices(a, indices);

//     BOOST_TEST(
//         a == expected
//     );

// }

// BOOST_AUTO_TEST_CASE(_isSyncmer) {
// /* should return true only for syncmers */
//     std::string a = "ABCDEF";

//     BOOST_TEST(
//         is_syncmer(a, 3, true) == true
//     );
//     BOOST_TEST(
// 		is_syncmer(a, 3, false) == true
//     );
    
//     std::string b = "ZZZAAAZZZ";

//     BOOST_TEST(
// 		is_syncmer(b, 6, true) == false
//     );
//     BOOST_TEST(
// 		is_syncmer(b, 6, false) == true
//     );

//     std::string c = "Z";
 
//     BOOST_TEST(
// 		is_syncmer(c, 1, true) == true
//     );
// }

// BOOST_AUTO_TEST_CASE(_syncmerize) {
// /* should pick syncmers correctly */
//     std::string a = "AAAAAABB"; // k=3 s=2 open
//     /*             * AAA
//                    *  AAA
//                    *   AAA
//                    *    AAA
//                    *     AAB
//                    *      ABB
//     */
//     std::vector<kmer_t> ret_a = {
//         kmer_t{"AAA", 0, -1},
//         kmer_t{"AAA", 1, -1},
//         kmer_t{"AAA", 2, -1},
//         kmer_t{"AAA", 3, -1},
//         kmer_t{"AAB", 4, -1},
//         kmer_t{"ABB", 5, -1}                
//     };
//     std::vector<kmer_t> ret_a_padded = {
//         kmer_t{"AAA", 42, -1},
//         kmer_t{"AAA", 43, -1},
//         kmer_t{"AAA", 44, -1},
//         kmer_t{"AAA", 45, -1},
//         kmer_t{"AAB", 46, -1},
//         kmer_t{"ABB", 47, -1}                
//     };
//     BOOST_TEST(
// 		syncmerize(a, 3, 2, true, false, 0) == ret_a
//     );
//     BOOST_TEST(
// 		syncmerize(a, 3, 2, true, false, 42) == ret_a_padded
//     );

//     std::string b = "AABCBAAACC"; // k=6 s=3 open
//   /*               * AABCBA 
//                    *  ABCBAA
//                        BCBAAA
//                         CBAAAC
//                          BAAACC
//   */
//     std::vector<kmer_t> ret_b = {
//         kmer_t{"AABCBA", 0, -1},
//         kmer_t{"ABCBAA", 1, -1}
//     };
//     std::string c = "AABCBAAACC"; // k=6 s=3 closed
//   /*               * AABCBA 
//                    *  ABCBAA
//                    *   BCBAAA
//                         CBAAAC
//                          BAAACC
//   */
//     std::vector<kmer_t> ret_c = {
//         kmer_t{"AABCBA", 0, -1},
//         kmer_t{"ABCBAA", 1, -1},
//         kmer_t{"BCBAAA", 2, -1}
//     };

//     BOOST_TEST(
// 		syncmerize(b, 6, 3, true, false, 0) == ret_b
//     );

//     std::cout << "ret c: " << std::endl;
//     for (auto kmer : ret_c) {
//         std::cout << kmer.seq << " " << kmer.pos << std::endl;
//     }
//     std::cout << "syncmerize c: " << std::endl;
//     for (auto kmer : syncmerize(c, 6, 3, false, false, 0)) {
//         std::cout << kmer.seq  << " " << kmer.pos << std::endl;
//     }
//     BOOST_TEST(
// 		syncmerize(c, 6, 3, false, false, 0) == ret_c
//     );

//     std::string d = "--A-A-B--CBAAAC-C-"; // k=6 s=3 closed aligned
//   /*               * AABCBA 
//                    *  ABCBAA
//                    *   BCBAAA
//                         CBAAAC
//                          BAAACC
//   */
//      std::vector<kmer_t> ret_d = {
//         kmer_t{"AABCBA", 2, -1},
//         kmer_t{"ABCBAA", 4, -1},
//         kmer_t{"BCBAAA", 6, -1}
//     };
//     BOOST_TEST(
// 		syncmerize(d, 6, 3, false, true, 0) == ret_d
//     );

// }


// template<>
// struct boost::test_tools::tt_detail::print_log_value<std::pair<int,int>> {
//     void operator()(std::ostream& os, std::pair<int,int> const& value) {
//         os << '{' << value.first << ',' << value.second << '}';
//     }
// };

// BOOST_AUTO_TEST_CASE(_getRecomputePositions)
// {
//     int32_t k = 3;

//     std::pair<int32_t, int32_t> a = {6, 1};
//     std::pair<int32_t, int32_t> b = {22, 3};
//     std::pair<int32_t, int32_t> c = {31, 2};

//     //                 length:  0               3        2
//     //                          *               ***      **  
//     //      recompute range: xxxxxxx       xxxxxxxxxx  xxxxx
//     std::string gapped = "ABCD-EF-GH-IJKLM-N-O--P-QRSTUVWXYZ";
//     //                    0123456789111111111122222222223333
//     //                              012345678901234567890123

//     std::pair<int32_t, int32_t> expected_a = {3, 9};
//     std::pair<int32_t, int32_t> expected_b = {17, 26};
//     std::pair<int32_t, int32_t> expected_c = {29, 33};

//     BOOST_TEST(
//         getRecomputePositions(a, gapped, k) == expected_a
//     );
//     BOOST_TEST(
//         getRecomputePositions(b, gapped, k) == expected_b
//     );
//     BOOST_TEST(
//         getRecomputePositions(c, gapped, k) == expected_c
//     );
// }

// BOOST_AUTO_TEST_CASE(_alignedEndPos) {
//  /* should compute a syncmer's end pos with gaps */

//     std::string a = "--AAAAA--A-A-A-----A--AAA--";
//     //               012345678911111111112222222
//     //                         01234567890123456
//     //                            ^^^^^^^^^^^

//     BOOST_TEST(
//       alignedEndPos(13, 4, a) ==  23 
//     );
// }

// BOOST_AUTO_TEST_CASE(_discardSyncmers) {
// /*  should discard syncmers contained in any range in B 
//     unless they are going to be inserted. Then remove from
//     to_insert (based on sequence for now). Also, track variable syncmers.
// */
//     size_t k = 4;

//     std::vector<kmer_t> seeds = {
//         kmer_t{"AAAA", 0, -1}, // x
//         kmer_t{"ACCA", 2, -1}, //
//         kmer_t{"AABA", 7, -1}, // x
//         kmer_t{"AACA", 8, -1}, // x
//         kmer_t{"ABDA", 9, -1}, //
//         kmer_t{"ABBA", 15,-1}, // gapped length exceeds bound  
//         kmer_t{"CDEF", 23,-1}, 
//     };
//     std::string seq = "xxxxxxxxxxxxxxxABB---Axxxxxx"; // len 28

//     std::vector<std::pair<int32_t, int32_t>> B = {
//         {0, 4},
//         {7, 11},
//         {14, 20},
//         {22, 35}
//     };
//     std::vector<kmer_t> expected_seeds = {
//         kmer_t{"ACCA", 2, -1},
//         kmer_t{"ABDA", 9, -1},
//         kmer_t{"ABBA", 15,-1},
//         kmer_t{"CDEF", 23,-1}, 
//     };

//     std::unordered_map<std::string, kmer_t> to_insert = {
//         {"XYZW", kmer_t{"XYZW", 8, -1}},
//         {"CDEF", kmer_t{"CDEF", 23, -1}}
//     };
//     std::unordered_map<std::string, kmer_t> expected_to_insert = {
//         {"XYZW", kmer_t{"XYZW", 8, -1}}
//     };

//     seedIndex index;
    
//     std::vector<kmer_t> expected_deletions = {
//                 kmer_t{"AACA", 8, 3},
//                 kmer_t{"AABA", 7, 2},
//                 kmer_t{"AAAA", 0, 0},
//     };
//     std::string nid = "my_node_id";

//     std::unordered_map<std::string, bool> variable_syncmers;

//     std::unordered_map<std::string, bool> expected_variable_syncmers = {
//         {"AAAA", true},
//         {"AABA", true},
//         {"AACA", true}
//     };

//     discardSyncmers(seeds, B, seq, to_insert, variable_syncmers, index, nid, k);

//     BOOST_TEST(
//         seeds == expected_seeds
//     );
//     BOOST_TEST(
//         to_insert == expected_to_insert
//     );

//     BOOST_TEST(
//         index.deletions[nid] == expected_deletions
//     );
//     BOOST_TEST(index.deletions[nid][0].pos == 8);
//     BOOST_TEST(index.deletions[nid][0].idx == 3);
//     BOOST_TEST(index.deletions[nid][1].pos == 7);
//     BOOST_TEST(index.deletions[nid][1].idx == 2);
//     BOOST_TEST(index.deletions[nid][2].pos == 0);
//     BOOST_TEST(index.deletions[nid][2].idx == 0);

//     BOOST_TEST(
//         variable_syncmers == expected_variable_syncmers
//     );
// }

// BOOST_AUTO_TEST_CASE(_indexSyncmers) {
//     size_t k = 13;
//     size_t s = 4;
//     std::ifstream is("../src/test/test.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSyncmers(T, os, k, s);
// }

// BOOST_AUTO_TEST_CASE(_indexSyncmers) {
//     size_t k = 13;
//     size_t s = 4;
//     std::ifstream is("../src/test/test.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSyncmers(T, os, k, s);
// }

// BOOST_AUTO_TEST_CASE(accuracy) {
//     std::vector<std::string> files = {"../fastq/MT821778.1.1kreads.fastq_R1.fastq","../fastq/MT834260.1.1kreads.fastq_R1.fastq","../fastq/MT969389.1.1kreads.fastq_R1.fastq","../fastq/MW211005.1.1kreads.fastq_R1.fastq","../fastq/MW593420.1.1kreads.fastq_R1.fastq","../fastq/MW848017.1.1kreads.fastq_R1.fastq","../fastq/MW848911.1.1kreads.fastq_R1.fastq"};
//     size_t k = 13;
//     size_t s = 7;
//     std::ifstream is("../sars2k.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSyncmers(T, os, k, s);
    
//     is.close();
//     os.close();
//     PangenomeMAT::Node *root = T->root;
//     struct seedIndex index;
//     std::ifstream indexFile("./test.out");
//     PangenomeMAT::loadIndex(T->root, indexFile, index);
    
//     for (std::string f : files) {

//         std::vector<read_t> reads;
//         auto fastq_start = std::chrono::high_resolution_clock::now();
//         std::set<kmer_t> readSyncmers = syncmersFromFastq(f, reads, k, s);
//         auto fastq_end = std::chrono::high_resolution_clock::now();

//         std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

//         auto place_start = std::chrono::high_resolution_clock::now();

//         std::set<kmer_t> rootSyncmers = std::set<kmer_t>(index.rootSeeds.begin(), index.rootSeeds.end());

//         std::cerr << "\n";
//         std::cerr << "Placing sample... " << f << "\n";


//         struct dynamicJaccard dj;
    
//         dj.intersectionSize = intersection_size(rootSyncmers, readSyncmers);
//         dj.unionSize = rootSyncmers.size() + readSyncmers.size() - dj.intersectionSize;
//         dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
        
        
//         std::cout << "root seeds: " << rootSyncmers.size() << "\n";
//         std::cout << "read seeds: " << readSyncmers.size() << "\n";
//         for (const auto &k : readSyncmers) {
//             std::cout << k.seq << "\n";
//         }
//         std::cout << "initial jaccard: " << dj.jaccardIndex << "\n";

//         std::unordered_map<std::string, float> scores;
//         std::unordered_map<std::string, bool> readSyncmersMap;
//         for (const auto &k : readSyncmers) {
//             readSyncmersMap[k.seq] = true;
//         }

//         placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);

//         auto place_end = std::chrono::high_resolution_clock::now();

//         std::cout << "place time: " << std::chrono::duration_cast<std::chrono::milliseconds>(place_end - place_start).count() << "\n";

//         std::vector<std::pair<std::string, float>> v;
//         for ( const auto &p : scores ) {
//             v.push_back(std::make_pair(p.first, p.second));
//         } 
//         std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
//             return left.second > right.second;
//         });

//         std::string best_match = v[0].first;
//         for (const auto &s : v) {
//             std::cerr << s.first << ": " << s.second << "\n";
//         }

//     }
// }



vector< pair<string, string> > get_pruned_samples(string path) {
    vector< pair<string, string> > samples;
    for (const auto& file : filesystem::directory_iterator("../src/test/statsgenotype/pmat/")) {
        string pruned_tree_path = file.path().string();
        vector<string> fields;
        PangenomeMAT::stringSplit(pruned_tree_path, '/', fields);
        string pmat_filename = fields.back();
        fields.clear();
        PangenomeMAT::stringSplit(pmat_filename, '_', fields);
        string pruned_sample = fields[0];
        samples.push_back(make_pair(pruned_sample, pruned_tree_path));
    }

    return samples;
}



void pileup(string pluPath, string parentFastaPath, string samPath) {
    if(!filesystem::exists(pluPath)) {
        string cmd = "samtools mpileup " + samPath + " -f " + parentFastaPath + " -A > " + pluPath;
        system(cmd.c_str());
    }
}

void makeFasta(Tree* T, string nodeName, string path) {
    const std::string seq = T->getStringFromReference(nodeName, false);
    if(!filesystem::exists(path)) {
        ofstream faos(path);
        faos << '>' << nodeName << '\n';
        size_t linesize = 80;
        for (size_t i = 0; i < seq.size(); i += linesize) {
            faos << seq.substr(i, std::min(linesize, seq.size() - i)) << '\n';
        }
        faos.close();
    }
}

void simReads(string refFastaPath, string path, size_t numReads) {
    string numReadsString = to_string(numReads);
    if(!filesystem::exists(path + "_R1.fastq")) {
        string cmd = "iss generate --model NovaSeq --genomes " + refFastaPath + " -n " + numReadsString + " --output " + path + " --cpus 4";
        cout << cmd << endl;
        system(cmd.c_str());
    }
}

void mapReads(string samPath, string refFastaPath, string nodeName) {
    if(!filesystem::exists(samPath)) {
        string cmd = "minimap2 -ax sr " + refFastaPath + " ../src/test/statsgenotype/fastq_2k/" + nodeName + "*fastq_R* --heap-sort=yes | samtools sort -O sam > " + samPath;
        system(cmd.c_str());
    }
}

pair<size_t, size_t> getEdgeCoor(Tree* T, string identifier) {
    const string& seq = T->getStringFromReference(identifier);
    size_t count = 0;
    size_t start;
    size_t end;
    for (size_t i = 0; i < seq.size(); i++) {
        if (seq[i] != '-') {
            count += 1;
        }
        if (count >= 1000) {
            start = i;
            break;
        }
    }
    count = 0;
    for (size_t i = seq.size()-1; i >= 0; i--) {
        if (seq[i] != '-') {
            count += 1;
        }
        if (count >= 1000) {
            end = i;
            break;
        }
    }
    return make_pair(start, end);
}

int maskedSize(Tree* T, const vector<NucMut>& nucMutation, const pair<size_t, size_t>& edges) {
    vector<NucMut> masked;
    for (const auto& mut : nucMutation) {
        size_t globalCoor = T->getGlobalCoordinate(mut.primaryBlockId, mut.secondaryBlockId, mut.nucPosition, mut.nucGapPosition);
        if (globalCoor > edges.first && globalCoor < edges.second) {
            masked.push_back(mut);
        }
    }
    
    int size = 0;
    for (const auto& mut: masked) {
        if ((mut.mutInfo & 15) == 0) {
            size += (mut.mutInfo >> 4);
        } else {
            size += 1;
        }
    }
    return size;
}

size_t get_distance(Tree* T, Node* node_1, Node* node_2) {
    size_t n1_dist = 0;
    size_t n2_dist = 0;
    while (node_1->identifier != node_2->identifier) {
        if (node_1->level == node_2->level) {
            // n1_dist += node_1->nucMutation.size();
            // n2_dist += node_2->nucMutation.size();
            n1_dist += maskedSize(T, node_1->nucMutation, getEdgeCoor(T, node_1->identifier));
            n2_dist += maskedSize(T, node_2->nucMutation, getEdgeCoor(T, node_2->identifier));
            node_1 = node_1->parent;
            node_2 = node_2->parent;
        } else if (node_1->level > node_2->level) {
            // n1_dist += node_1->nucMutation.size();
            n1_dist += maskedSize(T, node_1->nucMutation, getEdgeCoor(T, node_1->identifier));
            node_1 = node_1->parent;
        } else if (node_2->level > node_1->level) {
            // n2_dist += node_2->nucMutation.size();
            n2_dist += maskedSize(T, node_2->nucMutation, getEdgeCoor(T, node_2->identifier));
            node_2 = node_2->parent;
        }
    }
    //cout << n1_dist << "\t" << n2_dist << endl;
    size_t distance = n1_dist + n2_dist;
    return distance;
}

size_t get_distance_branch(Tree* T, Node* node_1, Node* node_2) {
    size_t n1_dist = 0;
    size_t n2_dist = 0;
    while (node_1->identifier != node_2->identifier) {
        if (node_1->level == node_2->level) {
            // n1_dist += (node_1->nucMutation.size() == 0) ? 0 : 1;
            // n2_dist += (node_2->nucMutation.size() == 0) ? 0 : 1;
            n1_dist += (maskedSize(T, node_1->nucMutation, getEdgeCoor(T, node_1->identifier)) == 0) ? 0 : 1;
            n2_dist += (maskedSize(T, node_2->nucMutation, getEdgeCoor(T, node_2->identifier)) == 0) ? 0 : 1;
            node_1 = node_1->parent;
            node_2 = node_2->parent;
        } else if (node_1->level > node_2->level) {
            // n1_dist += (node_1->nucMutation.size() == 0) ? 0 : 1;
            n1_dist += (maskedSize(T, node_1->nucMutation, getEdgeCoor(T, node_1->identifier)) == 0) ? 0 : 1;
            node_1 = node_1->parent;
        } else if (node_2->level > node_1->level) {
            // n2_dist += (node_2->nucMutation.size() == 0) ? 0 : 1;
            n2_dist += (maskedSize(T, node_2->nucMutation, getEdgeCoor(T, node_2->identifier)) == 0) ? 0 : 1;
            node_2 = node_2->parent;
        }
    }
    //cout << n1_dist << "\t" << n2_dist << endl;
    size_t distance = n1_dist + n2_dist;
    return distance;
}



char getRandomChar(const std::vector<char>& charList) {
    // Create a random device and a generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a distribution in the range of the vector indices
    std::uniform_int_distribution<> distr(0, charList.size() - 1);

    // Get a random index
    int index = distr(gen);

    // Return the character at the random index
    return charList[index];
}

char getRandomCharWithWeights(const std::vector<char>& chars, const std::vector<int>& weights) {
    // Create a random device and generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a discrete distribution with the given weights
    std::discrete_distribution<> distr(weights.begin(), weights.end());

    // Get a random index based on the weights
    int index = distr(gen);

    // Return the character at the random index
    return chars[index];
}

double getMinDouble(const vector<double>& doubles) {
    if (doubles.empty()) {
        return -1.0;
    }

    double minDouble = doubles[0];
    for (size_t i = 1; i < doubles.size(); i++) {
        minDouble = min(minDouble, doubles[i]);
    }

    return minDouble;

}

char subNuc(char ref, const statsgenotype::mutationMatrices& mutmat) {
    vector<char> bases = {'A', 'C', 'G', 'T'};
    int refIdx = statsgenotype::getIndexFromNucleotide(ref);

    if (refIdx > 3) {
        return getRandomChar(bases);
    }
    
    vector<double> probs = mutmat.submat[refIdx];
    bases.erase(bases.begin() + refIdx);
    probs.erase(probs.begin() + refIdx);
    
    vector<int> wgts;
    double minProb = getMinDouble(probs);
    for (const auto& prob : probs) {
        double deci = pow(10, (minProb - prob) / 10);
        wgts.push_back(deci * 1000);
    }

    return getRandomCharWithWeights(bases, wgts);
}

void genMut2(
    const string& seq,  const statsgenotype::mutationMatrices& mutmat,
    const string& prefix, size_t beg, size_t end, size_t num, size_t len,
    vector<int> vtp, const string& outDir
){
    if(filesystem::exists((outDir + prefix + ".fa"))) {
        return;
    }

    vector<char> bases = {'A', 'C', 'G', 'T'};
    ofstream faos(outDir + prefix + ".fa");
    ofstream vros(outDir + prefix + ".rf");

    string nseq =  seq;
    vector< tuple<int, int, int> > muts;
    vector<string> vref;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distribPos(beg, seq.size() - end);
    uniform_int_distribution<> distribLen(1, len);
    int c = 0;
    while (muts.size() < num) {
        int curPos  = distribPos(gen);
        int varType = vtp[c % vtp.size()];
        int length = (varType == 1) ? 1 : distribLen(gen);
        
        bool posConflict = false;
        for (const auto& mut : muts) {
            if (get<1>(mut) == 4 && abs(curPos - get<0>(mut)) < (2 * get<2>(mut))) {
                posConflict = true;
                break;
            } else if (abs(curPos - get<0>(mut)) == 0) {
                posConflict = true;
                break;
            }

        }
        if (posConflict) {
            continue;
        }
        muts.emplace_back(make_tuple(curPos, varType, length));
        c++;
    }

    sort(muts.begin(), muts.end(), 
        [](const tuple<int, int, int>& a, const tuple<int, int, int>& b) {
            return get<0>(a) < std::get<0>(b);
        }
    );

    int offset = 0;
    for (size_t i = 0; i < num; i++) {
        int pos = get<0>(muts[i]);
        int var = get<1>(muts[i]);
        
        switch (var) {
            case 1:
                {
                    vros << to_string(pos + 1) + "\t";
                    char ref  = seq[pos];
                    char mut  = subNuc(ref, mutmat);
                    nseq[pos + offset] = mut;
                    vros << ref << "\t" << mut << "\n";
                    break;
                }
            case 2:
                {
                    vros << to_string(pos + 1) + "\t";
                    char ref = seq[pos];
                    string inss;
                    for (size_t j = 0; j < get<2>(muts[i]); j++) {
                        inss += getRandomChar(bases);
                    }
                    nseq = nseq.substr(0, pos + offset + 1) + inss + nseq.substr(pos + offset + 1, nseq.size() - (pos + offset + 1));
                    vros << ref << "\t" << ref << inss << "\n";
                    offset += get<2>(muts[i]);
                    break;
                }
            case 4:
                {
                    vros << to_string(pos + 1) + "\t";
                    string ref = seq.substr(pos, get<2>(muts[i]) + 1);
                    char   del = seq[pos];
                    nseq = nseq.substr(0, pos + offset + 1) + nseq.substr(pos + offset + get<2>(muts[i]) + 1, nseq.size() - (pos + offset + get<2>(muts[i]) + 1));
                    vros << ref << "\t" << del << "\n";
                    offset -= get<2>(muts[i]);
                    break;
                }
        }
    }
    
    int type = 0;
    for (auto tp : vtp) { type += tp; }

    string seqName = to_string(num) + to_string(len) + to_string(type);
    faos << '>' << seqName << '\n';
    size_t linesize = 80;
    for (size_t i = 0; i < nseq.size(); i += linesize) {
        faos << nseq.substr(i, std::min(linesize, nseq.size() - i)) << '\n';
    }
    
    vros.close();
    faos.close();
}

void genMut(
    const string& seq,  const statsgenotype::mutationMatrices& mutmat,
    const string& prefix, size_t beg, size_t end, size_t num, size_t len,
    vector<int> vtp
){
    if(filesystem::exists(("../src/test/statsgenotype/validation/muts/" + prefix + ".fa"))) {
        return;
    }

    vector<char> bases = {'A', 'C', 'G', 'T'};
    ofstream faos("../src/test/statsgenotype/validation/muts/" + prefix + ".fa");
    ofstream vros("../src/test/statsgenotype/validation/muts/" + prefix + ".rf");

    string nseq =  seq;
    vector< pair<int, int> > muts;
    vector<string> vref;

    // randomly generate mutations (2 * len + 100) bases apart
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(beg, seq.size() - end);
    int c = 0;
    while (muts.size() < num) {
        int curPos = distrib(gen);
        for (const auto mut : muts) {
            if (abs(curPos - mut.first) < (2 * len + 100)) {
                continue;
            }
        }
        muts.emplace_back(curPos, vtp[c % vtp.size()]);
        c++;
    }
    sort(muts.begin(), muts.end());

    int offset = 0;
    for (size_t i = 0; i < num; i++) {
        int pos = muts[i].first;
        int var = muts[i].second;
        
        switch (var) {
            case 1:
                {
                    for (size_t j = 0; j < len; j++) {
                        vros << to_string(pos + 1 + j) + "\t";
                        char ref  = seq[pos + j];
                        char mut  = subNuc(ref, mutmat);
                        nseq[pos + j + offset] = mut;
                        vros << ref << "\t" << mut << "\n";
                    }
                    break;
                }
            case 2:
                {
                    vros << to_string(pos + 1) + "\t";
                    char ref = seq[pos];
                    string inss;
                    for (size_t j = 0; j < len; j++) {
                        inss += getRandomChar(bases);
                    }
                    nseq = nseq.substr(0, pos + offset + 1) + inss + nseq.substr(pos + offset + 1, nseq.size() - (pos + offset + 1));
                    vros << ref << "\t" << ref << inss << "\n";
                    offset += len;
                    break;
                }
            case 4:
                {
                    vros << to_string(pos + 1) + "\t";
                    string ref = seq.substr(pos, len + 1);
                    char   del = seq[pos];
                    nseq = nseq.substr(0, pos + offset + 1) + nseq.substr(pos + offset + len + 1, nseq.size() - (pos + offset + len + 1));
                    vros << ref << "\t" << del << "\n";
                    offset -= len;
                    break;
                }
        }
    }
    
    int type = 0;
    for (auto tp : vtp) { type += tp; }

    string seqName = to_string(num) + to_string(len) + to_string(type);
    faos << '>' << seqName << '\n';
    size_t linesize = 80;
    for (size_t i = 0; i < nseq.size(); i += linesize) {
        faos << nseq.substr(i, std::min(linesize, nseq.size() - i)) << '\n';
    }
    
    vros.close();
    faos.close();
}

map<size_t, pair<string, string> > parseSim(ifstream& simIn) {
    map<size_t, pair<string, string> > simVar;
    string line;
    while(getline(simIn, line)) {
        vector<string> fields;
        PangenomeMAT::stringSplit(line, '\t', fields);
        size_t pos = stoul(fields[0]);
        string ref = fields[1];
        string var = fields[2];
        simVar[pos] = make_pair(ref, var);
    }
    return simVar;
}

map<size_t, pair< vector<string>, vector<int> > > parseVcf(ifstream& vcfIn) {
    map<size_t, pair< vector<string>, vector<int> > > vcfVar;
    string line;
    while(getline(vcfIn, line)) {
        if (line[0] == '#') { continue; }
        vector<string> gts;
        vector<int>    gps;
        vector<string> fields;
        PangenomeMAT::stringSplit(line, '\t', fields);
        size_t pos = stoul(fields[1]);
        gts.push_back(fields[3]);

        vector<string> altAlleles;
        PangenomeMAT::stringSplit(fields[4], ',', altAlleles);
        for(const auto& alt : altAlleles) { gts.push_back(alt); }

        vector<string> sampleInfo;
        PangenomeMAT::stringSplit(fields[9], ':', sampleInfo);
        vector<string> gpsStr;
        PangenomeMAT::stringSplit(sampleInfo[2], ',', gpsStr);
        for (const auto& gp : gpsStr) {
            if (gp == ".") {
                gps.push_back(999);
            } else {
                gps.push_back(stoul(gp));
            }
        }

        vcfVar[pos] = make_pair(gts, gps);
    }
    return vcfVar;
}

vector<size_t> parseDepths(ifstream& pluIn, size_t genomeSize) {
    vector<size_t> depths(genomeSize + 1);
    string line;
    while(getline(pluIn, line)) {
        vector<string> fields;
        PangenomeMAT::stringSplit(line, '\t', fields);
        size_t pos = stoul(fields[1]);
        size_t dep = stoul(fields[3]);
        depths[pos] = dep;
    }
    return depths;
}

size_t size_tDiff(size_t a, size_t b) {
    if (a > b) {
        return a - b;
    }
    
    return b - a;
}

bool isVisited(const unordered_set<size_t>& uoset, size_t position) {
    if (uoset.find(position) != uoset.end()) {
        return true;
    }
    return false;
}

bool isoAlignment(
    const pair<size_t, pair<string, string> >& simVar, 
    const pair<size_t, pair< vector<string>, vector<int> > >& vcfVar, 
    const string& refSeq
){
    bool ins =  (simVar.second.second.size() > simVar.second.first.size()) ? true : false;
    
    size_t gtIdx;
    for (size_t i = 0; i < vcfVar.second.second.size(); i++) {
        if (vcfVar.second.second[i] == 0) {
            gtIdx = i;
            break;
        }
    }
    
    size_t varSize = (ins) ? simVar.second.second.size() - 1 : simVar.second.first.size() - 1;

    if (size_tDiff(simVar.first, vcfVar.first) > varSize * 2) {
        return false;
    }
    


    size_t beg = min(simVar.first - 1, vcfVar.first - 1);
    size_t end = max(simVar.first - 1, vcfVar.first - 1);

    string simLocSeq;
    string varLocSeq;
    if (ins) {
        end += simVar.second.second.size();
        simLocSeq = refSeq.substr(beg, simVar.first - 1 - beg) + simVar.second.second       + refSeq.substr(simVar.first, end - beg - (simVar.first - 1 - beg) - simVar.second.second.size() + 1);
        varLocSeq = refSeq.substr(beg, vcfVar.first - 1 - beg) + vcfVar.second.first[gtIdx] + refSeq.substr(vcfVar.first, end - beg - (vcfVar.first - 1 - beg) - simVar.second.second.size() + 1);
    } else {
        end += simVar.second.first.size();
        simLocSeq = refSeq.substr(beg, simVar.first - 1 - beg) + simVar.second.second       + refSeq.substr(simVar.first + simVar.second.first.size() - 1, end - beg - (simVar.first - 1 - beg) - simVar.second.first.size() + 1);
        varLocSeq = refSeq.substr(beg, vcfVar.first - 1 - beg) + vcfVar.second.first[gtIdx] + refSeq.substr(vcfVar.first + simVar.second.first.size() - 1, end - beg - (vcfVar.first - 1 - beg) - simVar.second.first.size() + 1);
    }


    if (simLocSeq == varLocSeq) {
        return true;
    }

    return false;
}

struct TestStats {
    TestStats(int tp, int fp, int tn, int fn, int refSize) {
        truePositive  = tp;
        falsePositive = fp;
        trueNegative  = tn;
        falseNegative = fn;

        trueNegative += (refSize - tp - fp - tn - fn);
    }

    int truePositive;
    int falsePositive;
    int trueNegative;
    int falseNegative;
};

TestStats runTest(
    Tree* T, const string& refSeq, string mutPrefix, string refPath, size_t coverage,
    string workingDir, map<size_t, vector<int> >& byBaseCoverage
) {
    // workingDir = "../src/test/statsgenotype/validation/";
    
    // make reads
    string numReads = to_string(200 * coverage);
    string readsCmd = "iss generate --model NovaSeq --genomes ../src/test/statsgenotype/validation/muts/" + mutPrefix +
                      ".fa -n " + numReads + " --output ../src/test/statsgenotype/validation/reads/" +
                      mutPrefix + " --cpus 4";
    //cout << readsCmd << endl;
    if (!filesystem::exists(("../src/test/statsgenotype/validation/reads/" + mutPrefix + "_R1.fastq"))) {
        system(readsCmd.c_str());
    }
    
    // map reads
    string mapCmd = "minimap2 -ax sr " + refPath + " ../src/test/statsgenotype/validation/reads/" +
                    mutPrefix + "_R*.fastq --heap-sort=yes | samtools sort -O sam > " +
                    workingDir + "sam/" + mutPrefix + ".sam";
    // string mapCmd = "bwa mem -t 4 " + refPath + " ../src/test/statsgenotype/validation/reads/" +
    //                 mutPrefix + "_R*.fastq | samtools sort -O sam > " +
    //                 "../src/test/statsgenotype/validation/sam/" + mutPrefix + ".sam";
    // cout << mapCmd << endl;
    if (!filesystem::exists((workingDir + "sam/" + mutPrefix + ".sam"))) {
        system(mapCmd.c_str());
    }

    // pile up
    string pluCmd = "samtools mpileup " + workingDir+ "sam/"  + mutPrefix + ".sam " +
                    "-f " + refPath + " > " + workingDir + "pileup/" + mutPrefix + ".pileup";
    //cout << pluCmd << endl;
    if (!filesystem::exists((workingDir + "pileup/" + mutPrefix + ".pileup"))) {
        system(pluCmd.c_str());
    }
    
    // vcf
    if (!filesystem::exists((workingDir + "vcf/" + mutPrefix + ".vcf"))) {
        ifstream fin((workingDir + "pileup/" + mutPrefix + ".pileup"));
        // ifstream mmin("../mutation_matrices/empty.mutmat");
        ifstream mmin("../mutation_matrices/sars_4k.mutmat");

        streambuf* originalCoutBuffer = std::cout.rdbuf();
        stringstream buffer;
        cout.rdbuf(buffer.rdbuf());

        T->printSamplePlacementVCF(fin, &mmin);
        
        cout.rdbuf(originalCoutBuffer);
        string capturedVCF = buffer.str();

        // save vcf to file
        ofstream vcfout((workingDir + "vcf/" + mutPrefix + ".vcf"));
        vcfout << capturedVCF;

        fin.close();
        mmin.close();
        vcfout.close();
    }

    // compare
    ifstream simIn((workingDir + "muts/" + mutPrefix + ".rf"));
    ifstream vcfIn((workingDir + "vcf/"  + mutPrefix + ".vcf"));
    ifstream pluIn((workingDir + "pileup/"  + mutPrefix + ".pileup"));
    map<size_t, pair<string, string> > simVars = parseSim(simIn);
    map<size_t, pair< vector<string>, vector<int> > > vcfVars = parseVcf(vcfIn);
    vector<size_t> depths = parseDepths(pluIn, refSeq.size());

    simIn.close();
    vcfIn.close();

    int truePositive  = 0;
    int falsePositive = 0;
    int trueNegative  = 0;
    int falseNegative = 0;
    unordered_set<size_t> visited;

    for (const auto& simVar : simVars) {
        auto pos = simVar.first;
        if (isVisited(visited, pos)) {
            continue;
        } else {
            visited.insert(pos);
        }

        bool isIndel = true;
        if (simVar.second.first.size() == simVar.second.second.size()) {
            isIndel = false;
        }

        if (byBaseCoverage.find(depths[pos]) == byBaseCoverage.end()) {
            byBaseCoverage[depths[pos]].resize(4);
        }

        if (vcfVars.find(pos) != vcfVars.end()) {
            /*
            Same position. Check if variants match as well.
            */
            size_t gtIdx;
            for (size_t i = 0; i < vcfVars[pos].second.size(); i++) {
                if (vcfVars[pos].second[i] == 0) {
                    gtIdx = i;
                    break;
                }
            }

            if ((simVar.second.first == vcfVars[pos].first[0]) && (simVar.second.second == vcfVars[pos].first[gtIdx])) {
                truePositive++;
                byBaseCoverage[depths[pos]][0]++;
            } else {
                falseNegative++;
                byBaseCoverage[depths[pos]][1]++;
                byBaseCoverage[depths[pos]][3]++;
            }
        } else if (isIndel) {
            /*
            if no matching position found and is indel, check for near indels that give same alignment.
            */
            bool matchFound = false;
            for (const auto& vcfVar : vcfVars) {
                if (isVisited(visited, vcfVar.first)) {
                    continue;
                }
                if (isoAlignment(simVar, vcfVar, refSeq)) {
                    truePositive++;
                    if (byBaseCoverage.find(depths[vcfVar.first]) == byBaseCoverage.end()) {
                        byBaseCoverage[depths[vcfVar.first]].resize(4);
                    }
                    byBaseCoverage[depths[vcfVar.first]][0]++;
                    visited.insert(vcfVar.first);
                    matchFound = true;
                }
            }
            if (!matchFound) {
                falseNegative++;
                byBaseCoverage[depths[pos]][3]++;
            }
        } else {
            falseNegative++;
            byBaseCoverage[depths[pos]][3]++;
        }
        
    }

    for (const auto& vcfVar : vcfVars) {
        auto pos = vcfVar.first;
        
        if (isVisited(visited, pos)) {
            continue;
        }
        if (byBaseCoverage.find(depths[pos]) == byBaseCoverage.end()) {
            byBaseCoverage[depths[pos]].resize(4);
        }

        if (vcfVar.second.second[0] == 0) {
            trueNegative++;
            byBaseCoverage[depths[pos]][2]++;
        } else {
            falsePositive++;
            byBaseCoverage[depths[pos]][1]++;
        }
        
        visited.insert(pos);
    }

    
    for(size_t i = 1; i < refSeq.size() + 1; i++) {
        if (visited.find(i) == visited.end()) {
            if (byBaseCoverage.find(depths[i]) == byBaseCoverage.end()) {
                byBaseCoverage[depths[i]].resize(4);
            }
            byBaseCoverage[depths[i]][2]++;
        }
    }

    TestStats stats(truePositive, falsePositive, trueNegative, falseNegative, refSeq.size());
    return stats;
}

int distanceToRoot(Tree* T, Node* node) {
    int distance = 0;
    while (true) {
        // cout << node->identifier << "\t";
        auto edges = getEdgeCoor(T, node->identifier);
        int masked = maskedSize(T, node->nucMutation, edges);
        // cout << masked << endl;
        distance += masked;
        if (node->parent->identifier == T->root->identifier) {
            break;
        }
        node = node->parent;
    }
    return distance;
}

/*
Things to test for:
    - variant number         1 - 100; SNP only, INDEL only, mix
    - variant type           SNP vs INS vs DEL
    - variant length         1 - 20
    - read depth             1x - 10x; 100 -> 1000
*/

BOOST_AUTO_TEST_CASE(validationTesting) {
    using namespace std;
    std::ifstream is("/home/azhang/rotations/rotation_2/pangenome-mat/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    Tree *T = new Tree(inputStream);

    std::ifstream mmi("../mutation_matrices/empty.mutmat");
    auto mutmat = statsgenotype::mutationMatrices();
    int idx = 0;
    string line;
    while(getline(mmi, line)) {
        vector<double> probs;
        vector<string> fields;
        stringSplit(line, ' ', fields);
        for (const auto& f : fields) {
            probs.push_back(stod(f));
        }
        if (idx < 4) {
            mutmat.submat[idx] = move(probs);
        } else if (idx == 4) {
            mutmat.insmat = move(probs);
        } else {
            mutmat.delmat = move(probs);
        }
        idx++;
    }



    /*
    vector< pair<string, string> > refsamples = {
        make_pair("ON824526.1", "20"), make_pair("ON836678.1", "20"), make_pair("ON852975.1", "20"), make_pair("ON875122.1", "20"), make_pair("ON823489.1", "20"),
        make_pair("ON823489.1", "40"), make_pair("ON870142.1", "40"), make_pair("ON869493.1", "40"), make_pair("ON870165.1", "40"), make_pair("ON870142.1", "40"),
        make_pair("ON833982.1", "60"), make_pair("ON869964.1", "60"), make_pair("ON811181.1", "60"), make_pair("ON845091.1", "60"), make_pair("ON872856.1", "60"),
        make_pair("LC603288.1", "80"), make_pair("ON876177.1", "80"), make_pair("OM188340.1", "80"), make_pair("ON831637.1", "80"), make_pair("ON870955.1", "80"),
        make_pair("ON833626.1", "100"), make_pair("ON872899.1", "100"), make_pair("ON833710.1", "100"), make_pair("ON832661.1", "100"), make_pair("ON833719.1", "100")
    };

    for (const auto& refsample : refsamples) {
        string reference = refsample.first;
        string referencePath = "/home/azhang/rotations/rotation_2/panmap/talk/20ksamples/fasta/" + reference + ".fa";
        string rootPath = "/home/azhang/rotations/rotation_2/panmap/talk/20kroot.fa";
        string refSeq = T->getStringFromReference(reference, false);
        string rootSeq = T->getStringFromReference(T->root->identifier, false);
        
        //cout << "20k tree root: " << T->root->identifier << endl;
        //makeFasta(T, T->root->identifier, "/home/azhang/rotations/rotation_2/panmap/talk/20ksamples/root.fa");
        filesystem::create_directories("../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference);
        filesystem::create_directories("../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference + "/pileup");
        filesystem::create_directories("../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference + "/sam");
        filesystem::create_directories("../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference + "/vcf");
        map<size_t, vector<int> > byBaseCoverageTest;
        vector<size_t> depths = {1, 5, 10, 15, 20};
        // for (size_t depth : depths) {
        //     for (size_t it = 0; it < 10; it++) {
        //         string prefix = reference + "_prior_20_1_1_" + to_string(depth) + "_" + to_string(it);
        //         genMut2(refSeq, mutmat, prefix, 200, 200, 20, 1, {1});
        //         auto testStats  = runTest(T, refSeq,  prefix, referencePath, depth, "../src/test/statsgenotype/validation/", byBaseCoverageTest);
        //         auto testStats2 = runTest(T, rootSeq, prefix, rootPath,      depth, "../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference + "/", byBaseCoverageTest);
        //     }
        // }
        byBaseCoverageTest.clear();
        for (size_t depth : depths) {
            for (size_t it = 0; it < 10; it++) {
                string prefix = reference + "_20_9_6_" + to_string(depth) + "_" + to_string(it);
                genMut2(refSeq, mutmat, prefix, 200, 200, 20, 9, {2, 4});
                auto testStats  = runTest(T, refSeq,  prefix, referencePath, depth, "../src/test/statsgenotype/validation/", byBaseCoverageTest);
                auto testStats2 = runTest(T, rootSeq, prefix, rootPath,      depth, "../src/test/statsgenotype/validation/to_root/" + refsample.second + "/" + reference + "/", byBaseCoverageTest);
            }
        }
    }

    */
    /*
    depth of coverage
    distance from root
    map to root vs map to reference
    */

    
    // SUB only
    /*
    {
    string reference = "OM878109.1";
    string referencePath = "../src/test/statsgenotype/validation/OM878109.1.fa";
    string refSeq = T->getStringFromReference(reference, false);
    // makeFasta(T, reference, referencePath);
    map<size_t, vector<int> > byBaseCoverageSub;
    for (size_t varLen = 1; varLen < 2; varLen++) {
        vector< tuple<double, double, double, double, double, double> > byCoverage;
        size_t numIt = 20;
        for (size_t coverage = 1; coverage <= 10; coverage++) {
            vector<TestStats> iterations;
            for (size_t it = 0; it < numIt; it++) {
                string mutPrefix = "10_" + to_string(varLen) + "_1_" + to_string(coverage) + "_" + to_string(it);
                genMut(refSeq, mutmat, mutPrefix, 200, 200, 10, varLen, {1});
                auto curRunStats = runTest(T, refSeq, mutPrefix, referencePath, coverage, "../src/test/statsgenotype/validation/", byBaseCoverageSub);
                iterations.push_back(curRunStats);
            }
            
            int sumTP = 0;
            int sumFP = 0;
            int sumTN = 0;
            int sumFN = 0;
            for (const auto& it : iterations) {
                sumTP += it.truePositive;
                sumFP += it.falsePositive;
                sumTN += it.trueNegative;
                sumFN += it.falseNegative;
            }
            double avgTP = static_cast<double>(sumTP) / iterations.size();
            double avgFP = static_cast<double>(sumFP) / iterations.size();
            double avgTN = static_cast<double>(sumTN) / iterations.size();
            double avgFN = static_cast<double>(sumFN) / iterations.size();
            double avgTPR = avgTP / (avgTP + avgFN);
            double avgFPR = avgFP / (avgFP + avgTN);
            
            byCoverage.emplace_back(make_tuple(avgTP, avgFP, avgTN, avgFN, avgTPR, avgFPR));
        }

        ofstream osf("../src/test/statsgenotype/validation/statsout/1_" + to_string(varLen) + "_" + to_string(numIt) + ".out");
        osf << "@it="  << numIt  << "\n";
        osf << "@len=" << varLen << "\n";
        osf << "coverage\tTP\tFP\tTN\tFN\tTPR\tFPR" << endl;
        for(size_t i = 0; i < byCoverage.size(); i++) {
            osf << i + 1 << "\t" << get<0>(byCoverage[i]) << "\t"
                << get<1>(byCoverage[i]) << "\t"
                << get<2>(byCoverage[i]) << "\t"
                << get<3>(byCoverage[i]) << "\t"
                << get<4>(byCoverage[i]) << "\t"
                << get<5>(byCoverage[i]) << "\n";
        }
        osf.close();
    }

    ofstream bcof("../src/test/statsgenotype/validation/statsout/sub_by_coverage.out");
    for (const auto& stats : byBaseCoverageSub) {
        bcof << stats.first << "\t"
             << stats.second[0] << "\t"
             << stats.second[1] << "\t"
             << stats.second[2] << "\t"
             << stats.second[3] << "\t";
        if (stats.second[0] == 0) {
            bcof << 0 << "\n";
        } else if ((stats.second[0] + stats.second[3]) == 0) {
            bcof << 1 << "\n";
        } else {
            bcof << static_cast<double>(stats.second[0]) / (stats.second[0] + stats.second[3]) << "\n";
        }
    }
    bcof.close();
    }
    */
    
    /*
    // INS only
    {
    map<size_t, vector<int> > byBaseCoverageIns;
    for (size_t varLen = 1; varLen <= 15; varLen += 2) {
        vector< tuple<double, double, double, double, double, double> > byCoverage;
        size_t numIt = 20;
        for (size_t coverage = 1; coverage <= 10; coverage++) {
            vector<TestStats> iterations;
            for (size_t it = 0; it < numIt; it++) {
                string mutPrefix = "10_" + to_string(varLen) + "_2_" + to_string(coverage) + "_" + to_string(it);
                genMut(refSeq, mutmat, mutPrefix, 200, 200, 10, varLen, {2});
                auto curRunStats = runTest(T, refSeq, mutPrefix, referencePath, coverage, "../src/test/statsgenotype/validation/", byBaseCoverageIns);
                //updateBaseCoverageStats(mutPrefix, byBaseCoverageSub);
                iterations.push_back(curRunStats);
            }
            
            int sumTP = 0;
            int sumFP = 0;
            int sumTN = 0;
            int sumFN = 0;
            for (const auto& it : iterations) {
                sumTP += it.truePositive;
                sumFP += it.falsePositive;
                sumTN += it.trueNegative;
                sumFN += it.falseNegative;
            }
            double avgTP = static_cast<double>(sumTP) / iterations.size();
            double avgFP = static_cast<double>(sumFP) / iterations.size();
            double avgTN = static_cast<double>(sumTN) / iterations.size();
            double avgFN = static_cast<double>(sumFN) / iterations.size();
            double avgTPR = avgTP / (avgTP + avgFN);
            double avgFPR = avgFP / (avgFP + avgTN);
            
            byCoverage.emplace_back(make_tuple(avgTP, avgFP, avgTN, avgFN, avgTPR, avgFPR));
        }

        ofstream osf("../src/test/statsgenotype/validation/statsout/2_" + to_string(varLen) + "_" + to_string(numIt) + ".out");
        osf << "@it="  << numIt  << "\n";
        osf << "@len=" << varLen << "\n";
        osf << "coverage\tTP\tFP\tTN\tFN\tTPR\tFPR" << endl;
        for(size_t i = 0; i < byCoverage.size(); i++) {
            osf << i + 1 << "\t" << get<0>(byCoverage[i]) << "\t"
                << get<1>(byCoverage[i]) << "\t"
                << get<2>(byCoverage[i]) << "\t"
                << get<3>(byCoverage[i]) << "\t"
                << get<4>(byCoverage[i]) << "\t"
                << get<5>(byCoverage[i]) << "\n";
        }
        osf.close();
    }
    ofstream bcof("../src/test/statsgenotype/validation/statsout/ins_by_coverage.out");
    for (const auto& stats : byBaseCoverageIns) {
        bcof << stats.first << "\t"
             << stats.second[0] << "\t"
             << stats.second[1] << "\t"
             << stats.second[2] << "\t"
             << stats.second[3] << "\t";
        if (stats.second[0] == 0) {
            bcof << 0 << "\n";
        } else if ((stats.second[0] + stats.second[3]) == 0) {
            bcof << 1 << "\n";
        } else {
            bcof << static_cast<double>(stats.second[0]) / (stats.second[0] + stats.second[3]) << "\n";
        }
    }
    bcof.close();
    }

    // DEl only

    {
    map<size_t, vector<int> > byBaseCoverageDel;
    for (size_t varLen = 1; varLen <= 15; varLen += 2) {
        vector< tuple<double, double, double, double, double, double> > byCoverage;
        size_t numIt = 20;
        for (size_t coverage = 1; coverage <= 10; coverage++) {
            vector<TestStats> iterations;
            for (size_t it = 0; it < numIt; it++) {
                string mutPrefix = "10_" + to_string(varLen) + "_4_" + to_string(coverage) + "_" + to_string(it);
                genMut(refSeq, mutmat, mutPrefix, 200, 200, 10, varLen, {4});
                auto curRunStats = runTest(T, refSeq, mutPrefix, referencePath, coverage, "../src/test/statsgenotype/validation/", byBaseCoverageDel);
                iterations.push_back(curRunStats);
            }
            
            int sumTP = 0;
            int sumFP = 0;
            int sumTN = 0;
            int sumFN = 0;
            for (const auto& it : iterations) {
                sumTP += it.truePositive;
                sumFP += it.falsePositive;
                sumTN += it.trueNegative;
                sumFN += it.falseNegative;
            }
            double avgTP = static_cast<double>(sumTP) / iterations.size();
            double avgFP = static_cast<double>(sumFP) / iterations.size();
            double avgTN = static_cast<double>(sumTN) / iterations.size();
            double avgFN = static_cast<double>(sumFN) / iterations.size();
            double avgTPR = avgTP / (avgTP + avgFN);
            double avgFPR = avgFP / (avgFP + avgTN);
            
            byCoverage.emplace_back(make_tuple(avgTP, avgFP, avgTN, avgFN, avgTPR, avgFPR));
        }

        ofstream osf("../src/test/statsgenotype/validation/statsout/4_" + to_string(varLen) + "_" + to_string(numIt) + ".out");
        osf << "@it="  << numIt  << "\n";
        osf << "@len=" << varLen << "\n";
        osf << "coverage\tTP\tFP\tTN\tFN\tTPR\tFPR" << endl;
        for(size_t i = 0; i < byCoverage.size(); i++) {
            osf << i + 1 << "\t" << get<0>(byCoverage[i]) << "\t"
                << get<1>(byCoverage[i]) << "\t"
                << get<2>(byCoverage[i]) << "\t"
                << get<3>(byCoverage[i]) << "\t"
                << get<4>(byCoverage[i]) << "\t"
                << get<5>(byCoverage[i]) << "\n";
        }
        osf.close();
    }
    ofstream bcof("../src/test/statsgenotype/validation/statsout/del_by_coverage.out");
    for (const auto& stats : byBaseCoverageDel) {
        bcof << stats.first << "\t"
             << stats.second[0] << "\t"
             << stats.second[1] << "\t"
             << stats.second[2] << "\t"
             << stats.second[3] << "\t";
        if (stats.second[0] == 0) {
            bcof << 0 << "\n";
        } else if ((stats.second[0] + stats.second[3]) == 0) {
            bcof << 1 << "\n";
        } else {
            bcof << static_cast<double>(stats.second[0]) / (stats.second[0] + stats.second[3]) << "\n";
        }
    }
    bcof.close();
    }
    */
}

BOOST_AUTO_TEST_CASE(genReadsForTesting) {
    using namespace std;
    std::ifstream is("/home/azhang/rotations/rotation_2/pangenome-mat/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    Tree *T = new Tree(inputStream);

    std::ifstream mmi("../mutation_matrices/sars_4k.mutmat");
    auto mutmat = statsgenotype::mutationMatrices();
    int idx = 0;
    string line;
    while(getline(mmi, line)) {
        vector<double> probs;
        vector<string> fields;
        stringSplit(line, ' ', fields);
        for (const auto& f : fields) {
            probs.push_back(stod(f));
        }
        if (idx < 4) {
            mutmat.submat[idx] = move(probs);
        } else if (idx == 4) {
            mutmat.insmat = move(probs);
        } else {
            mutmat.delmat = move(probs);
        }
        idx++;
    }
    string reference = "node_3";
    string refSeq = T->getStringFromReference(reference, false);
    makeFasta(T, reference, "/home/azhang/rotations/rotation_2/panmap/src/test/data/sim_variants/node_3.fa");
    for (size_t i = 0; i < 100; i++) {
        string prefix = "node_3_variant_" + to_string(i);
        string outDir = "/home/azhang/rotations/rotation_2/panmap/src/test/data/sim_variants/";
        genMut2(refSeq, mutmat, prefix, 100, 100, 15, 9,
            {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4}, outDir);
        string faFile = outDir + prefix + ".fa";
        string rfFile = outDir + prefix + ".rf";
        fs::rename(fs::path(faFile), fs::path(outDir + "fasta/" + prefix + ".fa"));
        fs::rename(fs::path(rfFile), fs::path(outDir + "vcftrue/" + prefix + ".rf"));
        simReads(outDir + "fasta/" + prefix + ".fa", outDir + "reads/" + prefix + "_reads", 2000);

    }
    
}


// BOOST_AUTO_TEST_CASE(_indexSyncmers) {
//     size_t k = 16;
//     size_t s = 5;
//     std::ifstream is("../mammal_mito_refseq.pmat");
// //    std::string fastqPath = "";
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSyncmers(T, os, k, s);
    
//     is.close();
//     os.close();

    // PangenomeMAT::Node *root = T->root;
    // std::vector<read_t> reads;

    // struct seedIndex index;
    // std::ifstream indexFile("./test.out");

    // PangenomeMAT::loadIndex(T->root, indexFile, index);

    // auto fastq_start = std::chrono::high_resolution_clock::now();
    // std::set<kmer_t> readSyncmers = syncmersFromFastq(fastqPath, reads, k, s);
    // auto fastq_end = std::chrono::high_resolution_clock::now();

    // std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

    // auto place_start = std::chrono::high_resolution_clock::now();

    // std::set<kmer_t> rootSyncmers = std::set<kmer_t>(index.rootSeeds.begin(), index.rootSeeds.end());

    // std::cerr << "\n";
    // std::cerr << "Placing sample...\n";


    // struct dynamicJaccard dj;
 
    // dj.intersectionSize = intersection_size(rootSyncmers, readSyncmers);
    // dj.unionSize = rootSyncmers.size() + readSyncmers.size() - dj.intersectionSize;
    // dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
    
    
    // std::cout << "root seeds: " << rootSyncmers.size() << "\n";
    // std::cout << "read seeds: " << readSyncmers.size() << "\n";
    // for (const auto &k : readSyncmers) {
    //     std::cout << k.seq << "\n";
    // }
    // std::cout << "initial jaccard: " << dj.jaccardIndex << "\n";

    // std::unordered_map<std::string, float> scores;
    // std::unordered_map<std::string, bool> readSyncmersMap;
    // for (const auto &k : readSyncmers) {
    //     readSyncmersMap[k.seq] = true;
    // }

    // placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);

    // auto place_end = std::chrono::high_resolution_clock::now();

    // std::cout << "place time: " << std::chrono::duration_cast<std::chrono::milliseconds>(place_end - place_start).count() << "\n";

    // std::vector<std::pair<std::string, float>> v;
    // for ( const auto &p : scores ) {
    //     v.push_back(std::make_pair(p.first, p.second));
    // } 
    // std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
    //     return left.second > right.second;
    // });

    // std::string best_match = v[0].first;
    // for (const auto &s : v) {
    //     std::cerr << s.first << ": " << s.second << "\n";
    // }

// }

