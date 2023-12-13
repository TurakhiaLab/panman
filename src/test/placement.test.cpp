 #define BOOST_TEST_MODULE Index construction
#include <boost/test/included/unit_test.hpp>
#include "../PangenomeMAT.hpp"
#include <stack>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <json/json.h>
#include <filesystem>
#include <regex>




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

size_t get_distance(Tree* T, Node* node_1, Node* node_2) {
    size_t n1_dist = 0;
    size_t n2_dist = 0;
    while (node_1->identifier != node_2->identifier) {
        if (node_1->level == node_2->level) {
            n1_dist += node_1->nucMutation.size();
            n2_dist += node_2->nucMutation.size();
            node_1 = node_1->parent;
            node_2 = node_2->parent;
        } else if (node_1->level > node_2->level) {
            n1_dist += node_1->nucMutation.size();
            node_1 = node_1->parent;
        } else if (node_2->level > node_1->level) {
            n2_dist += node_2->nucMutation.size();
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
            n1_dist += (node_1->nucMutation.size() == 0) ? 0 : 1;
            n2_dist += (node_2->nucMutation.size() == 0) ? 0 : 1;
            node_1 = node_1->parent;
            node_2 = node_2->parent;
        } else if (node_1->level > node_2->level) {
            n1_dist += (node_1->nucMutation.size() == 0) ? 0 : 1;
            node_1 = node_1->parent;
        } else if (node_2->level > node_1->level) {
            n2_dist += (node_2->nucMutation.size() == 0) ? 0 : 1;
            node_2 = node_2->parent;
        }
    }
    //cout << n1_dist << "\t" << n2_dist << endl;
    size_t distance = n1_dist + n2_dist;
    return distance;
}

BOOST_AUTO_TEST_CASE(genotypeUncertainties) {
    using namespace std;


    std::ifstream is("../sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    Tree *T = new Tree(inputStream);

    vector< pair<string, string> > samples = get_pruned_samples("../src/test/statsgenotype/pmat/");
    for (const auto& sample : samples) {
        // get pruned sample name, tree path, and fastq files
        string pruned_sample = sample.first;
        string pruned_tree_path = sample.second;
        string pruned_sample_fastq_1 = "../src/test/statsgenotype/fastq_2k/" + pruned_sample + ".1kreads.fastq_R1.fastq";
        string pruned_sample_fastq_2 = "../src/test/statsgenotype/fastq_2k/" + pruned_sample + ".1kreads.fastq_R2.fastq";
        cout << pruned_sample << "\t" << pruned_tree_path << "\t" << pruned_sample_fastq_1 << pruned_sample_fastq_2 << endl;

        // read pruned tree
        ifstream ptis(pruned_tree_path);
        boost::iostreams::filtering_streambuf< boost::iostreams::input> inptPMATBuffer;
        inptPMATBuffer.push(boost::iostreams::gzip_decompressor());
        inptPMATBuffer.push(ptis);
        istream ptinputStream(&inptPMATBuffer);
        Tree *pT = new Tree(ptinputStream);
        ptis.close();

        // index pruned tree
        size_t k = 13;
        size_t s = 7;
        string pruned_tree_index_path = "../src/test/statsgenotype/indices/" + pruned_sample + ".index.out";
        if (!filesystem::exists(pruned_tree_index_path)) {
            ofstream ptos(pruned_tree_index_path);
            indexSyncmers(pT, ptos, k, s);
            ptos.close();
        }
        
        
        struct seedIndex index;
        std::ifstream indexFile(pruned_tree_index_path);
        PangenomeMAT::loadIndex(pT->root, indexFile, index);

        // place index and print SAM

        // PangenomeMAT::placeSample(pT, pruned_sample_fastq_1, index, k, s);
        // break;


        std::string best_match;
        PangenomeMAT::placeSample(pT, pruned_sample_fastq_1, index, k, s, best_match);
        
        cout << pruned_sample << "\t" << T->allNodes[pruned_sample]->parent->identifier << "\t" << best_match << endl;
        
        const std::string ref_node_seq = T->getStringFromReference(best_match, false);
        std::string reference_node_fa_path = "../src/test/statsgenotype/fasta/" + best_match + ".fa";
        if(!filesystem::exists(reference_node_fa_path)) {
            ofstream faos(reference_node_fa_path);
            faos << '>' << best_match << '\n';
            size_t linesize = 80;
            for (size_t i = 0; i < ref_node_seq.size(); i += linesize) {
                faos << ref_node_seq.substr(i, std::min(linesize, ref_node_seq.size() - i)) << '\n';
            }
            faos.close();
        }

        std::string sam_path = "../src/test/statsgenotype/sam/" + pruned_sample + ".sam";
        std::string minimap_cmd_str = "minimap2 -ax sr " + reference_node_fa_path + " " + pruned_sample_fastq_1 + " " + pruned_sample_fastq_2 + " --heap-sort=yes | samtools sort -O sam > " + sam_path;
        cout << minimap_cmd_str << endl;
        if (!filesystem::exists(sam_path)) {
            const char* minimap_cmd = minimap_cmd_str.c_str();
            cout << minimap_cmd_str << endl;
            system(minimap_cmd);
        }

        std::string plu_path = "../src/test/statsgenotype/pileup/" + pruned_sample + ".pileup";
        std::string pileup_cmd_str = "samtools mpileup " + sam_path + " -f " + reference_node_fa_path + " -A > " + plu_path;
        cout << pileup_cmd_str << endl;
        if (!filesystem::exists(plu_path)) {
            const char* pileup_cmd = pileup_cmd_str.c_str();
            cout << pileup_cmd_str << endl;
            system(pileup_cmd);
        }
        
        std::ifstream fin(plu_path);
        std::ifstream min("../mutation_matrices/sars_4k.mutmat");
        pT->printSamplePlacementVCF(fin, &min);

        // cout << endl;
        // cout << "----- Nodes: " << best_match << "\t" << T->allNodes[pruned_sample]->parent->identifier << endl;
        // cout << "----- Level: " << T->allNodes[best_match]->level << "\t" << T->allNodes[pruned_sample]->parent->level << endl;
        // auto dist = get_distance(T, T->allNodes[best_match], T->allNodes[pruned_sample]->parent);
        // cout << "----- Dist:  " << dist << endl;
        // cout << endl;
        // cout << "@"
        //      << pruned_sample << "\t"
        //      << T->allNodes[pruned_sample]->parent->identifier << "\t"
        //      << best_match << "\t"
        //      << get_distance(T, T->allNodes[best_match], T->allNodes[pruned_sample]->parent) << "\t"
        //      << get_distance_branch(T, T->allNodes[best_match], T->allNodes[pruned_sample]->parent) << endl;
        break;
    }
    std::ifstream testplu("../src/test/statsgenotype/pileup/test/test.pileup");
    std::ifstream min("../mutation_matrices/sars_4k.mutmat");
    T->printSamplePlacementVCF(testplu, &min);

    // std::ifstream ecis("../ecoli_100.pmat");
    // boost::iostreams::filtering_streambuf< boost::iostreams::input> ecinPMATBuffer;
    // ecinPMATBuffer.push(boost::iostreams::gzip_decompressor());
    // ecinPMATBuffer.push(ecis);
    // std::istream ecinputStream(&ecinPMATBuffer);
    // Tree *ecT = new Tree(ecinputStream);

    // cout << "NZ_CP007594.1 -> node_63" << endl;
    // Node* sampleNode = ecT->allNodes["NZ_CP007594.1"];
    // cout << sampleNode->parent->identifier << endl;
    // cout << sampleNode->nucMutation.size() << endl;
    // cout << sampleNode->blockMutation.size() << endl;
    // for (const auto& mut : sampleNode->nucMutation) {
    //     cout << mut.nucPosition << "\t" << (mut.mutInfo >> 4) << "\t" << (mut.mutInfo & 15) << endl;
    // }

    is.close();
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

