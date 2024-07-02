#include "panmanUtils.hpp"

std::string getNucleotideSequenceFromBlockCoordinates(std::tuple< int, int, int, int > start,
        std::tuple< int, int, int, int > end, const sequence_t& sequence,
        const blockExists_t& blockExists, const blockStrand_t& blockStrand) {

    std::string ntSequence;

    for(size_t i = (size_t)std::get<0>(start); i <= (size_t)std::get<0>(end); i++) {
        if(blockStrand[i].first) {
            size_t startMainNuc = (i == (size_t)std::get<0>(start))? (size_t)std::get<2>(start): 0;

            for(size_t j = startMainNuc; j < sequence[i].first.size(); j++) {
                int startGapNuc = (i == (size_t)std::get<0>(start) && j == (size_t)startMainNuc)?
                                  std::get<3>(start): 0;
                if(startGapNuc != -1) {
                    for(size_t k = (size_t)std::get<3>(start);
                            k < sequence[i].first[j].second.size(); k++) {
                        if(std::make_tuple((int)i, -1, (int)j, (int)k) == end) {
                            return ntSequence;
                        }
                        if(blockExists[i].first) {
                            ntSequence += sequence[i].first[j].second[k];
                        } else {
                            ntSequence += '-';
                        }
                    }
                }
                if(std::make_tuple((int)i, -1, (int)j, -1) == end) {
                    return ntSequence;
                }
                if(blockExists[i].first) {
                    ntSequence += sequence[i].first[j].first;
                } else {
                    ntSequence += '-';
                }
            }
        } else {
            size_t startMainNuc = (i == (size_t)std::get<0>(start))? (size_t)std::get<2>(start):
                                  sequence[0].first.size()-1;
            for(size_t j = startMainNuc; j+1 > 0; j--) {
                if(std::make_tuple((int)i, -1, (int)j, -1) == end) {
                    return ntSequence;
                }
                if(blockExists[i].first) {
                    ntSequence += sequence[i].first[j].first;
                } else {
                    ntSequence += '-';
                }
                size_t startGapNuc = (i == (size_t)(std::get<0>(start)) && j == startMainNuc) ?
                                     (size_t)std::get<3>(start): sequence[i].first[j].second.size() - 1;
                for(size_t k = startGapNuc; k+1 > 0; k--) {
                    if(std::make_tuple((int)i, -1, (int)j, k) == end) {
                        return ntSequence;
                    }
                    if(blockExists[i].first) {
                        ntSequence += sequence[i].first[j].second[k];
                    } else {
                        ntSequence += '-';
                    }
                }
            }
        }
    }

    return "ERROR";
}

std::tuple< std::vector< std::string >, std::vector< size_t >, std::vector< size_t > >
getAminoAcidSequence(std::tuple< int, int, int, int > start,
                     std::tuple< int, int, int, int > end, const sequence_t& sequence,
                     const blockExists_t& blockExists,
                     const blockStrand_t& blockStrand) {

    std::string ntSequence = getNucleotideSequenceFromBlockCoordinates(start, end, sequence,
                             blockExists, blockStrand);

    if(ntSequence == "ERROR") {
        panmanUtils::printError("Error in translating input coordinates to PanMAT coordinates."
                                " Coordinates may be out of range.");
        return std::make_tuple(std::vector< std::string >({"ERROR"}), std::vector<size_t>(),
                               std::vector<size_t>());
    }

    for(size_t i = 0; i < ntSequence.size(); i++) {
        if(ntSequence[i] != 'A' && ntSequence[i] != 'G' && ntSequence[i] != 'T'
                && ntSequence[i] != 'C') {
            ntSequence[i] = '-';
        }
    }

    std::unordered_map< std::string, std::string > nucToAA = {
        {"TTT", "Phe"},
        {"TTC", "Phe"},
        {"TTA", "Leu"},
        {"TTG", "Leu"},
        {"CTT", "Leu"},
        {"CTC", "Leu"},
        {"CTA", "Leu"},
        {"CTG", "Leu"},
        {"ATT", "Ile"},
        {"ATC", "Ile"},
        {"ATA", "Ile"},
        {"ATG", "Met"},
        {"GTT", "Val"},
        {"GTC", "Val"},
        {"GTA", "Val"},
        {"GTG", "Val"},
        {"TCT", "Ser"},
        {"TCC", "Ser"},
        {"TCA", "Ser"},
        {"TCG", "Ser"},
        {"CCT", "Pro"},
        {"CCC", "Pro"},
        {"CCA", "Pro"},
        {"CCG", "Pro"},
        {"ACT", "Thr"},
        {"ACC", "Thr"},
        {"ACA", "Thr"},
        {"ACG", "Thr"},
        {"GCT", "Ala"},
        {"GCC", "Ala"},
        {"GCA", "Ala"},
        {"GCG", "Ala"},
        {"TAT", "Tyr"},
        {"TAC", "Tyr"},
        {"TAA", "*"},
        {"TAG", "*"},
        {"CAT", "His"},
        {"CAC", "His"},
        {"CAA", "Gln"},
        {"CAG", "Gln"},
        {"AAT", "Asn"},
        {"AAC", "Asn"},
        {"AAA", "Lys"},
        {"AAG", "Lys"},
        {"GAT", "Asp"},
        {"GAC", "Asp"},
        {"GAA", "Glu"},
        {"GAG", "Glu"},
        {"TGT", "Cys"},
        {"TGC", "Cys"},
        {"TGA", "*"},
        {"TGG", "Trp"},
        {"CGT", "Arg"},
        {"CGC", "Arg"},
        {"CGA", "Arg"},
        {"CGG", "Arg"},
        {"AGT", "Ser"},
        {"AGC", "Ser"},
        {"AGA", "Arg"},
        {"AGG", "Arg"},
        {"GGT", "Gly"},
        {"GGC", "Gly"},
        {"GGA", "Gly"},
        {"GGG", "Gly"}
    };

    std::vector< std::string > aaSequence;
    std::vector< size_t > starts, ends;
    std::string currentString;
    for(size_t i = 0; i < ntSequence.size(); i++) {
        if(ntSequence[i] != '-') {
            if(currentString.size() == 0) {
                starts.push_back(i);
            }
            currentString += ntSequence[i];
        }
        if(currentString.size() == 3) {
            ends.push_back(i);
            aaSequence.push_back(nucToAA[currentString]);
            currentString = "";
        }
    }
    while(starts.size() > ends.size()) {
        starts.pop_back();
    }

    assert(starts.size() == ends.size() && starts.size() == aaSequence.size());

    return std::make_tuple(aaSequence, starts, ends);

}

void panmanUtils::Tree::extractAminoAcidTranslations(std::ostream& fout,
        int64_t start, int64_t end) {
    sequence_t referenceSequence;
    blockExists_t referenceBlockExists;
    blockStrand_t referenceBlockStrand;

    if(end <= start) {
        printError("End coordinate must be greater than start");
        return;
    }

    // Get reference sequence in the PanMAT coordinate system
    getSequenceFromReference(referenceSequence, referenceBlockExists, referenceBlockStrand,
                             root->identifier);

    // get PanMAT coordinates from global coordinates
    std::tuple< int, int, int, int > panMATStart = globalCoordinateToBlockCoordinate(start,
            referenceSequence, referenceBlockExists, referenceBlockStrand);
    std::tuple< int, int, int, int > panMATEnd = globalCoordinateToBlockCoordinate(end,
            referenceSequence, referenceBlockExists, referenceBlockStrand);

    if(std::get<0>(panMATStart) == -1 || std::get<0>(panMATEnd) == -1) {
        printError("Error in translating input coordinates to PanMAT coordinates in reference"
                   " sequence. Coordinates may be out of range");
        return;
    }

    auto aaSeq = getAminoAcidSequence(panMATStart, panMATEnd, referenceSequence,
                                      referenceBlockExists, referenceBlockStrand);

    std::cout << std::get<0>(aaSeq).size() << " " << std::get<1>(aaSeq).size() << " " << std::get<2>(aaSeq).size() << std::endl;

    std::vector< std::string > referenceAASequence = std::get<0>(aaSeq);
    std::vector< size_t > referenceStarts = std::get<1>(aaSeq);
    std::vector< size_t > referenceEnds = std::get<2>(aaSeq);

    if(referenceAASequence.size() == 0) {
        return;
    }

    tbb::concurrent_unordered_map< std::string, std::string > aaMutations;

    tbb::parallel_for_each(allNodes, [&](auto u) {
        sequence_t altSequence;
        blockExists_t altBlockExists;
        blockStrand_t altBlockStrand;

        // Get alternate sequence in the PanMAT coordinate system
        getSequenceFromReference(altSequence, altBlockExists, altBlockStrand,
                                 u.first);

        // get PanMAT coordinates from global coordinates
        std::tuple< int, int, int, int > altPanMATStart = globalCoordinateToBlockCoordinate(start,
                altSequence, altBlockExists, altBlockStrand);
        std::tuple< int, int, int, int > altPanMATEnd = globalCoordinateToBlockCoordinate(end,
                altSequence, altBlockExists, altBlockStrand);

        if(std::get<0>(altPanMATStart) == -1 || std::get<0>(altPanMATEnd) == -1) {
            printError("Error in translating input coordinates to PanMAT coordinates in sequence "
                       + u.first + ". Coordinates may be out of range");
            return;
        }

        auto aaSeq = getAminoAcidSequence(altPanMATStart, altPanMATEnd, altSequence,
                                          altBlockExists, altBlockStrand);

        std::vector< std::string > altAASequence = std::get<0>(aaSeq);
        std::vector< size_t > altStarts = std::get<1>(aaSeq);
        std::vector< size_t > altEnds = std::get<2>(aaSeq);

        std::vector< std::pair< size_t, std::string > > matches;
        std::vector< std::pair< size_t, std::string > > insertions;
        std::vector< size_t > deletions;

        size_t refItr = 0, altItr = 0;

        while(altItr < altStarts.size() && refItr < referenceStarts.size()) {
            if(altStarts[altItr] > referenceEnds[refItr]) {
                // deletions
                deletions.push_back(refItr);
                refItr++;
            } else if(altStarts[altItr] < referenceStarts[refItr]) {
                insertions.push_back(std::make_pair(refItr, altAASequence[altItr]));
                altItr++;
            } else {
                // match found
                matches.push_back(std::make_pair(refItr, altAASequence[altItr]));
                altItr++;
                refItr++;
            }
        }
        // Insert remaining Amino Acids
        while(altItr < altStarts.size()) {
            insertions.push_back(std::make_pair(refItr, altAASequence[altItr]));
            altItr++;
        }
        // Delete remaining Amino Acids
        while(refItr < referenceStarts.size()) {
            deletions.push_back(refItr);
            refItr++;
        }

        for(auto i: matches) {
            if(referenceAASequence[i.first] != i.second) {
                aaMutations[u.first] += "S:"+std::to_string(i.first)+":"+i.second+";";
            }
        }
        for(auto i: insertions) {
            aaMutations[u.first] += "I:"+std::to_string(i.first)+":"+i.second+";";
        }
        for(auto i: deletions) {
            aaMutations[u.first] += "D:"+std::to_string(i)+";";
        }
    });

    fout << "node_id\taa_mutations" << "\n";
    for(auto u: aaMutations) {
        fout << u.first << "\t" << u.second << "\n";
    }
}
