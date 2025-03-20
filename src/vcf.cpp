#include "panmanUtils.hpp"

void panmanUtils::Tree::printVCFParallel(std::string reference, std::ostream& fout) {

    std::string referenceSequence = getStringFromReference(reference);

    if(referenceSequence == "Error: Reference sequence with matching name not found!") {
        std::cerr << referenceSequence << std::endl;
        return;
    }

    size_t recordID = 0;

    std::mutex vcfMapMutex;
    std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

    tbb::parallel_for_each(allNodes, [&](auto& n) {
        if(n.second->children.size() == 0 && n.first != reference) {
            std::string altSequence = getStringFromReference(n.first);
            if(altSequence.length() != referenceSequence.length()) {
                std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
                return;
            }

            std::string currentRefString, currentAltString;
            int currentCoordinate = 1;

            int diffStart = 1;

            for(size_t i = 0; i < referenceSequence.length(); i++) {

                if(referenceSequence[i] == '-' && altSequence[i] == '-') {
                    continue;
                } else if(referenceSequence[i] != '-' && altSequence[i] == '-') {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }

                    currentRefString += referenceSequence[i];
                } else if(referenceSequence[i] == '-' && altSequence[i] != '-') {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }

                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] != altSequence[i]) {
                    if(currentRefString == "" && currentAltString == "") {
                        diffStart = currentCoordinate;
                    }
                    if(currentRefString == currentAltString) {
                        currentRefString = "";
                        currentAltString = "";
                        diffStart = currentCoordinate;
                    }
                    currentRefString += referenceSequence[i];
                    currentAltString += altSequence[i];
                } else if(referenceSequence[i] == altSequence[i]) {
                    if(currentRefString == currentAltString) {
                        // Reset
                        diffStart = currentCoordinate;
                        currentRefString = "";
                        currentRefString += referenceSequence[i];
                        currentAltString = currentRefString;
                    } else {
                        // Create VCF record at position i
                        if(currentRefString == "") {
                            currentRefString += referenceSequence[i];
                            currentAltString += altSequence[i];
                            diffStart = currentCoordinate;
                            vcfMapMutex.lock();
                            vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                            vcfMapMutex.unlock();
                            diffStart = currentCoordinate+1;
                            currentRefString = "";
                            currentAltString = "";
                        } else {
                            vcfMapMutex.lock();
                            vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                            vcfMapMutex.unlock();

                            // Reset
                            diffStart = currentCoordinate;
                            currentRefString = "";
                            currentRefString += referenceSequence[i];
                            currentAltString = currentRefString;
                        }
                    }
                }

                if(referenceSequence[i] != '-') {
                    currentCoordinate++;
                }
            }

            if(currentRefString != currentAltString) {
                vcfMapMutex.lock();
                vcfMap[diffStart][currentRefString][currentAltString].push_back(n.first);
                vcfMapMutex.unlock();

                // Reset
                diffStart = referenceSequence.size();
                currentRefString = "";
                currentAltString = currentRefString;
            }
        }
    });

    std::cout << vcfMap.size() << std::endl;

    std::mutex sequenceIdsMutex;
    std::map< std::string, size_t > sequenceIds;
    tbb::parallel_for_each(allNodes, [&](auto& u) {
        if(u.second->children.size() == 0 && u.first != reference) {
            sequenceIdsMutex.lock();
            sequenceIds[u.first] = 0;
            sequenceIdsMutex.unlock();
        }
    });


    fout << "##fileformat=VCFv" << VCF_VERSION << '\n';
    fout << "##fileDate=" << panmanUtils::getDate() << '\n';
    fout << "##source=PanMATv" << PMAT_VERSION << '\n';
    fout << "##reference=" << reference << '\n';
    fout << "#CHROM\t" << "POS\t" << "ID\t" << "REF\t" << "ALT\t" << "QUAL\t" << "FILTER\t" << "INFO\t" << "FORMAT\t";

    // fout << std::left << std::setw(20) << "#CHROM " << std::setw(20) << "POS " << std::setw(20) << "ID " << std::setw(20) << "REF " << std::setw(20) << "ALT " << std::setw(20) << "QUAL " << std::setw(20) << "FILTER " << std::setw(20) << "INFO " << std::setw(20) << "FORMAT ";
    for(auto u: sequenceIds) {
        if(u.first != sequenceIds.rbegin()->first) {
            fout << u.first + "\t";
        } else {
            fout << u.first;
        }
    }
    fout << '\n';

    for(auto u: vcfMap) {
        for(auto v: u.second) {
            if(v.first == "") {
                fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << ".\t";
            } else {
                fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << v.first << "\t";
            }

            std::map< std::string, size_t > tempSequenceIds = sequenceIds;

            int ctr = 1;
            std::string altStrings;

            for(auto w: v.second) {
                altStrings += (w.first == "" ? ".": w.first);
                altStrings += ",";
                for(auto uu: w.second) {
                    tempSequenceIds[uu] = ctr;
                }
                ctr++;
            }

            altStrings.pop_back();

            fout << altStrings << "\t.\t.\t.\t.\t";

            for(auto w: tempSequenceIds) {
                if(w.first != sequenceIds.rbegin()->first) {
                    fout << w.second << "\t";
                } else {
                    fout << w.second;
                }
            }

            fout << '\n';
        }
    }
}


void panmanUtils::Tree::printVCFParallel(panmanUtils::Node* refnode, std::string& fileName) {

    std::string reference = refnode->identifier;
    std::pair<std::vector<std::string>, std::vector<int>> referenceSequence = extractSingleSequence(refnode, true);

    if(reference == "") {
        std::cerr << "Reference not set correctly" << std::endl;
        return;
    }

    // std::string fileName="";
    
    size_t itr = 0;
    size_t batch_size = 1000;
    std::mutex sequenceIdsMutex;
    std::map< std::string, size_t > sequenceIds;
    tbb::parallel_for_each(allNodes, [&](auto& u) {
        if(u.second->children.size() == 0 && u.first != reference) {
            sequenceIdsMutex.lock();
            sequenceIds[u.first] = 0;
            sequenceIdsMutex.unlock();
        }
    });

    std::vector<panmanUtils::Node *> allNodes_vector;
    for(auto u: allNodes) {
        if (u.second->children.size() == 0 && u.second->identifier != reference)
            allNodes_vector.push_back(u.second);
    }

    for (size_t nc=0; nc<allNodes_vector.size(); nc+=batch_size)
    {
        size_t recordID = 0;
        std::ofstream outputFile;
        std::streambuf * buf;
        if(fileName != "") {
            outputFile.open("./info/" + fileName + "_" + std::to_string(itr++) + ".vcf");
            buf = outputFile.rdbuf();
        } else {
            buf = std::cout.rdbuf();
        }
        std::ostream fout (buf);
        std::mutex vcfMapMutex;
        std::map< int, std::map< std::string, std::map< std::string, std::vector< std::string > > > > vcfMap;

        fout << "##fileformat=VCFv" << VCF_VERSION << '\n';
        fout << "##fileDate=" << panmanUtils::getDate() << '\n';
        fout << "##source=PanMATv" << PMAT_VERSION << '\n';
        fout << "##reference=" << reference << '\n';
        fout << "#CHROM\t" << "POS\t" << "ID\t" << "REF\t" << "ALT\t" << "QUAL\t" << "FILTER\t" << "INFO\t" << "FORMAT\t";


        // fout << std::left << std::setw(20) << "#CHROM " << std::setw(20) << "POS " << std::setw(20) << "ID " << std::setw(20) << "REF " << std::setw(20) << "ALT " << std::setw(20) << "QUAL " << std::setw(20) << "FILTER " << std::setw(20) << "INFO " << std::setw(20) << "FORMAT ";
        for (size_t nj=nc; nj<std::min(nc+batch_size, allNodes_vector.size()); nj++)
        {
            panmanUtils::Node* n = allNodes_vector[nj];
            if(n->identifier != sequenceIds.rbegin()->first) {
                fout << n->identifier + "\t";
            } else {
                fout << n->identifier;
            }
        }
        fout << '\n';
        tbb::parallel_for(nc, std::min(nc+batch_size, allNodes_vector.size()), [&](size_t nj) {
            panmanUtils::Node* n = allNodes_vector[nj];
        // tbb::parallel_for_each(allNodes, [&](auto& n) {
            if(n->children.size() == 0 && n->identifier != refnode->identifier) {
                std::pair<std::vector<std::string>, std::vector<int>> altSequence = extractSingleSequence(n, true);
                size_t altTotalLen = 0;
                size_t refTotalLen = 0;
                for(auto u: referenceSequence.second) {
                    refTotalLen += u;
                }
                for(auto u: altSequence.second) {
                    altTotalLen += u;
                }
                if (altTotalLen != refTotalLen) {
                    std::cerr << "Logic error. String lengths don't match: " << refTotalLen << " " << altTotalLen << std::endl;
                    return;
                }
                // if(altSequence.length() != referenceSequence.length()) {
                //     std::cerr << "Logic error. String lengths don't match: " << referenceSequence.length() << " " << altSequence.length() << std::endl;
                //     return;
                // }

                std::string currentRefString, currentAltString;
                int currentCoordinate = 1;

                int diffStart = 1;

                std::vector<std::string> referenceString = std::get<0>(referenceSequence);
                std::vector<int> referenceLen = std::get<1>(referenceSequence);
                std::vector<std::string> altString = std::get<0>(altSequence);
                std::vector<int> altLen = std::get<1>(altSequence);
                
                for (int b=0; b < referenceLen.size();b++) {
                    std::string referenceChar = "";
                    std::string altChar = "";
                    
                    if (referenceString[b] == "-") {
                        for (int c=0; c<referenceLen[b];c++) {
                            referenceChar += "-";
                        }
                    } else referenceChar += referenceString[b];

                    if (altString[b] == "-") {
                        for (int c=0; c<altLen[b];c++) {
                            altChar += "-";
                        }
                    } else altChar += altString[b];

                    for(size_t i = 0; i < referenceChar.length(); i++) {

                        if(referenceChar[i] == '-' && altChar[i] == '-') {
                            continue;
                        } else if(referenceChar[i] != '-' && altChar[i] == '-') {
                            if(currentRefString == "" && currentAltString == "") {
                                diffStart = currentCoordinate;
                            }

                            currentRefString += referenceChar[i];
                        } else if(referenceChar[i] == '-' && altChar[i] != '-') {
                            if(currentRefString == "" && currentAltString == "") {
                                diffStart = currentCoordinate;
                            }

                            currentAltString += altChar[i];
                        } else if(referenceChar[i] != altChar[i]) {
                            if(currentRefString == "" && currentAltString == "") {
                                diffStart = currentCoordinate;
                            }
                            if(currentRefString == currentAltString) {
                                currentRefString = "";
                                currentAltString = "";
                                diffStart = currentCoordinate;
                            }
                            currentRefString += referenceChar[i];
                            currentAltString += altChar[i];
                        } else if(referenceChar[i] == altChar[i]) {
                            if(currentRefString == currentAltString) {
                                // Reset
                                diffStart = currentCoordinate;
                                currentRefString = "";
                                currentRefString += referenceChar[i];
                                currentAltString = currentRefString;
                            } else {
                                // Create VCF record at position i
                                if(currentRefString == "") {
                                    currentRefString += referenceChar[i];
                                    currentAltString += altChar[i];
                                    diffStart = currentCoordinate;
                                    vcfMapMutex.lock();
                                    vcfMap[diffStart][currentRefString][currentAltString].push_back(n->identifier);
                                    vcfMapMutex.unlock();
                                    diffStart = currentCoordinate+1;
                                    currentRefString = "";
                                    currentAltString = "";
                                } else {
                                    vcfMapMutex.lock();
                                    vcfMap[diffStart][currentRefString][currentAltString].push_back(n->identifier);
                                    vcfMapMutex.unlock();

                                    // Reset
                                    diffStart = currentCoordinate;
                                    currentRefString = "";
                                    currentRefString += referenceChar[i];
                                    currentAltString = currentRefString;
                                }
                            }
                        }

                        if(referenceChar[i] != '-') {
                            currentCoordinate++;
                        }
                    }
                }
                if(currentRefString != currentAltString) {
                    vcfMapMutex.lock();
                    vcfMap[diffStart][currentRefString][currentAltString].push_back(n->identifier);
                    vcfMapMutex.unlock();

                    // Reset
                    diffStart = refTotalLen;
                    currentRefString = "";
                    currentAltString = currentRefString;
                }
            }
        });


        

        for(auto u: vcfMap) {
            for(auto v: u.second) {
                if(v.first == "") {
                    fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << ".\t";
                } else {
                    fout << reference << "\t" << u.first << "\t" << recordID++ << "\t" << v.first << "\t";
                }

                std::map< std::string, size_t > tempSequenceIds = sequenceIds;

                int ctr = 1;
                std::string altStrings;

                for(auto w: v.second) {
                    altStrings += (w.first == "" ? ".": w.first);
                    altStrings += ",";
                    for(auto uu: w.second) {
                        tempSequenceIds[uu] = ctr;
                    }
                    ctr++;
                }

                altStrings.pop_back();

                fout << altStrings << "\t.\t.\t.\t.\t";

                for(auto w: tempSequenceIds) {
                    if(w.first != sequenceIds.rbegin()->first) {
                        fout << w.second << "\t";
                    } else {
                        fout << w.second;
                    }
                }

                fout << '\n';
            }
        }
    }

    
}
