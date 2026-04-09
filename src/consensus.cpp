#include "panmanUtils.hpp"

void panmanUtils::Tree::printConsensus(std::ostream& fout){
    std::mutex printMutex;

    tbb::parallel_for_each(blocks.begin(), blocks.end(), [&](const auto& eachBlock) {
    // for(size_t i = 0; i < blocks.size(); i++) {

        std::string consensus;

        int32_t primaryBlockId = ((int32_t)eachBlock.primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)eachBlock.secondaryBlockId);

        for(size_t j = 0; j < eachBlock.consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < panmanUtils::consensusSymbolsPerPackedWord(alphabet); k++) {
                const int nucCode = static_cast<int>(panmanUtils::packedConsensusSymbolAt(eachBlock.consensusSeq[j], k, alphabet));

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = panmanUtils::getSymbolFromCode(nucCode, alphabet);
                consensus.push_back(nucleotide);
            }
            if(endFlag) {
                break;
            }
        }
        std::lock_guard<std::mutex> guard(printMutex);
        fout << ">" << primaryBlockId << "\n" << consensus << "\n";
    });

}

void panmanUtils::Tree::printPseduoRoot(std::ostream& fout){
    std::mutex printMutex;

    fout << ">pseudoRoot" << "\n";
    for(size_t i = 0; i < blocks.size(); i++) {
        auto eachBlock = blocks[i];
        std::string consensus;

        int32_t primaryBlockId = ((int32_t)eachBlock.primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)eachBlock.secondaryBlockId);

        for(size_t j = 0; j < eachBlock.consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < panmanUtils::consensusSymbolsPerPackedWord(alphabet); k++) {
                const int nucCode = static_cast<int>(panmanUtils::packedConsensusSymbolAt(eachBlock.consensusSeq[j], k, alphabet));

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = panmanUtils::getSymbolFromCode(nucCode, alphabet);
                consensus.push_back(nucleotide);
            }
            if(endFlag) {
                break;
            }
        }
        std::lock_guard<std::mutex> guard(printMutex);
        fout << consensus << "\n";
    }
    fout << "\n";

}