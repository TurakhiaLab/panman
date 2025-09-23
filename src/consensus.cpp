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
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (((eachBlock.consensusSeq[j]) >> (4*(7 - k))) & 15);

                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
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