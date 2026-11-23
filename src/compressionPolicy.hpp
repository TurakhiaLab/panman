#pragma once
#include <cstddef>
#include "alphabetCodec.hpp"
#include <cstdint>
#include <set>
#include <vector>

namespace panmanUtils {
constexpr uint64_t largeConsensusThreshold=10000000000ULL;
inline int automaticCompressionLevel(uint64_t bases) {
    return bases>=largeConsensusThreshold?3:9;
}
// Match the writer's per-tree consensus dictionary. Compare existing vectors
// through pointers so counting does not copy potentially large consensuses.
template<class Blocks> uint64_t storedConsensusBases(const Blocks& blocks,Alphabet alphabet) {
    auto less=[](const std::vector<uint32_t>* a,const std::vector<uint32_t>* b){return *a<*b;};
    std::set<const std::vector<uint32_t>*,decltype(less)> seen(less);
    uint64_t total=0;
    for(const auto& block:blocks)if(seen.insert(&block.consensusSeq).second) {
        uint64_t length=block.blockLength>0?uint64_t(block.blockLength):0;
        if(!length) { // Legacy blocks can lack an explicit length.
            bool done=false;
            for(uint32_t word:block.consensusSeq) {
                for(size_t k=0;k<consensusSymbolsPerPackedWord(alphabet);++k) {
                    if(!packedConsensusSymbolAt(word,k,alphabet)){done=true;break;}
                    ++length;
                }
                if(done)break;
            }
        }
        total+=length;
    }
    return total;
}
}
