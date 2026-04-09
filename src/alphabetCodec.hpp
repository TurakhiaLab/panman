#pragma once

#include <cstdint>

namespace panmanUtils {

enum class Alphabet : uint8_t {
    DNA = 0,
    PROTEIN = 1
};

using symbol_code_t = uint8_t;
using packed_word_t = uint32_t;
using mutation_payload_t = uint32_t;

struct AlphabetCodec {
    Alphabet alphabet;
    uint8_t bitsPerSymbol;
    uint8_t stateCount;
    char unknownChar;
    char gapChar;

    constexpr uint32_t symbolMask() const {
        return (1u << bitsPerSymbol) - 1u;
    }

    constexpr uint8_t codesPerWord() const {
        return 32 / bitsPerSymbol;
    }
};

inline const AlphabetCodec& getCodec(Alphabet alphabet) {
    static const AlphabetCodec dnaCodec{
        Alphabet::DNA,
        4,   // 4-bit IUPAC
        16,  // 0..15
        'N',
        '-'
    };
    static const AlphabetCodec proteinCodec{
        Alphabet::PROTEIN,
        5,   // 20 aa + unknown + gap
        22,  // 0..21
        'X',
        '-'
    };
    return (alphabet == Alphabet::PROTEIN) ? proteinCodec : dnaCodec;
}

inline Alphabet& activeAlphabetRef() {
    static Alphabet activeAlphabet = Alphabet::DNA;
    return activeAlphabet;
}

inline Alphabet getActiveAlphabet() {
    return activeAlphabetRef();
}

inline void setActiveAlphabet(Alphabet alphabet) {
    activeAlphabetRef() = alphabet;
}

/** Number of 5-bit/4-bit symbols that fit in the 24-bit mutation payload window. */
inline uint8_t mutationPayloadCapacity(Alphabet a) {
    const auto& c = getCodec(a);
    return static_cast<uint8_t>(24 / c.bitsPerSymbol);
}

inline uint8_t mutationPayloadCapacity() {
    return mutationPayloadCapacity(getActiveAlphabet());
}

/** Bits of the 24-bit wire slot that actually hold symbols (rest are padding when 24 % bps != 0). */
inline uint8_t mutationPayloadStorageBits(Alphabet a) {
    const auto& c = getCodec(a);
    return static_cast<uint8_t>((24u / c.bitsPerSymbol) * c.bitsPerSymbol);
}

inline uint8_t mutationPayloadStorageBits() {
    return mutationPayloadStorageBits(getActiveAlphabet());
}

/** Extract symbol `symbolIndex` from an in-memory mutation payload (high slot = index 0). */
inline uint32_t packedMutationSymbolAt(uint32_t payload, size_t symbolIndex, Alphabet a) {
    const auto& c = getCodec(a);
    const uint8_t cap = static_cast<uint8_t>(24 / c.bitsPerSymbol);
    const uint8_t shift = static_cast<uint8_t>(
        c.bitsPerSymbol * (cap - 1 - static_cast<uint8_t>(symbolIndex)));
    return (payload >> shift) & c.symbolMask();
}

inline uint32_t packedMutationSymbolAt(uint32_t payload, size_t symbolIndex) {
    return packedMutationSymbolAt(payload, symbolIndex, getActiveAlphabet());
}

/** First (only) symbol for SNP-style payloads. */
inline uint32_t singleMutationSymbol(uint32_t payload, Alphabet a) {
    return packedMutationSymbolAt(payload, 0, a);
}

/** Extract symbol `idx` from one packed consensus word (32 bits). */
inline uint32_t packedConsensusSymbolAt(uint32_t packedWord, size_t idx, Alphabet a) {
    const auto& c = getCodec(a);
    const uint8_t cpw = c.codesPerWord();
    const uint8_t shift = static_cast<uint8_t>(
        c.bitsPerSymbol * (cpw - 1 - static_cast<uint8_t>(idx)));
    return (packedWord >> shift) & c.symbolMask();
}

inline uint32_t packedConsensusSymbolAt(uint32_t packedWord, size_t idx) {
    return packedConsensusSymbolAt(packedWord, idx, getActiveAlphabet());
}

inline uint8_t consensusSymbolsPerPackedWord(Alphabet a) {
    return getCodec(a).codesPerWord();
}

inline uint8_t consensusSymbolsPerPackedWord() {
    return consensusSymbolsPerPackedWord(getActiveAlphabet());
}

/** Bits to store in the wire format mutInfo upper bits (capnp/proto). */
inline uint32_t compactedMutationPayloadForWire(uint32_t packedMutation, uint8_t len,
                                                Alphabet a) {
    const auto& c = getCodec(a);
    const uint8_t usedBits = len * c.bitsPerSymbol;
    const uint8_t storageBits = mutationPayloadStorageBits(a);
    if(usedBits >= storageBits) {
        return packedMutation;
    }
    return (packedMutation >> (storageBits - usedBits));
}

inline uint32_t compactedMutationPayloadForWire(uint32_t packedMutation, uint8_t len) {
    return compactedMutationPayloadForWire(packedMutation, len, getActiveAlphabet());
}

/** Inverse of compactedMutationPayloadForWire: wire bits -> in-memory layout for packedMutationSymbolAt. */
inline uint32_t expandMutationPayloadFromWire(uint32_t wirePayload, uint8_t len, Alphabet a) {
    const auto& c = getCodec(a);
    const uint8_t usedBits = len * c.bitsPerSymbol;
    const uint8_t storageBits = mutationPayloadStorageBits(a);
    if(usedBits >= storageBits) {
        return wirePayload;
    }
    return wirePayload << (storageBits - usedBits);
}

inline uint32_t expandMutationPayloadFromWire(uint32_t wirePayload, uint8_t len) {
    return expandMutationPayloadFromWire(wirePayload, len, getActiveAlphabet());
}

}  // namespace panmanUtils
