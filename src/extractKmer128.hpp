#pragma once

// Extract n bases at the fiven position from a LongBaseSequenceView
// and store them in the first n bases of a 128-base kmer.
// Uses bit operations for speed, so it can be used in
// performance critical code.

#include <cstdint.hpp>

namespace shasta {
    class LongBaseSequenceView;
    template<class Int> class ShortBaseSequence;

    // Extract n bits from x, starting at position xPosition,
    // and store them in y, starting at position yPosition,
    // leaving the remaining bits of y unchanged.
    // Here, xPosition and yPosition are counted with
    // 0 at the most significant bit and moving towards the least
    // significant bit.
    // The above is repeated for x[0]->y[0], x[1]->y[1].
    void extractBits128(
        const uint64_t* x,
        uint64_t xPosition, // 0 = MSB
        uint64_t n,
        __uint128_t* y,
        uint64_t yPosition  // 0 = MSB
        );

    void extractKmer128(
        const LongBaseSequenceView&,
        uint64_t position,
        uint64_t n,
        ShortBaseSequence<__uint128_t>&);

    void testExtractKmer128();
}
