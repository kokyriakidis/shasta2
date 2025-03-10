#include "extractKmer128.hpp"
#include "LongBaseSequence.hpp"
#include "ShortBaseSequence.hpp"

// ************** TAKE OUT AFTER DEBUGGING.
#include <iostream.hpp>

// Extract n bases at the fiven position from a LongBaseSequenceView
// and store them in the first n bases of a 128-base kmer.
// Uses bit operations for speed, so it can be used in
// performance critical code.

// In both the LongBaseSequenceView and the ShortBaseSequence, base 0
// is the most significant bit.

// We have to process up to 3 pairs of 2 words in the LongBaseSequenceView.

void shasta::extractKmer128(
    const LongBaseSequenceView& v,
    uint64_t position,
    uint64_t n,
    ShortBaseSequence<__uint128_t>& kmer)
{
    // Sanity checks.
    SHASTA_ASSERT(n <= 128);
    SHASTA_ASSERT(position + n <= v.baseCount);

    // Start by clearing the kmer.
    kmer.data[0] = __uint128_t(0);
    kmer.data[1] = __uint128_t(0);

    // Access the first two LongBaseSequenceView words containing the k-mer we want.
    const uint64_t i0 = (position >> 6) << 1;
    const uint64_t i1 = i0 + 1;
    array<uint64_t, 2> ww01 = {v.begin[i0], v.begin[i1]};

    // The starting position of the k-mer in the first two words.
    const uint64_t position01 = position & 63;

    // The number of k-mer bases in the first two words.
    const uint64_t length01 = min(n, 64 - position01);


    // Store these length01 bases at the beginning of the Kmer.
    uint64_t positionInKmer = 0;
    extractBits128(&(ww01[0]), position01, length01, &(kmer.data[0]), positionInKmer);

    // Update the number of bases we still need and the current position in the Kmer.
    n -= length01;
    positionInKmer += length01;

    // If we don't need more bases, we are done.
    if(n == 0) {
        return;
    }

    // We need to get more bases from the next two words.
    const uint64_t length23 = min(n, 64UL);
    const uint64_t i2 = i1 + 1;
    const uint64_t i3 = i2 + 1;
    array<uint64_t, 2> ww23 = {v.begin[i2], v.begin[i3]};
    extractBits128(&(ww23[0]), 0, length23, &(kmer.data[0]), positionInKmer);

    // Update the number of bases we still need and the current position in the Kmer.
    n -= length23;
    positionInKmer += length23;

    // If we don't need more bases, we are done.
    if(n == 0) {
        return;
    }

    // We need to get the remaining bases from the next two words.
    // The number of k-mer bases in the second two words
    const uint64_t length45 = min(n, 64UL);
    const uint64_t i4 = i3 + 1;
    const uint64_t i5 = i4 + 1;
    array<uint64_t, 2> ww45 = {v.begin[i4], v.begin[i5]};
    extractBits128(&(ww45[0]), 0, length45, &(kmer.data[0]), positionInKmer);

    // There can't be any more bases to add.
    n -= length45;
    SHASTA_ASSERT(n == 0);


}



// Extract n bits from x, starting at position xPosition,
// and store them in y, starting at position yPosition,
// leaving the remaining bits of y unchanged.
// Here, xPosition and yPosition are counted with
// 0 at the most significant bit and moving towards the least
// significant bit.
// THIS IS DONE SEPARATELY ON x[0], y[0] AND x[1], y[1].
void shasta::extractBits128(
    const uint64_t* x,
    uint64_t xPosition, // 0 = MSB
    uint64_t n,
    __uint128_t* y,
    uint64_t yPosition  // 0 = MSB
    )
{
    using std::cout;
    using std::endl;

    uint64_t x0 = x[0];
    uint64_t x1 = x[1];
    __uint128_t& y0 = y[0];
    __uint128_t& y1 = y[1];


    SHASTA_ASSERT(xPosition + n <= 64);

    // Shift x right so the n bits are the least significant.
    const uint64_t xShift = 64 - xPosition - n;
    x0 >>= xShift;
    x1 >>= xShift;

    // Copy to an __uint128_t.
    __uint128_t z0 = __uint128_t(x0);
    __uint128_t z1 = __uint128_t(x1);

    // Shift left so the n bits are in the right place.
    const uint64_t zShift = 128 - yPosition - n;
    z0 <<= zShift;
    z1 <<= zShift;

    // Copy these n bits to y without changing the remaining bits.
    const __uint128_t one = __uint128_t(1UL);
    const __uint128_t zMask = ((one << n) - one) << zShift;
    const __uint128_t yMask = ~zMask;
    y0 = (y0 & yMask) | (z0 & zMask);
    y1 = (y1 & yMask) | (z1 & zMask);

}



void shasta::testExtractKmer128()
{
    const string longBaseSequenceString =
        "TCTAGAGTGTCTAATATAGCACTTCTCTTGTGATAAGGGTCCCATAAATA"
        "TTTTAATATAAGTGAATGAATAAATAACCCCATAAATGAAAGAAATCTAG"
        "AATACAAAATTAAATAGATTCAAAAATGAGTTGTGAACTCCTTCTAAGAT"
        "AAGATGAATGAATTGAAACTTTAAATATGAATTTGAGAGGTCAAAATACC"
        "TTCCCATAATCATTAAATAATTAAATTTGACTAGTAATAAATATGGAGAA"
        "CCTCTGGAGTCAAATATTTTGTGAAGTTTGAAAAACTGAGAAGGGTATAT"
        "TTCTTATACCTTAAAGTACCCTTTAAGATAAACCTACTATATTAGAAATA"
        "ATTTATTCTTAGGTTACTTCATGAACATTCACCATGGTTATTTCTGTTAC"
        "ACAAAAGTATCAATTAGCATTCAGTTGCAAAGTAAAACTAAACAAATAAA"
        "AATAGTCCATTGCTTTAAAATGTTTACCTTTGAGTGGTACTATAAGAGAG"
        "CAGAAGAAACTATAGGTGAACATTAATTTACAAGCAAAAGAAGATTCACA"
        "CAAATAAGATTACTAGACTTCATTAATAAAAGCACTTTTTAAAATAACAT"
        "TAATGTATTCACAAAACAAACTGGAAAATAATACAACATATATGACAGTC";

    const uint64_t N = longBaseSequenceString.size();
    cout << "N " << N << endl;

    LongBaseSequence sequence(N);
    for(uint64_t i=0; i<N; i++) {
        sequence.set(i, Base::fromCharacter(longBaseSequenceString[i]));
    }


    for(uint64_t startPosition = 0; startPosition<N; startPosition++) {

        for(uint64_t n=1; n<=min(128UL, N - startPosition); n++) {
            cout << startPosition << " " << n << endl;
            SHASTA_ASSERT(n <= 128);

            ShortBaseSequence128 kmer;
            extractKmer128(sequence, startPosition, n, kmer);

            for(uint64_t i=0; i<n; i++) {
                SHASTA_ASSERT(kmer[i] == sequence[startPosition+i]);
            }
            for(uint64_t i=n; i<128; i++) {
                SHASTA_ASSERT(kmer[i] == Base::fromInteger(uint64_t(0)));
            }
        }
    }

}


