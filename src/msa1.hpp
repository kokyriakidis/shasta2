#pragma once

#include "Base.hpp"

#include "cstdint.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class AlignedBase;

    class ExtendedBase;
    class ExtendedBaseInitializer;
    inline ostream& operator<<(ostream&, ExtendedBase);
    inline string toString(const vector<ExtendedBase>&);
    inline vector<ExtendedBase> vectorOfExtendedBasesFromString(const string&);

    // The default homopolymer threshold used by encodeExtended.
    // A homopolymer run longer than this has its length treated as unreliable.
    // This matches maxAnchorRepeatLength[0], whose default is 6 (see Options.hpp).
    const uint64_t defaultHomopolymerThreshold = 6;

    // Encode a sequence in the extended alphabet.
    // A homopolymer run of base X and length L becomes:
    // - L copies of plain X, if L <= threshold.
    // - a single polyX, if L > threshold.
    // runLengths is parallel to encoded and holds the number of bases each
    // symbol stands for: always 1 for a plain symbol, and L for a polyX.
    // So the two vectors together are a lossless representation of the input
    // sequence, and decodeExtended needs nothing else.
    void encodeExtended(
        const vector<Base>& sequence,
        uint64_t threshold,
        vector<ExtendedBase>& encoded,
        vector<uint64_t>& runLengths);

    // The inverse of encodeExtended. Gaps in encoded are skipped.
    void decodeExtended(
        const vector<ExtendedBase>& encoded,
        const vector<uint64_t>& runLengths,
        vector<Base>& sequence);

    // How the consensus length of a long homopolymer run is chosen from the
    // lengths observed in the reads that cover it.
    //
    // Median is the default. On reads simulated from a known sequence, with 1%
    // substitutions, 0.5% indels and homopolymer length noise, mean edit
    // distance from the truth over 30 trials at 15x was:
    //
    //     threshold   Mode   Median   theseus majority voting
    //         8       2.23    0.80            1.70
    //         6       2.23    0.93            1.70
    //         5       2.23    0.90            1.70
    //         4       2.23    0.90            1.70
    //
    // So the mode is worse than the majority voting it replaces, and the median
    // is about twice as good as it. The mode is unstable when the observed
    // lengths are spread thinly over several values, which is what a homopolymer
    // at ordinary coverage looks like: in one real 19 read case the two
    // candidate lengths were supported by 7 and 6 reads.
    //
    // Note the simulated length noise is symmetric, which favours the median.
    // Real homopolymer error is biased toward under-calling, so the true gap is
    // probably smaller than the numbers above. Mode remains available.
    enum class RunLengthEstimator {

        // The most frequent length, by total weight. Ties go to the shorter run.
        // Note that column-wise majority voting cannot produce this: a column of
        // a left justified run block is occupied by a majority exactly when over
        // half the reads are at least that long, which is the median.
        Mode,

        // The weighted median length. The default, see above.
        Median
    };


    // Compute a consensus from an alignment over the extended alphabet.
    //
    // This is the step that makes long homopolymer run lengths a voted quantity
    // rather than a by-product of where the aligner placed its gaps. Each column
    // is resolved by weighted majority over the symbols, exactly as before. A
    // column that resolves to a poly symbol then has its run length chosen
    // separately, by voting on the TRUE run lengths of the reads that occupy it.
    //
    // runLengths[i] is parallel to alignment[i]: it holds, for each column, the
    // number of bases that row's symbol stands for at that column, or 0 where
    // the row has a gap. Build it with alignedRunLengths().
    void extendedConsensus(
        const vector< vector<ExtendedBase> >& alignment,
        const vector< vector<uint64_t> >& runLengths,
        const vector<uint64_t>& weights,
        RunLengthEstimator estimator,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The aligned consensus, in the extended alphabet, one entry per column.
        vector<ExtendedBase>& alignedExtendedConsensus,

        // The consensus run length at each column. 1 for a plain symbol, the
        // voted length for a poly symbol, 0 for a gap.
        vector<uint64_t>& consensusRunLengths);


    // Spread the run lengths of an unaligned encoding over the columns of an
    // alignment row, so that they line up with the aligned symbols.
    //
    // alignedRow is one row of an alignment of encoded. Dropping its gaps gives a
    // contiguous window of encoded, and not necessarily all of it: theseus trims
    // the overhang of a sequence fixed on one side only, so a left fixed sequence
    // loses symbols off its right end and a right fixed one off its left. The
    // window is located here, and the columns it does not cover get run length 0.
    void alignedRunLengths(
        const vector<ExtendedBase>& alignedRow,
        const vector<ExtendedBase>& encoded,
        const vector<uint64_t>& encodedRunLengths,
        vector<uint64_t>& alignedRunLengths);


    // Expand an alignment over the extended alphabet into one over plain bases,
    // by widening each column to the largest run length seen in it.
    // The consensus is widened the same way, so all rows and the aligned
    // consensus keep the same length, as callers require.
    void expandExtendedAlignment(
        const vector< vector<ExtendedBase> >& alignment,
        const vector< vector<uint64_t> >& runLengths,
        const vector<ExtendedBase>& alignedExtendedConsensus,
        const vector<uint64_t>& consensusRunLengths,
        vector< vector<AlignedBase> >& expandedAlignment,
        vector<AlignedBase>& expandedAlignedConsensus);


    void testMsa1ExtendedBase();
    void testMsa1Consensus();


    void msa1(

        // The input sequences fixed on both sides, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

        // The input sequences fixed on the left only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

        // The input sequences fixed on the right only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The computed alignment.
        // Each element of the vector correspond to one of the input sequences,
        // in the same order.
        // These all have the same length, which equals the length of the aligned consensus.
        vector< vector<AlignedBase> >& alignment,

        // The aligned consensus.
        vector<AlignedBase>& alignedConsensus,

        // Consensus and alignedConsensus are always computed.
        // Alignment is only computed if this set to true.
        bool computeAlignment,

        // The homopolymer threshold used to encode the sequences.
        // See the comment on defaultHomopolymerThreshold. The default of 6
        // matches maxAnchorRepeatLength[0], but note that on real data a run of
        // exactly 6 then stays plain while longer runs collapse, which costs a
        // spurious mismatch. See testMsa1ExtendedBase.
        uint64_t threshold = defaultHomopolymerThreshold,

        // How the consensus length of a long homopolymer run is chosen.
        RunLengthEstimator estimator = RunLengthEstimator::Median
    );
}



// Class used only to store a static look up table
// used by ExtendedBase::fromCharacter to convert
// characters to extended bases.
class shasta2::ExtendedBaseInitializer{
public:
    ExtendedBaseInitializer();
    static array<uint8_t, 256> table;
    static ExtendedBaseInitializer singleton;
};



// Class ExtendedBase represents a symbol of the extended alphabet used by msa1.
//
// Aligning ONT reads in the plain ACGT alphabet misplaces bases that border a
// long homopolymer run. A gap costs gap_open + gap_extend while a mismatch
// costs only gap_open, so the aligner prefers to explain a run length
// difference as two substitutions rather than two indels, and slides the
// bordering base into the run. The result is an alignment column holding two
// different bases: a substitution that is not real.
//
// Plain run-length encoding (see rle.hpp) removes that ambiguity but overshoots.
// It discards the length of every run, including short ones, and short run
// lengths are reliable: in ONT reads about 97% of runs are 4 bases or shorter,
// and those lengths agree across reads. Plain RLE also makes AC and AAC
// indistinguishable.
//
// The extended alphabet keeps short run lengths and discards only long ones.
// A run no longer than a threshold encodes to one plain symbol per base, so its
// length still takes part in the alignment. A run longer than the threshold
// encodes to a single poly symbol, whatever its length, so its length takes no
// part in the alignment at all. Length differences in long runs are therefore
// invisible to the aligner and are decided separately, by voting on runLengths.
//
// For example, with a threshold of 5 and CTC A{12} G A{11} GTT:
//     CTCAAAAAAAAAAAAGAAAAAAAAAAAGTT  ->  CTC polyA G polyA GTT
// and every read over this locus encodes to that same symbol string no matter
// what its two A run lengths are.
//
// A poly symbol also has a length: the number of bases in the run it stands
// for. That length is deliberately NOT a member of this class. It is returned
// alongside the symbols, in the runLengths argument of encodeExtended.
//
// Keeping it out of the class is what makes the length invisible to alignment.
// An ExtendedBase holds nothing but its symbol value, so two poly symbols for
// the same base always compare equal however long their runs, and an aligner
// working on a vector<ExtendedBase> cannot turn a run length difference into a
// mismatch. The lengths are used only at the end, to decide the consensus
// length of each homopolymer run. Adding a length member here, or making
// operator== consider one, would silently undo that.
//
// Represented as a 1-byte integer:
// - 0, 1, 2, 3 are the plain bases A, C, G, T.
// - 4, 5, 6, 7 are the poly symbols polyA, polyC, polyG, polyT.
// - 8 is a gap.
// The low two bits are the base and bit 2 is the poly flag, so base() and
// isPoly() are cheap and the plain values agree with class Base.
class shasta2::ExtendedBase {
public:

    // The byte value is always one of 0 through 8.
    uint8_t value;

    // The default constructor constructs A.
    ExtendedBase() : value(0) {}

    static const uint8_t gapValue = 8;

    bool isValid() const
    {
        return value <= gapValue;
    }

    // We use static member functions instead of constructors.
    // This is safer due to the possibility of unwanted
    // conversions between characters and integers,
    // or confusion between the value stored (0 through 8) and
    // the representing character (A, C, G, T, a, c, g, t, -).

    // Construct from a character.
    // Throw an exception if the character does not represent a valid symbol.
    static ExtendedBase fromCharacter(char c)
    {
        ExtendedBase extendedBase;
        extendedBase.value = ExtendedBaseInitializer::table[uint8_t(c)];
        if(extendedBase.value == 255) {
            string message = "Invalid extended base character: " + to_string(c);
            if(std::isprint(c)) {
                message += ' ';
                message += c;
            }
            throw runtime_error(message);
        }
        return extendedBase;
    }

    // Construct from an integer.
    // This does not check validity.
    static ExtendedBase fromInteger(uint8_t i)
    {
        ExtendedBase extendedBase;
        extendedBase.value = i;
        return extendedBase;
    }
    static ExtendedBase fromInteger(uint16_t i)
    {
        return fromInteger(uint8_t(i));
    }
    static ExtendedBase fromInteger(uint32_t i)
    {
        return fromInteger(uint8_t(i));
    }
    static ExtendedBase fromInteger(uint64_t i)
    {
        return fromInteger(uint8_t(i));
    }

    // Construct the plain symbol for a Base.
    static ExtendedBase fromBase(Base base)
    {
        return fromInteger(base.value);
    }

    // Construct the poly symbol for a Base.
    static ExtendedBase polyFromBase(Base base)
    {
        return fromInteger(uint8_t(base.value + 4));
    }

    // Return a gap.
    static ExtendedBase gap()
    {
        return fromInteger(gapValue);
    }

    // Return true if this is a poly symbol, that is, if it stands for the tail
    // of a homopolymer run longer than the threshold used to encode it.
    bool isPoly() const
    {
        return (value >= 4) and (value < gapValue);
    }

    bool isGap() const
    {
        return value == gapValue;
    }

    // Return the base this symbol stands for, ignoring the poly flag.
    // This asserts if the current value is a gap.
    Base base() const
    {
        SHASTA2_ASSERT(not isGap());
        return Base::fromInteger(uint8_t(value & 3));
    }

    // Return the character representing this symbol.
    // Plain symbols use upper case and poly symbols use lower case, so that a
    // vector<ExtendedBase> can be written out and read back as a string.
    char character() const
    {
        switch(value) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return 'a';
        case 5: return 'c';
        case 6: return 'g';
        case 7: return 't';
        case gapValue: return '-';
        default:
            throw runtime_error("Invalid extended base value " + to_string(value));
        }
    }

    // Return the complement of this symbol, preserving the poly flag,
    // or the gap if it is already a gap.
    ExtendedBase complement() const
    {
        if(isGap()) {
            return gap();
        }
        return fromInteger(uint8_t((value & 4) | (3 - (value & 3))));
    }

    bool operator==(ExtendedBase that) const
    {
        return value == that.value;
    }
    bool operator!=(ExtendedBase that) const
    {
        return value != that.value;
    }
    bool operator<(ExtendedBase that) const
    {
        return value < that.value;
    }

    // The html color used to represent this symbol.
    // A poly symbol uses a lighter shade of the color of its base, so that the
    // two are recognizably related in an alignment display.
    string htmlColor() const
    {
        switch(value) {
        case 0: return "#ff6666";
        case 1: return "#6666ff";
        case 2: return "#ffff66";
        case 3: return "#66ff66";
        case 4: return "#ffb3b3";
        case 5: return "#b3b3ff";
        case 6: return "#ffffb3";
        case 7: return "#b3ffb3";
        case gapValue: return "";
        default:
            throw runtime_error("Invalid extended base value " + to_string(value));
        }
    }
};



inline std::ostream& shasta2::operator<<(
    std::ostream& s,
    shasta2::ExtendedBase extendedBase)
{
    s << extendedBase.character();
    return s;
}



inline std::string shasta2::toString(const vector<ExtendedBase>& v)
{
    string s;
    for(const ExtendedBase e: v) {
        s.push_back(e.character());
    }
    return s;
}



inline std::vector<shasta2::ExtendedBase> shasta2::vectorOfExtendedBasesFromString(const string& s)
{
    vector<ExtendedBase> v;
    for(const char c: s) {
        v.push_back(ExtendedBase::fromCharacter(c));
    }
    return v;
}



