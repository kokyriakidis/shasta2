#pragma once

#include "Base.hpp"
#include "SHASTA2_ASSERT.hpp"

#include "cstdint.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class AlignedBase;

    // ExtendedBase is the 8 symbol alphabet itself: A, C, G, T and the four poly
    // symbols. AlignedExtendedBase is the same plus a gap, for use in an
    // alignment. This mirrors Base and AlignedBase, and for the same reason: a
    // sequence cannot contain a gap, so a type that cannot represent one is the
    // honest type for a sequence.
    class ExtendedBase;
    class ExtendedBaseInitializer;
    class AlignedExtendedBase;
    class AlignedExtendedBaseInitializer;
    inline ostream& operator<<(ostream&, ExtendedBase);
    inline ostream& operator<<(ostream&, AlignedExtendedBase);

    // A sequence encoded in the extended alphabet: each symbol paired with the
    // number of bases it stands for. That is 1 for a plain symbol and the run
    // length for a poly symbol.
    //
    // The length is carried alongside the symbol rather than inside it, so that
    // two poly symbols for the same base compare equal however long their runs.
    // If the length were part of the symbol, the 19 read locus this was
    // developed on would encode to 11 distinct symbol strings instead of 1, and
    // the aligner would treat a run length difference as a mismatch. That is the
    // exact defect the alphabet exists to remove.
    using ExtendedSequence = vector< pair<ExtendedBase, uint64_t> >;

    // One row of an alignment of ExtendedSequences, or an aligned consensus.
    // A gap has length 0.
    using AlignedExtendedSequence = vector< pair<AlignedExtendedBase, uint64_t> >;

    inline string toString(const ExtendedSequence&);
    inline string toString(const AlignedExtendedSequence&);
    inline string toString(const vector<ExtendedBase>&);
    inline string toString(const vector<AlignedExtendedBase>&);
    inline vector<ExtendedBase> vectorOfExtendedBasesFromString(const string&);
    inline vector<AlignedExtendedBase> vectorOfAlignedExtendedBasesFromString(const string&);

    // The default homopolymer threshold used by encodeExtended.
    // A homopolymer run longer than this has its length treated as unreliable
    // and collapses to a poly symbol; a run this long or shorter keeps one plain
    // symbol per base, so its length still takes part in the alignment.
    //
    // Four, because that is where the measured reliability boundary is. Across
    // 91 positions of a 19 read locus, runs of 1, 2, 3 and 4 bases never varied
    // between reads, while runs of 11 and 12 always did. Runs of 4 or shorter
    // are 97% of all runs in these reads, so nearly all length information is
    // kept, and only the unreliable tail is discarded.
    //
    // The threshold must sit BELOW the shortest run whose length is unreliable.
    // A run of exactly the threshold length stays plain while a run one base
    // longer collapses to a single symbol, so if noise can move a run across the
    // threshold, the aligner sees a one base difference as a difference of
    // threshold-1 symbols, and opens a large gap instead of taking one mismatch.
    // Theseus cannot avoid this: it compares characters, so 'a' and 'A' are
    // simply different, and its scalar penalty has no way to price them as
    // similar.
    //
    // Measured on reads simulated from a sequence containing runs of exactly 4,
    // 5, 6, 7 and 12, mean edit distance from the truth:
    //
    //     threshold      3     4     5     6     7     8
    //     runs >=5 noisy 1.95  1.95  3.90  6.80  3.85  3.80
    //     runs >=6 noisy 1.95  1.95  3.90  6.80  3.85    -
    //     runs >=4 noisy 3.27  2.80  4.15  7.35  4.08    -
    //
    // Four is best or tied for best in all three, including the case where runs
    // of 4 are themselves unreliable. Six, the value maxAnchorRepeatLength[0]
    // uses, is the worst of those tried: it lands directly on a run length that
    // is present and noisy. Note maxAnchorRepeatLength exists to keep anchors
    // unique, which is a different question from which run lengths can be
    // trusted, so there is no reason for the two to agree.
    const uint64_t defaultHomopolymerThreshold = 4;

    // Encode a sequence in the extended alphabet.
    // A homopolymer run of base X and length L becomes:
    // - L copies of plain X, each standing for 1 base, if L <= threshold.
    // - a single polyX standing for L bases, if L > threshold.
    // The result is a lossless representation of the input sequence.
    void encodeExtended(
        const vector<Base>& sequence,
        uint64_t threshold,
        ExtendedSequence& encoded);

    // The inverse of encodeExtended.
    void decodeExtended(
        const ExtendedSequence& encoded,
        vector<Base>& sequence);

    // How the consensus length of a long homopolymer run is chosen from the
    // lengths observed in the reads that cover it.
    //
    // Mode is the default, on the strength of the one case here where the true
    // run lengths are known. At a locus whose true lengths are 12 and 11, the 19
    // reads report the first run as:
    //
    //     length  6   7  10  11  12  14  15
    //     reads   1   1   2   6   7   1   1
    //
    // The mode is 12 and is correct. The median is 11 and is wrong, because 10
    // of the 19 reads fall below the truth: rare large deletions (-6, -5) plus
    // six single base under-calls. That is the shape of homopolymer error in ONT
    // reads. Usually correct, occasionally a large deletion, so the distribution
    // is peaked at the truth but has more than half its mass at or below it.
    // Any median type estimator lands one base low on such a distribution, while
    // the mode finds the peak.
    //
    // Which estimator wins is genuinely sensitive to that shape, and simulation
    // was not able to settle it. Drawing length noise from the distribution
    // measured at this run gives mode 1.00 and median 1.15 mean edit distance
    // from the truth; pooling it with the second run, whose errors are far more
    // symmetric, reverses that to mode 0.40 and median 0.23. Both beat the
    // majority voting they replace. Since the simulated answer flips with an
    // assumption that cannot be pinned down from the data available, the real
    // measured case decides it, and that case says mode.
    //
    // Worth revisiting with more loci of known truth. Median stays available.
    enum class RunLengthEstimator {

        // The most frequent length, by total weight. Ties go to the shorter run.
        // The default, see above. Note that column-wise majority voting cannot
        // produce this: a column of a left justified run block is occupied by a
        // majority exactly when over half the reads are at least that long,
        // which is the median.
        Mode,

        // The weighted median length. Biased low when the reads under-call, as
        // they do on long homopolymers.
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
    // Note the two votes happen in this order, and the expansion back to bases
    // happens after both. Expanding first and then voting on the resulting
    // columns cannot give the mode: a column of an expanded run block is
    // occupied by a majority exactly when over half the reads are at least that
    // long, which is the median. It would also let a base bordering a run share
    // a column with the run again, which is the defect being fixed.
    //
    // Build the alignment rows with attachRunLengths().
    void extendedConsensus(
        const vector<AlignedExtendedSequence>& alignment,
        const vector<uint64_t>& weights,
        RunLengthEstimator estimator,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The aligned consensus: one symbol per column, with its voted run
        // length. A gap column has a gap with length 0.
        AlignedExtendedSequence& alignedConsensus);


    // Which ends of a sequence are anchored. This determines which end theseus
    // is allowed to trim, and so which part of the encoding an alignment row
    // corresponds to.
    enum class Anchoring {

        // Fixed on both sides. The row covers the whole encoding.
        BothSides,

        // Fixed on the left only. Theseus trims the overhang off the right, so
        // the row is a prefix of the encoding.
        LeftOnly,

        // Fixed on the right only. Theseus trims the overhang off the left, so
        // the row is a suffix of the encoding.
        RightOnly
    };


    // Put back the run lengths that were held out of the alignment.
    //
    // Theseus is given symbols only, so the rows it returns carry no lengths.
    // This pairs each symbol of a row with the number of bases it stands for,
    // taken from the encoding that was handed to the aligner. A gap gets 0.
    //
    // Dropping the gaps from alignedRow gives a contiguous window of encoded,
    // and not necessarily all of it, because theseus trims the overhang of a
    // sequence fixed on one side only. The window is derived from `anchoring`,
    // NOT by searching the encoding for the row's symbols. Searching is wrong on
    // a repeat: if the encoding is A a C A a C with poly run lengths 10 and 20,
    // and the row is the second copy A a C because the first was trimmed, a
    // search matches the first copy and pairs the row with 10 instead of 20.
    // Tandem repeats are exactly what this code is for, so that is not a remote
    // possibility.
    void attachRunLengths(
        const vector<AlignedExtendedBase>& alignedRow,
        const ExtendedSequence& encoded,
        Anchoring anchoring,
        AlignedExtendedSequence& row);


    // Expand an alignment over the extended alphabet into one over plain bases,
    // by widening each column to the largest run length seen in it.
    // The consensus is widened the same way, so all rows and the aligned
    // consensus keep the same length, as callers require.
    void expandExtendedAlignment(
        const vector<AlignedExtendedSequence>& alignment,
        const AlignedExtendedSequence& alignedConsensus,
        vector< vector<AlignedBase> >& expandedAlignment,
        vector<AlignedBase>& expandedAlignedConsensus);


    // Cheap test for the problem pattern, run on the READS.
    //
    // Returns true if any of the sequences contains two homopolymer runs longer
    // than the threshold separated by a single base. That is the pattern an
    // aligner misplaces, because it can explain the length difference as two
    // substitutions on the separating base rather than as two indels.
    //
    // This is deliberately run on the reads and not on the consensus. The defect
    // can vote the separating base away: at the locus this was developed on the
    // two placements were complementary and exactly one G survived, but with the
    // base split over three or more columns no column need reach a plurality and
    // the consensus would lose it entirely. The cases most in need of repair are
    // exactly the ones the consensus would hide. The reads have not been through
    // any vote, so they still show the pattern.
    //
    // Its purpose is to answer, before an alignment exists, whether it is worth
    // computing one. In production all three local assembly paths set
    // computeAlignment = bool(html), so with html off there is no alignment to
    // repair and one has to be asked for.
    // What makes a region worth repairing.
    enum class Msa1Trigger {

        // Only the problem pattern: two homopolymer runs longer than the
        // threshold separated by a single base. This is the case where an
        // aligner demonstrably misplaces the separating base, because it can
        // explain the run length difference as two substitutions on that base
        // instead of as two indels.
        PatternOnly,

        // Any homopolymer run longer than the threshold. This also covers a run
        // with no bordering base, two runs separated by more than one base, and
        // two long runs side by side, none of which the pattern test sees. It
        // fires far more often, so it repairs far more of the assembly.
        AnyLongRun
    };


    // One sequence. This is the primitive the others are built on. It walks the
    // runs keeping only the two previous run lengths, so it allocates nothing
    // and stops at the first hit.
    bool msa1PatternPresent(
        const vector<Base>& sequence,
        uint64_t threshold = defaultHomopolymerThreshold);

    // True if the sequence contains any homopolymer run longer than the
    // threshold. This is the test for Msa1Trigger::AnyLongRun.
    bool msa1LongRunPresent(
        const vector<Base>& sequence,
        uint64_t threshold = defaultHomopolymerThreshold);

    // Either test, chosen by the trigger.
    bool msa1TriggerPresent(
        const vector<Base>& sequence,
        Msa1Trigger trigger,
        uint64_t threshold = defaultHomopolymerThreshold);

    // Any of several sequences. Note these take DISTINCT sequences: a caller
    // that has repeated each sequence by its coverage, as the abpoa paths do,
    // should pass the distinct ones instead and save the repeats.
    bool msa1PatternPresent(
        const vector< vector<Base> >& sequences,
        uint64_t threshold = defaultHomopolymerThreshold);

    // Same, for sequences that carry a coverage.
    bool msa1PatternPresent(
        const vector< pair<vector<Base>, uint64_t> >& sequences,
        uint64_t threshold = defaultHomopolymerThreshold);


    // A range of alignment columns.
    class Msa1Region {
    public:
        uint64_t begin = 0;
        uint64_t end = 0;      // one past the last column
        uint64_t size() const { return end - begin; }
    };


    // Find the regions of an alignment that are worth repairing.
    //
    // A column is discordant if any row differs from the aligned consensus
    // there. That is exactly what LocalAssembly7::writeAlignment paints pink.
    // Discordant columns are clustered, clusters closer together than
    // mergeDistance are merged, each is widened by flank, and a cluster is kept
    // only if the sequences inside it actually show the problem pattern.
    //
    // The returned regions are disjoint and in increasing order.
    void msa1FindBadRegions(
        const vector< vector<AlignedBase> >& alignment,
        const vector<AlignedBase>& alignedConsensus,
        Msa1Trigger trigger,
        uint64_t threshold,
        uint64_t flank,
        uint64_t mergeDistance,
        vector<Msa1Region>& regions);


    // Repair the bad regions of an alignment and its consensus, in place.
    //
    // This is the entry point. It does NOT replace the alignment: abpoa or
    // theseus run exactly as before, and what they produce is passed here. Only
    // the columns inside a bad region are touched, so sequence accuracy
    // everywhere else is unaffected. On input with no bad region, which is the
    // common case, all three arguments come back unchanged.
    //
    // alignment, alignedConsensus and consensus are all modified in place.
    // Returns the number of regions repaired, which is usually 0.
    uint64_t msa1(

        // The alignment computed by abpoa or theseus, one row per input
        // sequence. All rows must have the same length as alignedConsensus.
        vector< vector<AlignedBase> >& alignment,

        // The aligned consensus, one entry per alignment column.
        vector<AlignedBase>& alignedConsensus,

        // The consensus sequence and its coverage. This is alignedConsensus
        // with the gaps removed, and is kept consistent with it.
        vector< pair<Base, uint64_t> >& consensus,

        // The weight of each row, normally its coverage. May be empty, in which
        // case every row weighs 1.
        const vector<uint64_t>& weights,

        // What makes a region worth repairing. See Msa1Trigger.
        Msa1Trigger trigger = Msa1Trigger::AnyLongRun,

        // The homopolymer threshold. See defaultHomopolymerThreshold.
        uint64_t threshold = defaultHomopolymerThreshold,

        // How the consensus length of a long homopolymer run is chosen.
        RunLengthEstimator estimator = RunLengthEstimator::Mode,

        // Columns of context included on each side of a bad region.
        uint64_t flank = 10,

        // Discordant columns closer together than this are treated as one region.
        uint64_t mergeDistance = 20);


    void testMsa1ExtendedBase();
    void testMsa1Consensus();
    void testMsa1Repair();
}



// Classes used only to store static look up tables used by the fromCharacter
// functions to convert characters to symbols.
class shasta2::ExtendedBaseInitializer{
public:
    ExtendedBaseInitializer();
    static array<uint8_t, 256> table;
    static ExtendedBaseInitializer singleton;
};
class shasta2::AlignedExtendedBaseInitializer{
public:
    AlignedExtendedBaseInitializer();
    static array<uint8_t, 256> table;
    static AlignedExtendedBaseInitializer singleton;
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
// invisible to the aligner and are decided separately, by voting.
//
// For example, with a threshold of 4 and CTC A{12} G A{11} GTT:
//     CTCAAAAAAAAAAAAGAAAAAAAAAAAGTT  ->  CTC polyA G polyA GTT
// and every read over this locus encodes to that same symbol string no matter
// what its two A run lengths are.
//
// The length a poly symbol stands for is deliberately NOT a member of this
// class. It travels beside the symbol, in an ExtendedSequence. Keeping it out is
// what makes the length invisible to alignment: an ExtendedBase holds nothing
// but its symbol value, so two poly symbols for the same base always compare
// equal however long their runs. Adding a length member here, or making
// operator== consider one, would silently undo that.
//
// There is no gap here, for the same reason class Base has none: a sequence
// cannot contain a gap. Use AlignedExtendedBase in an alignment.
//
// Represented as a 1-byte integer:
// - 0, 1, 2, 3 are the plain bases A, C, G, T.
// - 4, 5, 6, 7 are the poly symbols polyA, polyC, polyG, polyT.
// The low two bits are the base and bit 2 is the poly flag, so base() and
// isPoly() are cheap and the plain values agree with class Base.
class shasta2::ExtendedBase {
public:

    // The byte value is always one of 0 through 7.
    uint8_t value;

    // The default constructor constructs A.
    ExtendedBase() : value(0) {}

    bool isValid() const
    {
        return value < 8;
    }

    // We use static member functions instead of constructors.
    // This is safer due to the possibility of unwanted
    // conversions between characters and integers,
    // or confusion between the value stored (0 through 7) and
    // the representing character (A, C, G, T, a, c, g, t).

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

    // Construct from an integer. This does not check validity.
    static ExtendedBase fromInteger(uint8_t i)
    {
        ExtendedBase extendedBase;
        extendedBase.value = i;
        return extendedBase;
    }
    static ExtendedBase fromInteger(uint16_t i) { return fromInteger(uint8_t(i)); }
    static ExtendedBase fromInteger(uint32_t i) { return fromInteger(uint8_t(i)); }
    static ExtendedBase fromInteger(uint64_t i) { return fromInteger(uint8_t(i)); }

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

    // Return true if this is a poly symbol, that is, if it stands for a
    // homopolymer run longer than the threshold used to encode it.
    bool isPoly() const
    {
        return value >= 4;
    }

    // Return the base this symbol stands for, ignoring the poly flag.
    Base base() const
    {
        return Base::fromInteger(uint8_t(value & 3));
    }

    // Return the character representing this symbol.
    // Plain symbols use upper case and poly symbols use lower case, so that a
    // sequence of symbols can be written out and read back as a string. This is
    // how the encoding survives the round trip through theseus, which works on
    // strings.
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
        default:
            throw runtime_error("Invalid extended base value " + to_string(value));
        }
    }

    // Return the complement of this symbol, preserving the poly flag.
    ExtendedBase complement() const
    {
        return fromInteger(uint8_t((value & 4) | (3 - (value & 3))));
    }

    bool operator==(ExtendedBase that) const { return value == that.value; }
    bool operator!=(ExtendedBase that) const { return value != that.value; }
    bool operator<(ExtendedBase that) const { return value < that.value; }

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
        default:
            throw runtime_error("Invalid extended base value " + to_string(value));
        }
    }
};



// Class AlignedExtendedBase is ExtendedBase plus a gap, represented as 8.
// This is to ExtendedBase what AlignedBase is to Base.
class shasta2::AlignedExtendedBase {
public:

    // The byte value is always one of 0 through 8.
    uint8_t value;

    // The default constructor constructs A.
    AlignedExtendedBase() : value(0) {}

    static const uint8_t gapValue = 8;

    bool isValid() const
    {
        return value <= gapValue;
    }

    static AlignedExtendedBase fromCharacter(char c)
    {
        AlignedExtendedBase b;
        b.value = AlignedExtendedBaseInitializer::table[uint8_t(c)];
        if(b.value == 255) {
            string message = "Invalid aligned extended base character: " + to_string(c);
            if(std::isprint(c)) {
                message += ' ';
                message += c;
            }
            throw runtime_error(message);
        }
        return b;
    }

    // Construct from an integer. This does not check validity.
    static AlignedExtendedBase fromInteger(uint8_t i)
    {
        AlignedExtendedBase b;
        b.value = i;
        return b;
    }
    static AlignedExtendedBase fromInteger(uint16_t i) { return fromInteger(uint8_t(i)); }
    static AlignedExtendedBase fromInteger(uint32_t i) { return fromInteger(uint8_t(i)); }
    static AlignedExtendedBase fromInteger(uint64_t i) { return fromInteger(uint8_t(i)); }

    // Construct from an ExtendedBase.
    explicit AlignedExtendedBase(ExtendedBase e) : value(e.value) {}

    // Return a gap.
    static AlignedExtendedBase gap()
    {
        return fromInteger(gapValue);
    }

    bool isGap() const { return value == gapValue; }

    bool isPoly() const { return (value >= 4) and (value < gapValue); }

    // Convert to an ExtendedBase. This asserts if the current value is a gap.
    explicit operator ExtendedBase() const
    {
        SHASTA2_ASSERT(not isGap());
        return ExtendedBase::fromInteger(value);
    }

    // Return the base this symbol stands for. Asserts if this is a gap.
    Base base() const
    {
        SHASTA2_ASSERT(not isGap());
        return Base::fromInteger(uint8_t(value & 3));
    }

    char character() const
    {
        if(value == gapValue) {
            return '-';
        }
        return ExtendedBase::fromInteger(value).character();
    }

    AlignedExtendedBase complement() const
    {
        if(isGap()) {
            return gap();
        }
        return fromInteger(uint8_t((value & 4) | (3 - (value & 3))));
    }

    bool operator==(AlignedExtendedBase that) const { return value == that.value; }
    bool operator!=(AlignedExtendedBase that) const { return value != that.value; }
    bool operator<(AlignedExtendedBase that) const { return value < that.value; }

    string htmlColor() const
    {
        if(value == gapValue) {
            return "";
        }
        return ExtendedBase::fromInteger(value).htmlColor();
    }
};



inline std::ostream& shasta2::operator<<(std::ostream& s, shasta2::ExtendedBase e)
{
    s << e.character();
    return s;
}

inline std::ostream& shasta2::operator<<(std::ostream& s, shasta2::AlignedExtendedBase e)
{
    s << e.character();
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

inline std::string shasta2::toString(const vector<AlignedExtendedBase>& v)
{
    string s;
    for(const AlignedExtendedBase e: v) {
        s.push_back(e.character());
    }
    return s;
}

inline std::string shasta2::toString(const ExtendedSequence& v)
{
    string s;
    for(const auto& [e, runLength]: v) {
        s.push_back(e.character());
    }
    return s;
}

inline std::string shasta2::toString(const AlignedExtendedSequence& v)
{
    string s;
    for(const auto& [e, runLength]: v) {
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

inline std::vector<shasta2::AlignedExtendedBase>
    shasta2::vectorOfAlignedExtendedBasesFromString(const string& s)
{
    vector<AlignedExtendedBase> v;
    for(const char c: s) {
        v.push_back(AlignedExtendedBase::fromCharacter(c));
    }
    return v;
}
