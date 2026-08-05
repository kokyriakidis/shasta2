// Shasta2.
#include "msa1.hpp"
#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "theseusWrapper.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include "iostream.hpp"



// The look up tables used by the fromCharacter functions.
// Note these differ from BaseInitializer and AlignedBaseInitializer, which map
// lower case characters to the same values as upper case. Here lower case means
// a poly symbol, so 'a' and 'A' are deliberately different.
ExtendedBaseInitializer ExtendedBaseInitializer::singleton;
std::array<uint8_t, 256> ExtendedBaseInitializer::table;
ExtendedBaseInitializer::ExtendedBaseInitializer()
{
    fill(table.begin(), table.end(), 255);
    table['A'] = 0;
    table['C'] = 1;
    table['G'] = 2;
    table['T'] = 3;
    table['a'] = 4;
    table['c'] = 5;
    table['g'] = 6;
    table['t'] = 7;
}

AlignedExtendedBaseInitializer AlignedExtendedBaseInitializer::singleton;
std::array<uint8_t, 256> AlignedExtendedBaseInitializer::table;
AlignedExtendedBaseInitializer::AlignedExtendedBaseInitializer()
{
    fill(table.begin(), table.end(), 255);
    table['A'] = 0;
    table['C'] = 1;
    table['G'] = 2;
    table['T'] = 3;
    table['a'] = 4;
    table['c'] = 5;
    table['g'] = 6;
    table['t'] = 7;
    table['-'] = AlignedExtendedBase::gapValue;
}



// See msa1.hpp for comments.
void shasta2::encodeExtended(
    const vector<Base>& sequence,
    uint64_t threshold,
    ExtendedSequence& encoded)
{
    // A threshold of 0 would make every run, including a run of length 1,
    // encode to a bare poly symbol. That loses the length of every run and is
    // exactly the plain RLE behavior this alphabet exists to avoid.
    SHASTA2_ASSERT(threshold > 0);

    encoded.clear();

    // The encoding never has more symbols than the sequence has bases.
    encoded.reserve(sequence.size());

    const uint64_t n = sequence.size();
    for(uint64_t runBegin=0; runBegin<n; /* incremented below */) {
        const Base base = sequence[runBegin];

        // Find the end of the homopolymer run that begins here.
        uint64_t runEnd = runBegin + 1;
        while((runEnd < n) and (sequence[runEnd] == base)) {
            ++runEnd;
        }
        const uint64_t runLength = runEnd - runBegin;

        if(runLength > threshold) {

            // A long run becomes a single poly symbol standing for the whole
            // run. Its length is carried alongside and plays no part in any
            // alignment done on the symbols.
            encoded.push_back(make_pair(ExtendedBase::polyFromBase(base), runLength));

        } else {

            // A short run keeps one plain symbol per base, so its length is
            // still visible to an aligner. Short run lengths are reliable and
            // are real signal, unlike long ones.
            for(uint64_t i=0; i<runLength; i++) {
                encoded.push_back(make_pair(ExtendedBase::fromBase(base), 1UL));
            }
        }

        runBegin = runEnd;
    }
}



// See msa1.hpp for comments.
void shasta2::decodeExtended(
    const ExtendedSequence& encoded,
    vector<Base>& sequence)
{
    sequence.clear();

    uint64_t n = 0;
    for(const auto& [extendedBase, runLength]: encoded) {
        n += runLength;
    }
    sequence.reserve(n);

    for(const auto& [extendedBase, runLength]: encoded) {
        const Base base = extendedBase.base();
        for(uint64_t j=0; j<runLength; j++) {
            sequence.push_back(base);
        }
    }
    SHASTA2_ASSERT(sequence.size() == n);
}



// See msa1.hpp for comments.
void shasta2::attachRunLengths(
    const vector<AlignedExtendedBase>& alignedRow,
    const ExtendedSequence& encoded,
    Anchoring anchoring,
    AlignedExtendedSequence& row)
{
    // Count the non-gap symbols of the row.
    uint64_t rowSymbolCount = 0;
    for(const AlignedExtendedBase e: alignedRow) {
        if(not e.isGap()) {
            ++rowSymbolCount;
        }
    }
    SHASTA2_ASSERT(rowSymbolCount <= encoded.size());

    // Locate the part of the encoding that this row covers.
    // Which end theseus may have trimmed is known from the anchoring and is not
    // inferred: see the comment in msa1.hpp for why searching the encoding for
    // the row's symbols gives the wrong answer on a repeat.
    uint64_t offset = 0;
    switch(anchoring) {
    case Anchoring::BothSides:
        SHASTA2_ASSERT(rowSymbolCount == encoded.size());
        offset = 0;
        break;
    case Anchoring::LeftOnly:
        // Anchored on the left, trimmed on the right: a prefix.
        offset = 0;
        break;
    case Anchoring::RightOnly:
        // Anchored on the right, trimmed on the left: a suffix.
        offset = encoded.size() - rowSymbolCount;
        break;
    }

    row.clear();
    row.reserve(alignedRow.size());
    uint64_t i = 0;
    for(const AlignedExtendedBase e: alignedRow) {
        if(e.isGap()) {
            row.push_back(make_pair(e, 0UL));
        } else {
            const auto& [encodedBase, runLength] = encoded[offset + i];
            SHASTA2_ASSERT(e == AlignedExtendedBase(encodedBase));
            row.push_back(make_pair(e, runLength));
            ++i;
        }
    }
    SHASTA2_ASSERT(i == rowSymbolCount);
}



// See msa1.hpp for comments.
void shasta2::extendedConsensus(
    const vector<AlignedExtendedSequence>& alignment,
    const vector<uint64_t>& weights,
    RunLengthEstimator estimator,
    vector< pair<Base, uint64_t> >& consensus,
    AlignedExtendedSequence& alignedConsensus)
{
    consensus.clear();
    alignedConsensus.clear();

    const uint64_t n = alignment.size();
    if(n == 0) {
        return;
    }
    SHASTA2_ASSERT(weights.size() == n);

    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=0; i<n; i++) {
        SHASTA2_ASSERT(alignment[i].size() == alignmentLength);
    }

    alignedConsensus.resize(alignmentLength);

    // Scratch, reused by the run length vote at every poly column.
    vector<uint64_t> lengthWeight;

    // Loop over alignment columns.
    for(uint64_t j=0; j<alignmentLength; j++) {

        // Total weight for each symbol at this column, including the gap.
        array<uint64_t, AlignedExtendedBase::gapValue + 1> symbolWeight;
        fill(symbolWeight.begin(), symbolWeight.end(), 0);
        for(uint64_t i=0; i<n; i++) {
            symbolWeight[alignment[i][j].first.value] += weights[i];
        }

        // The consensus symbol is the one with the most weight.
        const auto it = std::ranges::max_element(symbolWeight);
        const AlignedExtendedBase consensusSymbol =
            AlignedExtendedBase::fromInteger(uint64_t(it - symbolWeight.begin()));

        if(consensusSymbol.isGap()) {
            alignedConsensus[j] = make_pair(consensusSymbol, 0UL);
            continue;
        }

        // The weight supporting this symbol.
        uint64_t coverage = *it;

        // The consensus run length. A plain symbol always stands for one base.
        // A poly symbol gets its length by voting on the TRUE run lengths of the
        // rows that carry the same poly symbol here. This is the whole point of
        // the extended alphabet: the length is decided here, not by counting how
        // many columns the aligner happened to fill.
        uint64_t consensusRunLength = 1;
        if(consensusSymbol.isPoly()) {

            // Gather the observed lengths with their weights.
            // Lengths are bounded by the sequence length, so a vector indexed by
            // length is both simpler and faster than a map. It is reused across
            // columns rather than reallocated at each one, because a long
            // sequence has many poly columns and each allocation would be
            // proportional to a run length.
            uint64_t maxObserved = 0;
            uint64_t totalWeight = 0;
            for(uint64_t i=0; i<n; i++) {
                if(alignment[i][j].first == consensusSymbol) {
                    maxObserved = max(maxObserved, alignment[i][j].second);
                    totalWeight += weights[i];
                }
            }
            SHASTA2_ASSERT(totalWeight > 0);

            lengthWeight.assign(maxObserved + 1, 0);
            for(uint64_t i=0; i<n; i++) {
                if(alignment[i][j].first == consensusSymbol) {
                    lengthWeight[alignment[i][j].second] += weights[i];
                }
            }

            if(estimator == RunLengthEstimator::Mode) {

                // The most frequent length by weight. Ties go to the shorter run,
                // which the strict > gives us because we scan in increasing order.
                uint64_t bestWeight = 0;
                for(uint64_t length=1; length<=maxObserved; length++) {
                    if(lengthWeight[length] > bestWeight) {
                        bestWeight = lengthWeight[length];
                        consensusRunLength = length;
                    }
                }
                coverage = bestWeight;

            } else {

                // The weighted median.
                uint64_t cumulative = 0;
                for(uint64_t length=1; length<=maxObserved; length++) {
                    cumulative += lengthWeight[length];
                    if(2 * cumulative >= totalWeight) {
                        consensusRunLength = length;
                        coverage = lengthWeight[length];
                        break;
                    }
                }
            }
            SHASTA2_ASSERT(consensusRunLength > 0);
        }

        alignedConsensus[j] = make_pair(consensusSymbol, consensusRunLength);

        // Store the run in the consensus.
        const Base base = consensusSymbol.base();
        for(uint64_t k=0; k<consensusRunLength; k++) {
            consensus.push_back(make_pair(base, coverage));
        }
    }
}



// See msa1.hpp for comments.
void shasta2::expandExtendedAlignment(
    const vector<AlignedExtendedSequence>& alignment,
    const AlignedExtendedSequence& alignedConsensus,
    vector< vector<AlignedBase> >& expandedAlignment,
    vector<AlignedBase>& expandedAlignedConsensus)
{
    expandedAlignment.clear();
    expandedAlignedConsensus.clear();

    const uint64_t n = alignment.size();
    if(n == 0) {
        return;
    }
    const uint64_t alignmentLength = alignment.front().size();
    SHASTA2_ASSERT(alignedConsensus.size() == alignmentLength);

    // The width of each column after expansion: the longest run seen there,
    // including the consensus, and at least 1.
    vector<uint64_t> width(alignmentLength, 1);
    uint64_t expandedLength = 0;
    for(uint64_t j=0; j<alignmentLength; j++) {
        for(uint64_t i=0; i<n; i++) {
            width[j] = max(width[j], alignment[i][j].second);
        }
        width[j] = max(width[j], alignedConsensus[j].second);
        expandedLength += width[j];
    }

    // Expand one row, padding on the right with gaps.
    // Every expanded row has the same length, known here, so it is reserved
    // once instead of being grown a base at a time.
    const auto expand = [&width, alignmentLength, expandedLength](
        const AlignedExtendedSequence& row,
        vector<AlignedBase>& expanded)
    {
        expanded.clear();
        expanded.reserve(expandedLength);
        for(uint64_t j=0; j<alignmentLength; j++) {
            const auto& [e, runLength] = row[j];
            for(uint64_t k=0; k<runLength; k++) {
                expanded.push_back(AlignedBase(e.base()));
            }
            for(uint64_t k=runLength; k<width[j]; k++) {
                expanded.push_back(AlignedBase::gap());
            }
        }
        SHASTA2_ASSERT(expanded.size() == expandedLength);
    };

    expandedAlignment.resize(n);
    for(uint64_t i=0; i<n; i++) {
        expand(alignment[i], expandedAlignment[i]);
    }
    expand(alignedConsensus, expandedAlignedConsensus);

    // All rows and the aligned consensus must have the same length.
    for(const vector<AlignedBase>& row: expandedAlignment) {
        SHASTA2_ASSERT(row.size() == expandedAlignedConsensus.size());
    }
}


void shasta2::msa1(

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
    uint64_t threshold,

    // How the consensus length of a long homopolymer run is chosen.
    RunLengthEstimator estimator
)
{
    SHASTA2_ASSERT(not fixedSequences.empty());

    // Encode all the input sequences in the extended alphabet.
    // A homopolymer run longer than the threshold becomes a single poly symbol
    // carrying its length, and only the symbols are shown to the aligner. This
    // is what stops a run length difference from being turned into a mismatch,
    // which is what misplaces the base bordering a long run.
    vector< pair<ExtendedSequence, uint64_t> > encodedFixed;
    vector< pair<ExtendedSequence, uint64_t> > encodedLeftFixed;
    vector< pair<ExtendedSequence, uint64_t> > encodedRightFixed;

    // How each sequence is anchored and what it weighs, in the same order the
    // rows will come back from the aligner: fixed, then left fixed, then right
    // fixed. The anchoring is needed to put the run lengths back, because it
    // says which end theseus may have trimmed.
    //
    // The encodings themselves are not copied here. They stay in the three
    // group vectors above and are reached through encodingOfRow below, which
    // maps a row index onto them. At the largest MSA this code sees that saves
    // a few megabytes of copying for nothing.
    vector<Anchoring> anchorings;
    vector<uint64_t> weights;

    const auto encodeGroup = [&](
        const vector< pair<vector<Base>, uint64_t> >& in,
        Anchoring anchoring,
        vector< pair<ExtendedSequence, uint64_t> >& out)
    {
        out.clear();
        out.reserve(in.size());
        for(const auto& [sequence, weight]: in) {
            ExtendedSequence encoded;
            encodeExtended(sequence, threshold, encoded);
            out.push_back(make_pair(std::move(encoded), weight));
            anchorings.push_back(anchoring);
            weights.push_back(weight);
        }
    };
    encodeGroup(fixedSequences, Anchoring::BothSides, encodedFixed);
    encodeGroup(leftFixedSequences, Anchoring::LeftOnly, encodedLeftFixed);
    encodeGroup(rightFixedSequences, Anchoring::RightOnly, encodedRightFixed);

    // Map a row index onto the encoding it came from.
    const auto encodingOfRow = [&](uint64_t i) -> const ExtendedSequence&
    {
        if(i < encodedFixed.size()) {
            return encodedFixed[i].first;
        }
        i -= encodedFixed.size();
        if(i < encodedLeftFixed.size()) {
            return encodedLeftFixed[i].first;
        }
        i -= encodedLeftFixed.size();
        return encodedRightFixed[i].first;
    };

    // Align the encoded sequences. Theseus sees symbols only.
    vector< vector<AlignedExtendedBase> > alignedSymbols;
    theseusExtended(encodedFixed, encodedLeftFixed, encodedRightFixed,
        alignedSymbols);
    SHASTA2_ASSERT(alignedSymbols.size() == anchorings.size());

    // Put the run lengths back, using the encoding that was handed to the
    // aligner rather than one recovered from the row: for a sequence fixed on
    // one side only, theseus trims the overhang, so the row covers only part of
    // the encoding, and the anchoring says which part.
    vector<AlignedExtendedSequence> extendedAlignment(alignedSymbols.size());
    for(uint64_t i=0; i<alignedSymbols.size(); i++) {
        attachRunLengths(alignedSymbols[i], encodingOfRow(i), anchorings[i],
            extendedAlignment[i]);
    }

    // Vote, column by column, with long run lengths voted separately.
    AlignedExtendedSequence alignedExtendedConsensus;
    extendedConsensus(extendedAlignment, weights, estimator,
        consensus, alignedExtendedConsensus);

    // Expand back to plain bases. This is done even when the alignment is not
    // requested, because alignedConsensus is always computed and has to have the
    // same length as the alignment rows. When it is not requested the expanded
    // rows are discarded, but they still have to be built to know that length.
    vector< vector<AlignedBase> > expandedAlignment;
    expandExtendedAlignment(extendedAlignment, alignedExtendedConsensus,
        expandedAlignment, alignedConsensus);

    if(computeAlignment) {
        alignment = std::move(expandedAlignment);
    } else {
        alignment.clear();
    }
}



// The 19 bubble arm sequences this code was developed against.
// They are reads over a single locus with structure CTC A{n1} G A{n2} GTT,
// inside a 2x tandem duplication. Only n1 and n2 vary between them, and the
// TRUE values are known: n1 = 12 and n2 = 11.
namespace shasta2 {
    static const vector<string> msa1TestSequences = {
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
    };

    // The alignment abpoa actually produces for those sequences. Four of its
    // columns hold a mixture of A and G, because the G between the two long A
    // runs was placed in two different columns.
    static const vector<string> msa1BadAlignmentRows = {
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---------AAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---------AAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC-------AAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC-------AAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
    };

    // Count alignment columns that hold more than one distinct base.
    static uint64_t msa1ImpureColumnCount(const vector<string>& rows)
    {
        const uint64_t alignmentLength = rows.front().size();
        uint64_t impureColumnCount = 0;
        for(uint64_t j=0; j<alignmentLength; j++) {
            array<bool, 4> seen = {false, false, false, false};
            for(const string& row: rows) {
                SHASTA2_ASSERT(row.size() == alignmentLength);
                const char c = row[j];
                if(c != '-') {
                    seen[Base::fromCharacter(c).value] = true;
                }
            }
            uint64_t distinctCount = 0;
            for(const bool b: seen) {
                if(b) {
                    ++distinctCount;
                }
            }
            if(distinctCount > 1) {
                ++impureColumnCount;
            }
        }
        return impureColumnCount;
    }

    static string msa1ToString(const vector< pair<Base, uint64_t> >& consensus)
    {
        string s;
        for(const auto& [base, coverage]: consensus) {
            s.push_back(base.character());
        }
        return s;
    }
}



// Test the extended alphabet and its encoding.
void shasta2::testMsa1ExtendedBase()
{
    // Check the symbol values are what the low-bits arithmetic in base() and
    // isPoly() assumes.
    for(uint64_t i=0; i<4; i++) {
        const Base base = Base::fromInteger(i);
        const ExtendedBase plain = ExtendedBase::fromBase(base);
        const ExtendedBase poly = ExtendedBase::polyFromBase(base);
        SHASTA2_ASSERT(plain.value == i);
        SHASTA2_ASSERT(poly.value == i + 4);
        SHASTA2_ASSERT(not plain.isPoly());
        SHASTA2_ASSERT(poly.isPoly());
        SHASTA2_ASSERT(plain.base() == base);
        SHASTA2_ASSERT(poly.base() == base);
        SHASTA2_ASSERT(plain.isValid());
        SHASTA2_ASSERT(poly.isValid());

        // The complement preserves the poly flag.
        SHASTA2_ASSERT(plain.complement().base() == base.complement());
        SHASTA2_ASSERT(poly.complement().base() == base.complement());
        SHASTA2_ASSERT(not plain.complement().isPoly());
        SHASTA2_ASSERT(poly.complement().isPoly());
        SHASTA2_ASSERT(plain.complement().complement() == plain);
        SHASTA2_ASSERT(poly.complement().complement() == poly);

        // AlignedExtendedBase agrees with ExtendedBase on the 8 real symbols.
        SHASTA2_ASSERT(AlignedExtendedBase(plain).value == plain.value);
        SHASTA2_ASSERT(AlignedExtendedBase(poly).value == poly.value);
        SHASTA2_ASSERT(not AlignedExtendedBase(plain).isGap());
        SHASTA2_ASSERT(AlignedExtendedBase(poly).isPoly());
        SHASTA2_ASSERT(ExtendedBase(AlignedExtendedBase(poly)) == poly);
    }

    // The gap lives in AlignedExtendedBase only. ExtendedBase cannot represent
    // one, which is the point: a sequence cannot contain a gap.
    {
        const AlignedExtendedBase gap = AlignedExtendedBase::gap();
        SHASTA2_ASSERT(gap.isGap());
        SHASTA2_ASSERT(not gap.isPoly());
        SHASTA2_ASSERT(gap.isValid());
        SHASTA2_ASSERT(gap.complement().isGap());
        SHASTA2_ASSERT(gap.character() == '-');
        SHASTA2_ASSERT(gap.htmlColor().empty());

        // '-' is not a valid ExtendedBase.
        bool threw = false;
        try {
            ExtendedBase::fromCharacter('-');
        } catch(const std::exception&) {
            threw = true;
        }
        SHASTA2_ASSERT(threw);

        // The alphabet has exactly 8 symbols.
        for(uint64_t value=0; value<8; value++) {
            SHASTA2_ASSERT(ExtendedBase::fromInteger(value).isValid());
        }
        SHASTA2_ASSERT(not ExtendedBase::fromInteger(8UL).isValid());
    }

    // character() and fromCharacter() round trip over all symbols.
    for(uint64_t value=0; value<8; value++) {
        const ExtendedBase e = ExtendedBase::fromInteger(value);
        SHASTA2_ASSERT(ExtendedBase::fromCharacter(e.character()) == e);
    }
    for(uint64_t value=0; value<=AlignedExtendedBase::gapValue; value++) {
        const AlignedExtendedBase e = AlignedExtendedBase::fromInteger(value);
        SHASTA2_ASSERT(AlignedExtendedBase::fromCharacter(e.character()) == e);
    }

    // An invalid character throws.
    {
        bool threw = false;
        try {
            ExtendedBase::fromCharacter('X');
        } catch(const std::exception&) {
            threw = true;
        }
        SHASTA2_ASSERT(threw);
    }

    // toString and the fromString helpers round trip.
    // This matters because it is how the encoding reaches theseus and comes
    // back: theseusExtended passes toString(encoded) in and parses the MSA rows
    // with vectorOfAlignedExtendedBasesFromString. If either mapping folded 'a'
    // into 'A', as AlignedBase::fromCharacter does, every poly symbol would
    // silently become a plain one.
    {
        SHASTA2_ASSERT(toString(vectorOfExtendedBasesFromString("ACGTacgt")) == "ACGTacgt");
        SHASTA2_ASSERT(toString(vectorOfAlignedExtendedBasesFromString("ACGTacgt-")) == "ACGTacgt-");
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a') != ExtendedBase::fromCharacter('A'));
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').isPoly());
        SHASTA2_ASSERT(not ExtendedBase::fromCharacter('A').isPoly());
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').base() ==
            ExtendedBase::fromCharacter('A').base());

        // AlignedBase, by contrast, would lose it. This documents exactly why
        // theseusExtended must not use vectorOfAlignedBasesFromString.
        SHASTA2_ASSERT(AlignedBase::fromCharacter('a') == AlignedBase::fromCharacter('A'));
    }



    // Round trip all 19 sequences, at a range of thresholds, and check that they
    // all encode to the SAME symbol string. That is the property that makes the
    // homopolymer ambiguity structurally impossible: if every read has the same
    // run structure there is nothing to align, and no scoring choice can
    // misplace the G between the two A runs.
    for(uint64_t threshold=1; threshold<=12; threshold++) {
        ExtendedSequence firstEncoded;
        bool uniform = true;
        uint64_t variablePositionCount = 0;

        for(uint64_t i=0; i<msa1TestSequences.size(); i++) {
            const vector<Base> sequence = vectorOfBasesFromString(msa1TestSequences[i]);

            ExtendedSequence encoded;
            encodeExtended(sequence, threshold, encoded);

            // The encoding is lossless.
            vector<Base> decoded;
            decodeExtended(encoded, decoded);
            SHASTA2_ASSERT(decoded == sequence);

            // A plain symbol always stands for exactly one base, and a poly
            // symbol for more than the threshold.
            for(const auto& [e, runLength]: encoded) {
                if(e.isPoly()) {
                    SHASTA2_ASSERT(runLength > threshold);
                } else {
                    SHASTA2_ASSERT(runLength == 1);
                }
            }

            // A poly symbol stands for a whole run on its own, so it is never
            // adjacent to a symbol with the same base.
            for(uint64_t j=1; j<encoded.size(); j++) {
                if(encoded[j].first.isPoly() or encoded[j-1].first.isPoly()) {
                    SHASTA2_ASSERT(encoded[j].first.base() != encoded[j-1].first.base());
                }
            }

            if(i == 0) {
                firstEncoded = encoded;
            } else if(toString(encoded) != toString(firstEncoded)) {
                uniform = false;
            } else {
                SHASTA2_ASSERT(encoded.size() == firstEncoded.size());
                variablePositionCount = 0;
                for(uint64_t j=0; j<encoded.size(); j++) {
                    if(encoded[j].second != firstEncoded[j].second) {
                        ++variablePositionCount;
                    }
                }
            }
        }

        cout << "Threshold " << threshold << ": " <<
            (uniform ? "all 19 sequences encode to the same " : "encodings DIFFER, first has ") <<
            firstEncoded.size() << " symbols." << endl;

        // The shortest A run in this set is 6 bases and the others at the same
        // two positions are 7 or longer, so a threshold below 6 collapses all of
        // them and every sequence encodes identically. At a threshold of 6 the
        // run of 6 is no longer above the threshold, stays plain, and the
        // guarantee is lost. Asserted in both directions.
        if(threshold < 6) {
            SHASTA2_ASSERT(uniform);
        }
        if(threshold == 6) {
            SHASTA2_ASSERT(not uniform);
        }

        // The default must keep the encoding uniform here. This is the property
        // the whole alphabet exists to provide, so a change to the default that
        // loses it should not pass silently.
        if(threshold == defaultHomopolymerThreshold) {
            SHASTA2_ASSERT(uniform);
        }
    }



    // A sequence with no run longer than the threshold encodes to itself, one
    // plain symbol per base, so behavior on sequences without long homopolymers
    // is unchanged by this alphabet.
    {
        const vector<Base> sequence = vectorOfBasesFromString("ACGTACGTAACCGGTT");
        ExtendedSequence encoded;
        encodeExtended(sequence, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == sequence.size());
        for(uint64_t i=0; i<encoded.size(); i++) {
            SHASTA2_ASSERT(not encoded[i].first.isPoly());
            SHASTA2_ASSERT(encoded[i].first.base() == sequence[i]);
            SHASTA2_ASSERT(encoded[i].second == 1);
        }
    }



    // Behavior of a single run as its length crosses the threshold.
    //
    // A run of length threshold encodes to threshold plain symbols, and a run of
    // length threshold+1 to a single poly symbol. So the encoding is
    // discontinuous at the threshold: one extra base changes the symbol string
    // by threshold-1 symbols, and to an aligner a run that straddles the
    // threshold looks like a large indel rather than a small one. That is the
    // price of collapsing a whole run to one symbol, and it is why the threshold
    // has to be kept below the shortest run whose length is unreliable.
    {
        const uint64_t threshold = defaultHomopolymerThreshold;
        for(uint64_t runLength=1; runLength<=threshold+3; runLength++) {
            const vector<Base> sequence(runLength, Base::fromCharacter('A'));
            ExtendedSequence encoded;
            encodeExtended(sequence, threshold, encoded);

            if(runLength > threshold) {
                SHASTA2_ASSERT(encoded.size() == 1);
                SHASTA2_ASSERT(encoded[0].first.isPoly());
                SHASTA2_ASSERT(encoded[0].second == runLength);
            } else {
                SHASTA2_ASSERT(encoded.size() == runLength);
                for(const auto& [e, ignore]: encoded) {
                    SHASTA2_ASSERT(not e.isPoly());
                }
            }

            vector<Base> decoded;
            decodeExtended(encoded, decoded);
            SHASTA2_ASSERT(decoded == sequence);
        }

        // All runs longer than the threshold encode to the same single SYMBOL,
        // regardless of length. This is the property that removes long run
        // lengths from the alignment entirely: the symbol strings are equal, and
        // only the lengths carried alongside differ.
        ExtendedSequence shortRun;
        ExtendedSequence longRun;
        encodeExtended(vector<Base>(threshold + 1, Base::fromCharacter('A')),
            threshold, shortRun);
        encodeExtended(vector<Base>(1000, Base::fromCharacter('A')),
            threshold, longRun);
        SHASTA2_ASSERT(toString(shortRun) == toString(longRun));
        SHASTA2_ASSERT(shortRun.size() == 1);
        SHASTA2_ASSERT(shortRun[0].second != longRun[0].second);

        // The discontinuity at the threshold.
        ExtendedSequence atThreshold;
        encodeExtended(vector<Base>(threshold, Base::fromCharacter('A')),
            threshold, atThreshold);
        SHASTA2_ASSERT(atThreshold.size() == threshold);
    }



    // Degenerate inputs.
    {
        ExtendedSequence encoded;
        vector<Base> decoded;

        // Empty sequence.
        encodeExtended(vector<Base>(), defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.empty());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded.empty());

        // Single base.
        const vector<Base> single = vectorOfBasesFromString("G");
        encodeExtended(single, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(not encoded[0].first.isPoly());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded == single);

        // One run spanning the whole sequence: a single poly symbol.
        const vector<Base> oneRun(100, Base::fromCharacter('T'));
        encodeExtended(oneRun, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(encoded[0].first.isPoly());
        SHASTA2_ASSERT(encoded[0].second == oneRun.size());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded == oneRun);
    }



    // Round trip random sequences over a range of thresholds.
    // A simple deterministic generator, so a failure is reproducible.
    {
        uint64_t random = 38495721;
        const auto next = [&random]() {
            random = random * 6364136223846793005ULL + 1442695040888963407ULL;
            return uint64_t(random >> 33);
        };

        for(uint64_t threshold=1; threshold<=8; threshold++) {
            for(uint64_t trial=0; trial<200; trial++) {
                vector<Base> sequence;
                const uint64_t runCount = next() % 40;
                Base previousBase = Base::fromInteger(uint64_t(next() % 4));
                for(uint64_t i=0; i<runCount; i++) {
                    // Consecutive runs must have different bases, so the run
                    // structure is what we intend it to be.
                    Base base = previousBase;
                    if(i > 0) {
                        base = Base::fromInteger(uint64_t((previousBase.value + 1 + next() % 3) % 4));
                        SHASTA2_ASSERT(base != previousBase);
                    }
                    const uint64_t runLength = 1 + next() % 20;
                    for(uint64_t j=0; j<runLength; j++) {
                        sequence.push_back(base);
                    }
                    previousBase = base;
                }

                ExtendedSequence encoded;
                encodeExtended(sequence, threshold, encoded);

                vector<Base> decoded;
                decodeExtended(encoded, decoded);
                SHASTA2_ASSERT(decoded == sequence);

                // The encoding never grows the sequence.
                SHASTA2_ASSERT(encoded.size() <= max(sequence.size(), 1UL));
            }
        }
    }

    cout << "testMsa1ExtendedBase passed." << endl;
}


// Test the consensus computed from an alignment over the extended alphabet.
//
// Theseus is not exercised here: extendedConsensus, attachRunLengths and
// expandExtendedAlignment are pure functions, so the alignment is supplied
// directly. That also lets the real, defective abpoa alignment be used as input.
void shasta2::testMsa1Consensus()
{
    // Confirm the input really is defective. If this stops being true the test
    // has lost its point.
    {
        const uint64_t impureColumnCount = msa1ImpureColumnCount(msa1BadAlignmentRows);
        cout << "The input alignment has " << impureColumnCount <<
            " columns containing more than one base." << endl;
        SHASTA2_ASSERT(impureColumnCount == 4);
    }

    // Recover the sequences from that alignment and encode them. They all encode
    // to the same symbol string, so they line up without any further work; this
    // stands in for what theseus would return. 11 distinct base sequences
    // collapse to 1 distinct encoding.
    const uint64_t threshold = defaultHomopolymerThreshold;
    vector<AlignedExtendedSequence> alignment;
    vector<uint64_t> weights;
    for(const string& row: msa1BadAlignmentRows) {
        string ungapped;
        for(const char c: row) {
            if(c != '-') {
                ungapped.push_back(c);
            }
        }
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(ungapped), threshold, encoded);

        // No gaps needed: every row has the same symbol string.
        AlignedExtendedSequence alignedRow;
        for(const auto& [e, runLength]: encoded) {
            alignedRow.push_back(make_pair(AlignedExtendedBase(e), runLength));
        }
        alignment.push_back(alignedRow);
        weights.push_back(1);
    }
    for(const AlignedExtendedSequence& row: alignment) {
        SHASTA2_ASSERT(toString(row) == toString(alignment.front()));
    }
    cout << "All " << alignment.size() << " sequences encode to the same " <<
        alignment.front().size() << " symbols." << endl;



    // Vote, using the default estimator.
    vector< pair<Base, uint64_t> > consensus;
    AlignedExtendedSequence alignedConsensus;
    extendedConsensus(alignment, weights, RunLengthEstimator::Mode,
        consensus, alignedConsensus);
    const string consensusString = msa1ToString(consensus);

    // The expected consensus. The two variable A runs come out at 12 and 11.
    //
    // These are the KNOWN TRUE run lengths at this locus, so this is a check
    // against ground truth and not merely against what some other tool produces.
    // The median gives 11 and 11 here and is wrong: 10 of the 19 reads report the
    // first run as shorter than 12, so the median lands one base low even though
    // 12 is the single most common length. See RunLengthEstimator.
    const string expectedConsensus =
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT"
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT";
    cout << "Consensus length " << consensusString.size() << "." << endl;
    SHASTA2_ASSERT(consensusString == expectedConsensus);

    // The consensus is the aligned consensus expanded.
    {
        string fromAligned;
        for(const auto& [e, runLength]: alignedConsensus) {
            if(e.isGap()) {
                SHASTA2_ASSERT(runLength == 0);
                continue;
            }
            for(uint64_t k=0; k<runLength; k++) {
                fromAligned.push_back(e.base().character());
            }
        }
        SHASTA2_ASSERT(fromAligned == consensusString);
    }

    // Coverage is the weight supporting the chosen length, and can never exceed
    // the total weight.
    {
        uint64_t totalWeight = 0;
        for(const uint64_t weight: weights) {
            totalWeight += weight;
        }
        for(const auto& [base, coverage]: consensus) {
            SHASTA2_ASSERT(coverage > 0);
            SHASTA2_ASSERT(coverage <= totalWeight);
        }
    }



    // The median gives a different, and here incorrect, answer at the first run:
    // 11 rather than the true 12, in each of the two tandem copies. Asserted so
    // the two estimators cannot silently converge, and so the direction of the
    // difference is recorded.
    {
        vector< pair<Base, uint64_t> > medianConsensus;
        AlignedExtendedSequence ignore;
        extendedConsensus(alignment, weights, RunLengthEstimator::Median,
            medianConsensus, ignore);
        const string medianString = msa1ToString(medianConsensus);
        cout << "Mode gives length " << consensusString.size() <<
            " (correct), median gives length " << medianString.size() << "." << endl;
        SHASTA2_ASSERT(medianString != consensusString);
        SHASTA2_ASSERT(medianString.size() == consensusString.size() - 2);
    }



    // Expanding must preserve every read and give every row the same length as
    // the aligned consensus.
    {
        vector< vector<AlignedBase> > expandedAlignment;
        vector<AlignedBase> expandedAlignedConsensus;
        expandExtendedAlignment(alignment, alignedConsensus,
            expandedAlignment, expandedAlignedConsensus);

        SHASTA2_ASSERT(expandedAlignment.size() == msa1BadAlignmentRows.size());
        vector<string> expandedRows;
        for(uint64_t i=0; i<expandedAlignment.size(); i++) {
            SHASTA2_ASSERT(expandedAlignment[i].size() == expandedAlignedConsensus.size());

            string ungapped;
            string full;
            for(const AlignedBase b: expandedAlignment[i]) {
                full.push_back(b.character());
                if(not b.isGap()) {
                    ungapped.push_back(b.character());
                }
            }
            expandedRows.push_back(full);

            string expected;
            for(const char c: msa1BadAlignmentRows[i]) {
                if(c != '-') {
                    expected.push_back(c);
                }
            }
            SHASTA2_ASSERT(ungapped == expected);
        }

        // No column of the expanded alignment holds two different bases. The
        // input had 4. This is structural: a poly symbol and the base next to it
        // are different symbols and get different columns, so they cannot be
        // confused.
        const uint64_t impureColumnCount = msa1ImpureColumnCount(expandedRows);
        cout << "The expanded alignment has " << impureColumnCount <<
            " columns containing more than one base, in " <<
            expandedAlignedConsensus.size() << " columns." << endl;
        SHASTA2_ASSERT(impureColumnCount == 0);
    }



    // Weights must matter, for both estimators. Give the two reads whose first
    // run is 10 and whose second is 12 overwhelming weight, and the consensus
    // must follow them rather than the other seventeen reads.
    {
        vector<uint64_t> skewedWeights = weights;
        skewedWeights[9] = 100;
        skewedWeights[10] = 100;

        for(const RunLengthEstimator estimator:
            {RunLengthEstimator::Median, RunLengthEstimator::Mode}) {
            vector< pair<Base, uint64_t> > skewedConsensus;
            AlignedExtendedSequence ignore;
            extendedConsensus(alignment, skewedWeights, estimator,
                skewedConsensus, ignore);
            const string expected =
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT" +
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT";
            SHASTA2_ASSERT(msa1ToString(skewedConsensus) == expected);
        }
    }



    // A real substitution must survive: the majority allele is kept, not
    // absorbed or dropped, and the consensus follows when the split moves.
    {
        vector<AlignedExtendedSequence> snpAlignment = alignment;

        uint64_t column = invalid<uint64_t>;
        for(uint64_t j=0; j<snpAlignment.front().size(); j++) {
            if(not snpAlignment.front()[j].first.isPoly()) {
                column = j;
                break;
            }
        }
        SHASTA2_ASSERT(column != invalid<uint64_t>);
        const AlignedExtendedBase original = snpAlignment.front()[column].first;
        const AlignedExtendedBase mutated =
            AlignedExtendedBase(ExtendedBase::fromBase(original.base().complement()));
        SHASTA2_ASSERT(mutated != original);

        for(uint64_t i=0; i<8; i++) {
            snpAlignment[i][column].first = mutated;
        }
        vector< pair<Base, uint64_t> > snpConsensus;
        AlignedExtendedSequence snpAlignedConsensus;
        extendedConsensus(snpAlignment, weights, RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus);
        SHASTA2_ASSERT(snpAlignedConsensus[column].first == original);
        SHASTA2_ASSERT(snpAlignedConsensus[column].second == 1);

        uint64_t positionInConsensus = 0;
        for(uint64_t j=0; j<column; j++) {
            positionInConsensus += snpAlignedConsensus[j].second;
        }
        SHASTA2_ASSERT(snpConsensus[positionInConsensus].second == 11);
        cout << "With a substitution in 8 of 19 rows, the majority allele is kept "
            "with coverage " << snpConsensus[positionInConsensus].second << "." << endl;

        // Flip it in 12 rows and the consensus must follow.
        for(uint64_t i=0; i<12; i++) {
            snpAlignment[i][column].first = mutated;
        }
        extendedConsensus(snpAlignment, weights, RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus);
        SHASTA2_ASSERT(snpAlignedConsensus[column].first == mutated);
    }



    // A column that is a gap in a majority of rows becomes a gap and contributes
    // nothing to the consensus.
    {
        vector<AlignedExtendedSequence> gapAlignment = alignment;
        for(AlignedExtendedSequence& row: gapAlignment) {
            row.insert(row.begin(), make_pair(AlignedExtendedBase::gap(), 0UL));
        }
        gapAlignment[0][0] = make_pair(
            AlignedExtendedBase::fromCharacter('C'), 1UL);

        vector< pair<Base, uint64_t> > gapConsensus;
        AlignedExtendedSequence gapAlignedConsensus;
        extendedConsensus(gapAlignment, weights, RunLengthEstimator::Mode,
            gapConsensus, gapAlignedConsensus);
        SHASTA2_ASSERT(gapAlignedConsensus[0].first.isGap());
        SHASTA2_ASSERT(gapAlignedConsensus[0].second == 0);
        SHASTA2_ASSERT(msa1ToString(gapConsensus) == consensusString);
    }



    // No regression when nothing is long enough to become a poly symbol: the
    // result must equal plain weighted column majority.
    {
        const vector<string> rows = {"ACGTACGTAC", "ACGTACGTAC", "ACGTTCGTAC"};
        vector<AlignedExtendedSequence> a;
        for(const string& row: rows) {
            AlignedExtendedSequence r;
            for(const char c: row) {
                r.push_back(make_pair(AlignedExtendedBase::fromCharacter(c), 1UL));
                SHASTA2_ASSERT(not r.back().first.isPoly());
            }
            a.push_back(r);
        }
        vector< pair<Base, uint64_t> > c;
        AlignedExtendedSequence ac;
        extendedConsensus(a, {1, 1, 1}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(msa1ToString(c) == "ACGTACGTAC");
        SHASTA2_ASSERT(c[4].second == 2);
        SHASTA2_ASSERT(c[0].second == 3);
    }



    // attachRunLengths puts the lengths back on the columns of a row.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "CTCaGa");

        const vector<AlignedExtendedBase> row =
            vectorOfAlignedExtendedBasesFromString("C-TC-aG-a");
        AlignedExtendedSequence attached;
        attachRunLengths(row, encoded, Anchoring::BothSides, attached);
        const vector<uint64_t> expected = {1, 0, 1, 1, 0, 12, 1, 0, 11};
        SHASTA2_ASSERT(attached.size() == expected.size());
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
        }
        SHASTA2_ASSERT(toString(attached) == "C-TC-aG-a");
    }



    // A row can cover only part of its encoding, because theseus trims the
    // overhang of a sequence fixed on one side only.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A') + string(7, 'T')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "CTCaGat");

        // Left fixed: the trailing poly symbol is trimmed.
        AlignedExtendedSequence attached;
        attachRunLengths(vectorOfAlignedExtendedBasesFromString("CTC-aGa"),
            encoded, Anchoring::LeftOnly, attached);
        const vector<uint64_t> expected = {1, 1, 1, 0, 12, 1, 11};
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
        }
    }
    {
        // Right fixed: the leading poly symbol is trimmed, so the lengths must
        // come from offset 1, not 0.
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            string(9, 'G') + "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "gCTCaGa");

        AlignedExtendedSequence attached;
        attachRunLengths(vectorOfAlignedExtendedBasesFromString("CTCaG-a"),
            encoded, Anchoring::RightOnly, attached);
        const vector<uint64_t> expected = {1, 1, 1, 12, 1, 0, 11};
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
            // The trimmed poly symbol's length of 9 must not appear anywhere.
            SHASTA2_ASSERT(attached[i].second != 9);
        }
    }



    // The window must come from the anchoring, not from searching the encoding
    // for the row's symbols.
    //
    // Here the encoding is a tandem repeat, A a C A a C, whose two poly runs have
    // different lengths, 10 and 20. The sequence is fixed on the right, so
    // theseus trimmed the first copy and the row is the second one, whose run is
    // 20. Searching for the row's symbols finds the FIRST copy and would pair the
    // row with 10 instead. Tandem repeats are exactly what this code is for, and
    // the 19 read locus it was developed on is itself a 2x tandem duplication.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "A" + string(10, 'G') + "C" + "A" + string(20, 'G') + "C"), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "AgCAgC");
        SHASTA2_ASSERT(encoded[1].second == 10);
        SHASTA2_ASSERT(encoded[4].second == 20);

        const vector<AlignedExtendedBase> row =
            vectorOfAlignedExtendedBasesFromString("AgC");

        AlignedExtendedSequence attached;
        attachRunLengths(row, encoded, Anchoring::RightOnly, attached);
        SHASTA2_ASSERT(attached[1].second == 20);

        // And the mirror image: fixed on the left, so the row is the first copy.
        attachRunLengths(row, encoded, Anchoring::LeftOnly, attached);
        SHASTA2_ASSERT(attached[1].second == 10);
    }



    // Degenerate inputs.
    {
        vector< pair<Base, uint64_t> > c;
        AlignedExtendedSequence ac;

        extendedConsensus({}, {}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.empty());

        // A single sequence: the consensus is that sequence.
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        AlignedExtendedSequence one;
        for(const auto& [e, runLength]: encoded) {
            one.push_back(make_pair(AlignedExtendedBase(e), runLength));
        }
        extendedConsensus({one}, {7}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(msa1ToString(c) == "CTC" + string(12, 'A') + "G" + string(11, 'A'));
        for(const auto& [base, coverage]: c) {
            SHASTA2_ASSERT(coverage == 7);
        }

        // Every row a gap at every column.
        const AlignedExtendedSequence allGapRow = {
            make_pair(AlignedExtendedBase::gap(), 0UL),
            make_pair(AlignedExtendedBase::gap(), 0UL)};
        extendedConsensus({allGapRow, allGapRow}, {1, 1},
            RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.size() == 2);
        SHASTA2_ASSERT(ac[0].first.isGap() and ac[1].first.isGap());
    }

    cout << "testMsa1Consensus passed." << endl;
}
