// Shasta2.
#include "msa1.hpp"
#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "theseusWrapper.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"



// The look up table used by ExtendedBase::fromCharacter.
// Note this differs from BaseInitializer and AlignedBaseInitializer, which map
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
    table['-'] = ExtendedBase::gapValue;
}



// See msa1.hpp for comments.
void shasta2::encodeExtended(
    const vector<Base>& sequence,
    uint64_t threshold,
    vector<ExtendedBase>& encoded,
    vector<uint64_t>& runLengths)
{
    // A threshold of 0 would make every run, including a run of length 1,
    // encode to a bare poly symbol. That loses the length of every run and is
    // exactly the plain RLE behavior this alphabet exists to avoid.
    SHASTA2_ASSERT(threshold > 0);

    encoded.clear();
    runLengths.clear();

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
            // run. Its length is kept in runLengths and plays no part in any
            // alignment done on the symbols.
            encoded.push_back(ExtendedBase::polyFromBase(base));
            runLengths.push_back(runLength);

        } else {

            // A short run keeps one plain symbol per base, so its length is
            // still visible to an aligner. Short run lengths are reliable and
            // are real signal, unlike long ones.
            for(uint64_t i=0; i<runLength; i++) {
                encoded.push_back(ExtendedBase::fromBase(base));
                runLengths.push_back(1);
            }
        }

        runBegin = runEnd;
    }

    SHASTA2_ASSERT(encoded.size() == runLengths.size());
}



// See msa1.hpp for comments.
void shasta2::decodeExtended(
    const vector<ExtendedBase>& encoded,
    const vector<uint64_t>& runLengths,
    vector<Base>& sequence)
{
    SHASTA2_ASSERT(encoded.size() == runLengths.size());

    sequence.clear();
    for(uint64_t i=0; i<encoded.size(); i++) {
        const ExtendedBase extendedBase = encoded[i];
        if(extendedBase.isGap()) {
            continue;
        }
        const Base base = extendedBase.base();
        for(uint64_t j=0; j<runLengths[i]; j++) {
            sequence.push_back(base);
        }
    }
}



// See msa1.hpp for comments.
void shasta2::alignedRunLengths(
    const vector<ExtendedBase>& alignedRow,
    const vector<ExtendedBase>& encoded,
    const vector<uint64_t>& encodedRunLengths,
    Anchoring anchoring,
    vector<uint64_t>& alignedRunLengthsArgument)
{
    SHASTA2_ASSERT(encoded.size() == encodedRunLengths.size());

    // Remove the gaps from the aligned row.
    vector<ExtendedBase> rowSymbols;
    for(const ExtendedBase e: alignedRow) {
        if(not e.isGap()) {
            rowSymbols.push_back(e);
        }
    }
    SHASTA2_ASSERT(rowSymbols.size() <= encoded.size());

    // Locate the part of the encoding that this row covers.
    //
    // For a sequence fixed on both sides the row covers all of it. For a
    // sequence fixed on one side only, theseus trims the overhang off the
    // unanchored end, so the row is a prefix or a suffix of the encoding.
    // Which one is known from the anchoring and is not inferred: see the
    // comment in msa1.hpp for why searching the encoding for the row's symbols
    // gives the wrong answer on a repeat.
    uint64_t offset = 0;
    switch(anchoring) {
    case Anchoring::BothSides:
        SHASTA2_ASSERT(rowSymbols.size() == encoded.size());
        offset = 0;
        break;
    case Anchoring::LeftOnly:
        // Anchored on the left, trimmed on the right: a prefix.
        offset = 0;
        break;
    case Anchoring::RightOnly:
        // Anchored on the right, trimmed on the left: a suffix.
        offset = encoded.size() - rowSymbols.size();
        break;
    }

    alignedRunLengthsArgument.clear();
    alignedRunLengthsArgument.resize(alignedRow.size(), 0);

    uint64_t i = 0;
    for(uint64_t j=0; j<alignedRow.size(); j++) {
        if(alignedRow[j].isGap()) {
            continue;
        }
        SHASTA2_ASSERT(alignedRow[j] == encoded[offset + i]);
        alignedRunLengthsArgument[j] = encodedRunLengths[offset + i];
        ++i;
    }
    SHASTA2_ASSERT(i == rowSymbols.size());
}



// See msa1.hpp for comments.
void shasta2::extendedConsensus(
    const vector< vector<ExtendedBase> >& alignment,
    const vector< vector<uint64_t> >& runLengths,
    const vector<uint64_t>& weights,
    RunLengthEstimator estimator,
    vector< pair<Base, uint64_t> >& consensus,
    vector<ExtendedBase>& alignedExtendedConsensus,
    vector<uint64_t>& consensusRunLengths)
{
    consensus.clear();
    alignedExtendedConsensus.clear();
    consensusRunLengths.clear();

    const uint64_t n = alignment.size();
    if(n == 0) {
        return;
    }
    SHASTA2_ASSERT(runLengths.size() == n);
    SHASTA2_ASSERT(weights.size() == n);

    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=0; i<n; i++) {
        SHASTA2_ASSERT(alignment[i].size() == alignmentLength);
        SHASTA2_ASSERT(runLengths[i].size() == alignmentLength);
    }

    alignedExtendedConsensus.resize(alignmentLength);
    consensusRunLengths.resize(alignmentLength, 0);

    // Loop over alignment columns.
    for(uint64_t j=0; j<alignmentLength; j++) {

        // Total weight for each symbol at this column, including the gap.
        array<uint64_t, ExtendedBase::gapValue + 1> symbolWeight;
        fill(symbolWeight.begin(), symbolWeight.end(), 0);
        for(uint64_t i=0; i<n; i++) {
            symbolWeight[alignment[i][j].value] += weights[i];
        }

        // The consensus symbol is the one with the most weight.
        const auto it = std::ranges::max_element(symbolWeight);
        const ExtendedBase consensusSymbol =
            ExtendedBase::fromInteger(uint64_t(it - symbolWeight.begin()));
        alignedExtendedConsensus[j] = consensusSymbol;

        if(consensusSymbol.isGap()) {
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
            // length is both simpler and faster than a map.
            uint64_t maxObserved = 0;
            for(uint64_t i=0; i<n; i++) {
                if(alignment[i][j] == consensusSymbol) {
                    maxObserved = max(maxObserved, runLengths[i][j]);
                }
            }
            vector<uint64_t> lengthWeight(maxObserved + 1, 0);
            uint64_t totalWeight = 0;
            for(uint64_t i=0; i<n; i++) {
                if(alignment[i][j] == consensusSymbol) {
                    lengthWeight[runLengths[i][j]] += weights[i];
                    totalWeight += weights[i];
                }
            }
            SHASTA2_ASSERT(totalWeight > 0);

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

        consensusRunLengths[j] = consensusRunLength;

        // Store the run in the consensus.
        const Base base = consensusSymbol.base();
        for(uint64_t k=0; k<consensusRunLength; k++) {
            consensus.push_back(make_pair(base, coverage));
        }
    }
}



// See msa1.hpp for comments.
void shasta2::expandExtendedAlignment(
    const vector< vector<ExtendedBase> >& alignment,
    const vector< vector<uint64_t> >& runLengths,
    const vector<ExtendedBase>& alignedExtendedConsensus,
    const vector<uint64_t>& consensusRunLengths,
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
    SHASTA2_ASSERT(alignedExtendedConsensus.size() == alignmentLength);
    SHASTA2_ASSERT(consensusRunLengths.size() == alignmentLength);

    // The width of each column after expansion: the longest run seen there,
    // including the consensus, and at least 1.
    vector<uint64_t> width(alignmentLength, 1);
    for(uint64_t j=0; j<alignmentLength; j++) {
        for(uint64_t i=0; i<n; i++) {
            width[j] = max(width[j], runLengths[i][j]);
        }
        width[j] = max(width[j], consensusRunLengths[j]);
    }

    // Expand each row, padding on the right with gaps.
    const auto expand = [&width, alignmentLength](
        const vector<ExtendedBase>& row,
        const vector<uint64_t>& rowRunLengths,
        vector<AlignedBase>& expanded)
    {
        expanded.clear();
        for(uint64_t j=0; j<alignmentLength; j++) {
            const ExtendedBase e = row[j];
            const uint64_t runLength = e.isGap() ? 0 : rowRunLengths[j];
            for(uint64_t k=0; k<runLength; k++) {
                expanded.push_back(AlignedBase(e.base()));
            }
            for(uint64_t k=runLength; k<width[j]; k++) {
                expanded.push_back(AlignedBase::gap());
            }
        }
    };

    expandedAlignment.resize(n);
    for(uint64_t i=0; i<n; i++) {
        expand(alignment[i], runLengths[i], expandedAlignment[i]);
    }
    expand(alignedExtendedConsensus, consensusRunLengths, expandedAlignedConsensus);

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
    // A homopolymer run longer than the threshold becomes a single poly symbol,
    // and its length moves into runLengths, where the aligner cannot see it.
    // This is what stops a run length difference from being turned into a
    // mismatch, which is what misplaces the base bordering a long run.
    vector< pair<vector<ExtendedBase>, uint64_t> > encodedFixed;
    vector< pair<vector<ExtendedBase>, uint64_t> > encodedLeftFixed;
    vector< pair<vector<ExtendedBase>, uint64_t> > encodedRightFixed;

    // The encodings, their run lengths, and how each is anchored, flattened in
    // the same order the rows will come back from the aligner: fixed, then left
    // fixed, then right fixed. The anchoring is carried along because it says
    // which end theseus may trim, which is needed to line the run lengths up
    // with the alignment columns.
    vector< vector<ExtendedBase> > encodings;
    vector< vector<uint64_t> > encodedRunLengths;
    vector<Anchoring> anchorings;
    vector<uint64_t> weights;

    const auto encodeGroup = [&](
        const vector< pair<vector<Base>, uint64_t> >& in,
        Anchoring anchoring,
        vector< pair<vector<ExtendedBase>, uint64_t> >& out)
    {
        out.clear();
        for(const auto& [sequence, weight]: in) {
            vector<ExtendedBase> encoded;
            vector<uint64_t> runLengths;
            encodeExtended(sequence, threshold, encoded, runLengths);
            out.push_back(make_pair(encoded, weight));
            encodings.push_back(encoded);
            encodedRunLengths.push_back(runLengths);
            anchorings.push_back(anchoring);
            weights.push_back(weight);
        }
    };
    encodeGroup(fixedSequences, Anchoring::BothSides, encodedFixed);
    encodeGroup(leftFixedSequences, Anchoring::LeftOnly, encodedLeftFixed);
    encodeGroup(rightFixedSequences, Anchoring::RightOnly, encodedRightFixed);

    // Align the encoded sequences.
    vector< vector<ExtendedBase> > extendedAlignment;
    theseusExtended(encodedFixed, encodedLeftFixed, encodedRightFixed,
        extendedAlignment);
    SHASTA2_ASSERT(extendedAlignment.size() == encodedRunLengths.size());

    // Line up each row's run lengths with the columns of the alignment.
    // The encoding passed here is the one we handed to the aligner, not one
    // recovered from the row: for a sequence fixed on one side only, theseus
    // trims the overhang, so the row covers only part of the encoding, and the
    // anchoring says which part.
    vector< vector<uint64_t> > alignedRunLengthsAllRows(extendedAlignment.size());
    for(uint64_t i=0; i<extendedAlignment.size(); i++) {
        alignedRunLengths(extendedAlignment[i], encodings[i], encodedRunLengths[i],
            anchorings[i], alignedRunLengthsAllRows[i]);
    }

    // Vote, column by column, with long run lengths voted separately.
    vector<ExtendedBase> alignedExtendedConsensus;
    vector<uint64_t> consensusRunLengths;
    extendedConsensus(extendedAlignment, alignedRunLengthsAllRows, weights,
        estimator, consensus, alignedExtendedConsensus, consensusRunLengths);

    // Expand back to plain bases. This is done even when the alignment is not
    // requested, because alignedConsensus is always computed and has to have the
    // same length as the alignment rows.
    vector< vector<AlignedBase> > expandedAlignment;
    expandExtendedAlignment(extendedAlignment, alignedRunLengthsAllRows,
        alignedExtendedConsensus, consensusRunLengths,
        expandedAlignment, alignedConsensus);

    if(computeAlignment) {
        alignment = expandedAlignment;
    } else {
        alignment.clear();
    }
}


// Test the extended alphabet and its encoding.
void shasta2::testMsa1ExtendedBase()
{
    // The 19 bubble arm sequences that motivated this alphabet.
    // They are reads over a single locus with structure
    // CTC A{n1} G A{n2} GTT, inside a 2x tandem duplication.
    // Only n1 and n2 vary between them.
    const vector<string> testSequences = {
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
        SHASTA2_ASSERT(not plain.isGap());
        SHASTA2_ASSERT(not poly.isGap());
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
    }
    const ExtendedBase gap = ExtendedBase::gap();
    SHASTA2_ASSERT(gap.isGap());
    SHASTA2_ASSERT(not gap.isPoly());
    SHASTA2_ASSERT(gap.isValid());
    SHASTA2_ASSERT(gap.complement().isGap());
    SHASTA2_ASSERT(gap.character() == '-');

    // character() and fromCharacter() round trip over all 9 symbols.
    for(uint64_t value=0; value<=ExtendedBase::gapValue; value++) {
        const ExtendedBase e = ExtendedBase::fromInteger(value);
        SHASTA2_ASSERT(ExtendedBase::fromCharacter(e.character()) == e);
    }

    // An invalid character throws.
    bool threw = false;
    try {
        ExtendedBase::fromCharacter('X');
    } catch(const std::exception&) {
        threw = true;
    }
    SHASTA2_ASSERT(threw);

    // toString() and vectorOfExtendedBasesFromString() round trip.
    {
        const string s = "ACGTacgt-";
        SHASTA2_ASSERT(toString(vectorOfExtendedBasesFromString(s)) == s);
    }



    // Round trip all 19 sequences, at a range of thresholds.
    // Also check that they all encode to the SAME symbol string. That is the
    // property that makes the homopolymer ambiguity structurally impossible:
    // if every read has the same run structure there is nothing to align, and
    // no scoring choice can misplace the G between the two A runs.
    for(uint64_t threshold=1; threshold<=12; threshold++) {
        vector<ExtendedBase> firstEncoded;
        vector<uint64_t> firstRunLengths;

        // Set to false if the sequences do not all encode to the same symbol string.
        bool uniform = true;

        // The positions at which runLengths differ between sequences.
        vector<uint64_t> variablePositions;

        for(uint64_t i=0; i<testSequences.size(); i++) {
            const vector<Base> sequence = vectorOfBasesFromString(testSequences[i]);

            vector<ExtendedBase> encoded;
            vector<uint64_t> runLengths;
            encodeExtended(sequence, threshold, encoded, runLengths);

            // The encoding is lossless.
            vector<Base> decoded;
            decodeExtended(encoded, runLengths, decoded);
            SHASTA2_ASSERT(decoded == sequence);

            // A plain symbol always stands for exactly one base, and a poly
            // symbol for at least one.
            for(uint64_t j=0; j<encoded.size(); j++) {
                if(encoded[j].isPoly()) {
                    SHASTA2_ASSERT(runLengths[j] >= 1);
                } else {
                    SHASTA2_ASSERT(runLengths[j] == 1);
                }
            }

            // A poly symbol stands for a whole run on its own, so it is never
            // adjacent to a symbol with the same base.
            for(uint64_t j=1; j<encoded.size(); j++) {
                if(encoded[j].isPoly() or encoded[j-1].isPoly()) {
                    SHASTA2_ASSERT(encoded[j].base() != encoded[j-1].base());
                }
            }

            if(i == 0) {
                firstEncoded = encoded;
                firstRunLengths = runLengths;
            } else if(encoded != firstEncoded) {
                uniform = false;
            } else {
                SHASTA2_ASSERT(runLengths.size() == firstRunLengths.size());
                for(uint64_t j=0; j<runLengths.size(); j++) {
                    if(runLengths[j] != firstRunLengths[j]) {
                        if(find(variablePositions.begin(), variablePositions.end(), j) ==
                            variablePositions.end()) {
                            variablePositions.push_back(j);
                        }
                    }
                }
            }
        }

        cout << "Threshold " << threshold << ": " <<
            (uniform ? "all 19 sequences encode to the same " : "encodings DIFFER, first has ") <<
            firstEncoded.size() << " symbols";
        if(uniform) {
            cout << ", " << variablePositions.size() << " of which have variable length";
        }
        cout << "." << endl;

        // The shortest A run in this set is 6 bases (sequence 16), and the
        // others at the same two positions are 7 or longer.
        //
        // At a threshold below 6 every one of those runs is longer than the
        // threshold, so all of them end in a poly symbol and all 19 sequences
        // encode to the SAME symbol string, differing only in runLengths at
        // exactly 4 positions: the two A runs, in each of the two tandem
        // copies. That is the property that makes the homopolymer ambiguity
        // structurally impossible here. There is nothing left to align, so no
        // scoring choice can misplace the G between the two runs.
        //
        // At a threshold of 6 the run of 6 is no longer above the threshold, so
        // it encodes as 6 plain symbols with no poly symbol while the others
        // encode as 6 plain symbols plus a poly symbol. The symbol strings
        // differ by one symbol and the guarantee is lost.
        //
        // This is the boundary effect that makes a fixed threshold inadequate,
        // and it happens at exactly the default of 6. It is asserted here in
        // both directions so that a later adaptive threshold has a concrete
        // case to satisfy.
        if(threshold < 6) {
            SHASTA2_ASSERT(uniform);
            SHASTA2_ASSERT(variablePositions.size() == 4);
        }
        if(threshold == 6) {
            SHASTA2_ASSERT(not uniform);
        }

        // The default must be one of the values that keeps the encoding uniform
        // here. This is the property the whole alphabet exists to provide, so a
        // change to the default that loses it should not pass silently.
        if(threshold == defaultHomopolymerThreshold) {
            SHASTA2_ASSERT(uniform);
        }
    }



    // A sequence with no run longer than the threshold encodes to itself, one
    // plain symbol per base. So behavior on sequences without long homopolymers
    // is unchanged by this alphabet.
    {
        const vector<Base> sequence = vectorOfBasesFromString("ACGTACGTAACCGGTT");
        vector<ExtendedBase> encoded;
        vector<uint64_t> runLengths;
        encodeExtended(sequence, defaultHomopolymerThreshold, encoded, runLengths);
        SHASTA2_ASSERT(encoded.size() == sequence.size());
        for(uint64_t i=0; i<encoded.size(); i++) {
            SHASTA2_ASSERT(not encoded[i].isPoly());
            SHASTA2_ASSERT(encoded[i].base() == sequence[i]);
            SHASTA2_ASSERT(runLengths[i] == 1);
        }
    }



    // The length attribute of a poly symbol must not take part in the alignment.
    //
    // This is why the length lives in runLengths and not inside ExtendedBase.
    // An ExtendedBase holds nothing but its symbol value, so two poly symbols
    // for the same base always compare equal however long the runs they stand
    // for. An aligner given the encoded vector therefore cannot see the lengths
    // at all, and cannot turn a run length difference into a mismatch. The
    // lengths are read back only at the end, to decide the consensus length of
    // each run.
    //
    // Asserted here because it is the property the whole alphabet exists for,
    // and it would be silently lost if ExtendedBase ever gained a length member.
    {
        const uint64_t threshold = 5;
        const vector<string> sameStructure = {
            "CTCAAAAAAAAAAAAGAAAAAAAAAAAGTT",   // A run of 12, then 11
            "CTCAAAAAAAAAAGAAAAAAAAAAAAGTT",    // A run of 10, then 12
            "CTCAAAAAAGAAAAAAAAAAGTT"           // A run of  6, then 10
        };

        vector<ExtendedBase> firstEncoded;
        for(uint64_t i=0; i<sameStructure.size(); i++) {
            vector<ExtendedBase> encoded;
            vector<uint64_t> runLengths;
            encodeExtended(vectorOfBasesFromString(sameStructure[i]), threshold,
                encoded, runLengths);

            // CTC polyA G polyA GTT.
            SHASTA2_ASSERT(toString(encoded) == "CTCaGaGTT");

            if(i == 0) {
                firstEncoded = encoded;
            } else {
                // The symbols are identical, so an aligner sees three copies of
                // the same thing. The run lengths differ, and are kept.
                SHASTA2_ASSERT(encoded == firstEncoded);
            }
        }

        // Two poly symbols for the same base are equal, whatever run they came from.
        vector<ExtendedBase> shortRun;
        vector<ExtendedBase> longRun;
        vector<uint64_t> shortRunLengths;
        vector<uint64_t> longRunLengths;
        encodeExtended(vector<Base>(6, Base::fromCharacter('A')), threshold,
            shortRun, shortRunLengths);
        encodeExtended(vector<Base>(600, Base::fromCharacter('A')), threshold,
            longRun, longRunLengths);
        SHASTA2_ASSERT(shortRun == longRun);
        SHASTA2_ASSERT(shortRun.size() == 1);
        SHASTA2_ASSERT(shortRunLengths[0] == 6);
        SHASTA2_ASSERT(longRunLengths[0] == 600);
    }



    // Behavior of a single run as its length crosses the threshold.
    //
    // A run of length threshold encodes to threshold plain symbols, and a run
    // of length threshold+1 encodes to a single poly symbol. So the encoding is
    // discontinuous at the threshold: one extra base changes the symbol string
    // by threshold-1 symbols, and to an aligner working on the symbols a run
    // that straddles the threshold looks like a large indel rather than a small
    // one. That is the price of collapsing a whole run to one symbol, and it is
    // the reason the threshold has to be chosen with care. It is asserted here
    // so the behavior is pinned down rather than discovered later.
    {
        const uint64_t threshold = defaultHomopolymerThreshold;
        for(uint64_t runLength=1; runLength<=threshold+3; runLength++) {
            const vector<Base> sequence(runLength, Base::fromCharacter('A'));
            vector<ExtendedBase> encoded;
            vector<uint64_t> runLengths;
            encodeExtended(sequence, threshold, encoded, runLengths);

            if(runLength > threshold) {
                // One poly symbol standing for the whole run, whatever its length.
                SHASTA2_ASSERT(encoded.size() == 1);
                SHASTA2_ASSERT(encoded[0].isPoly());
                SHASTA2_ASSERT(runLengths[0] == runLength);
            } else {
                // One plain symbol per base.
                SHASTA2_ASSERT(encoded.size() == runLength);
                for(uint64_t i=0; i<encoded.size(); i++) {
                    SHASTA2_ASSERT(not encoded[i].isPoly());
                }
            }

            vector<Base> decoded;
            decodeExtended(encoded, runLengths, decoded);
            SHASTA2_ASSERT(decoded == sequence);
        }

        // All runs longer than the threshold encode to the same single symbol,
        // regardless of length. This is the property that removes long run
        // lengths from the alignment entirely.
        vector<ExtendedBase> encodedShort;
        vector<ExtendedBase> encodedLong;
        vector<uint64_t> ignore;
        encodeExtended(vector<Base>(threshold + 1, Base::fromCharacter('A')),
            threshold, encodedShort, ignore);
        encodeExtended(vector<Base>(1000, Base::fromCharacter('A')),
            threshold, encodedLong, ignore);
        SHASTA2_ASSERT(encodedShort == encodedLong);
        SHASTA2_ASSERT(encodedShort.size() == 1);

        // The discontinuity at the threshold.
        vector<ExtendedBase> encodedAtThreshold;
        encodeExtended(vector<Base>(threshold, Base::fromCharacter('A')),
            threshold, encodedAtThreshold, ignore);
        SHASTA2_ASSERT(encodedAtThreshold.size() == threshold);
        SHASTA2_ASSERT(encodedShort.size() == 1);
    }



    // Degenerate inputs.
    {
        vector<ExtendedBase> encoded;
        vector<uint64_t> runLengths;
        vector<Base> decoded;

        // Empty sequence.
        encodeExtended(vector<Base>(), defaultHomopolymerThreshold, encoded, runLengths);
        SHASTA2_ASSERT(encoded.empty());
        SHASTA2_ASSERT(runLengths.empty());
        decodeExtended(encoded, runLengths, decoded);
        SHASTA2_ASSERT(decoded.empty());

        // Single base.
        const vector<Base> single = vectorOfBasesFromString("G");
        encodeExtended(single, defaultHomopolymerThreshold, encoded, runLengths);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(not encoded[0].isPoly());
        decodeExtended(encoded, runLengths, decoded);
        SHASTA2_ASSERT(decoded == single);

        // One run spanning the whole sequence: a single poly symbol.
        const vector<Base> oneRun(100, Base::fromCharacter('T'));
        encodeExtended(oneRun, defaultHomopolymerThreshold, encoded, runLengths);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(encoded[0].isPoly());
        SHASTA2_ASSERT(runLengths[0] == oneRun.size());
        decodeExtended(encoded, runLengths, decoded);
        SHASTA2_ASSERT(decoded == oneRun);

        // Gaps in the encoded sequence are skipped by decodeExtended.
        encoded.push_back(ExtendedBase::gap());
        runLengths.push_back(1);
        vector<Base> decodedWithGap;
        decodeExtended(encoded, runLengths, decodedWithGap);
        SHASTA2_ASSERT(decodedWithGap == oneRun);
    }



    // Round trip random sequences over a range of thresholds.
    // A simple deterministic generator, so a failure is reproducible.
    // Runs are drawn with a length distribution weighted toward short runs,
    // matching what the ONT reads actually look like, but reaching well past
    // any threshold tested.
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
                    // Make sure consecutive runs have different bases, so the
                    // run structure is what we intend it to be.
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

                vector<ExtendedBase> encoded;
                vector<uint64_t> runLengths;
                encodeExtended(sequence, threshold, encoded, runLengths);

                vector<Base> decoded;
                decodeExtended(encoded, runLengths, decoded);
                SHASTA2_ASSERT(decoded == sequence);

                // The encoding never grows the sequence.
                SHASTA2_ASSERT(encoded.size() <= max(sequence.size(), 1UL) );
            }
        }
    }

    cout << "testMsa1ExtendedBase passed." << endl;
}


// Test the consensus computed from an alignment over the extended alphabet.
//
// Theseus is not exercised here: extendedConsensus, alignedRunLengths and
// expandExtendedAlignment are pure functions, so the alignment is supplied
// directly. That also lets the real, defective alignment be used as input.
void shasta2::testMsa1Consensus()
{
    // The actual alignment abpoa produces for the 19 bubble arm sequences.
    // Four of its columns hold a mixture of A and G, because the G between the
    // two long A runs was placed in two different columns. Column-wise majority
    // voting over this alignment is contaminated by that. The point of this test
    // is that voting on run lengths is not.
    const vector<string> badAlignmentRows = {
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

    // Confirm the input really is defective: at least one column holds two
    // different bases. If this ever stops being true the test has lost its point.
    {
        const uint64_t alignmentLength = badAlignmentRows.front().size();
        uint64_t impureColumnCount = 0;
        for(uint64_t j=0; j<alignmentLength; j++) {
            array<bool, 4> seen = {false, false, false, false};
            for(const string& row: badAlignmentRows) {
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
        cout << "The input alignment has " << impureColumnCount <<
            " columns containing more than one base." << endl;
        SHASTA2_ASSERT(impureColumnCount == 4);
    }



    // Recover the sequences from the alignment, encode them, and build an
    // alignment over the extended alphabet by encoding each row in isolation.
    // The rows all encode to the same symbol string here, so they line up
    // without any further work. This stands in for what theseus would return.
    const uint64_t threshold = 5;
    vector< vector<ExtendedBase> > extendedAlignment;
    vector< vector<uint64_t> > alignedRunLengthsAllRows;
    vector<uint64_t> weights;
    for(const string& row: badAlignmentRows) {

        // Strip the gaps to recover the read.
        string ungapped;
        for(const char c: row) {
            if(c != '-') {
                ungapped.push_back(c);
            }
        }
        const vector<Base> sequence = vectorOfBasesFromString(ungapped);

        vector<ExtendedBase> encoded;
        vector<uint64_t> runLengths;
        encodeExtended(sequence, threshold, encoded, runLengths);

        extendedAlignment.push_back(encoded);
        alignedRunLengthsAllRows.push_back(runLengths);
        weights.push_back(1);
    }

    // All 19 encode to the same symbol string, so no alignment is needed.
    // 11 distinct base sequences collapse to 1 distinct encoding.
    for(const vector<ExtendedBase>& row: extendedAlignment) {
        SHASTA2_ASSERT(row == extendedAlignment.front());
    }
    cout << "All " << extendedAlignment.size() <<
        " sequences encode to the same " << extendedAlignment.front().size() <<
        " symbols." << endl;



    // Vote, using the default estimator.
    vector< pair<Base, uint64_t> > consensus;
    vector<ExtendedBase> alignedExtendedConsensus;
    vector<uint64_t> consensusRunLengths;
    extendedConsensus(extendedAlignment, alignedRunLengthsAllRows, weights,
        RunLengthEstimator::Mode,
        consensus, alignedExtendedConsensus, consensusRunLengths);

    // Build the consensus string and check it.
    string consensusString;
    for(const auto& [base, coverage]: consensus) {
        consensusString.push_back(base.character());
    }

    // The expected consensus. The two variable A runs come out at 12 and 11.
    //
    // These are the KNOWN TRUE run lengths at this locus, so this is a check
    // against ground truth and not merely against what some other tool produces.
    // The median gives 11 and 11 here and is wrong: 10 of the 19 reads report
    // the first run as shorter than 12, so the median lands one base low even
    // though 12 is the single most common length. See the comment on
    // RunLengthEstimator.
    const string expectedConsensus =
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT"
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT";
    cout << "Consensus length " << consensusString.size() << "." << endl;
    SHASTA2_ASSERT(consensusString == expectedConsensus);

    // The consensus is the aligned consensus with the gaps removed.
    {
        string fromAligned;
        for(uint64_t j=0; j<alignedExtendedConsensus.size(); j++) {
            const ExtendedBase e = alignedExtendedConsensus[j];
            if(e.isGap()) {
                continue;
            }
            for(uint64_t k=0; k<consensusRunLengths[j]; k++) {
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
    // 11 rather than the true 12, in each of the two tandem copies. Ten of the
    // nineteen reads report that run as shorter than it is, so the median lands
    // below the truth even though 12 is the most common single length.
    //
    // Asserted so that the two estimators cannot silently converge, and so that
    // the direction of the difference is recorded: if a future change makes the
    // median agree with the mode here, something has moved that needs looking at.
    {
        vector< pair<Base, uint64_t> > medianConsensus;
        vector<ExtendedBase> ignore1;
        vector<uint64_t> ignore2;
        extendedConsensus(extendedAlignment, alignedRunLengthsAllRows, weights,
            RunLengthEstimator::Median, medianConsensus, ignore1, ignore2);
        string medianString;
        for(const auto& [base, coverage]: medianConsensus) {
            medianString.push_back(base.character());
        }
        cout << "Mode gives length " << consensusString.size() <<
            " (correct), median gives length " << medianString.size() << "." << endl;
        SHASTA2_ASSERT(medianString != consensusString);
        SHASTA2_ASSERT(medianString.size() == consensusString.size() - 2);
    }



    // Expanding the alignment must preserve every read, and give every row the
    // same length as the aligned consensus.
    {
        vector< vector<AlignedBase> > expandedAlignment;
        vector<AlignedBase> expandedAlignedConsensus;
        expandExtendedAlignment(extendedAlignment, alignedRunLengthsAllRows,
            alignedExtendedConsensus, consensusRunLengths,
            expandedAlignment, expandedAlignedConsensus);

        SHASTA2_ASSERT(expandedAlignment.size() == badAlignmentRows.size());
        for(uint64_t i=0; i<expandedAlignment.size(); i++) {

            // Same length as the aligned consensus.
            SHASTA2_ASSERT(expandedAlignment[i].size() == expandedAlignedConsensus.size());

            // Removing the gaps gives the original read back.
            string ungapped;
            for(const AlignedBase b: expandedAlignment[i]) {
                if(not b.isGap()) {
                    ungapped.push_back(b.character());
                }
            }
            string expected;
            for(const char c: badAlignmentRows[i]) {
                if(c != '-') {
                    expected.push_back(c);
                }
            }
            SHASTA2_ASSERT(ungapped == expected);
        }

        // No column of the expanded alignment holds two different bases.
        // The input had 4 such columns. This is the defect being fixed, and it is
        // structural here: a poly symbol and the base next to it are different
        // symbols and get different columns, so they cannot be confused.
        uint64_t impureColumnCount = 0;
        for(uint64_t j=0; j<expandedAlignedConsensus.size(); j++) {
            array<bool, 4> seen = {false, false, false, false};
            for(const vector<AlignedBase>& row: expandedAlignment) {
                const AlignedBase b = row[j];
                if(not b.isGap()) {
                    seen[b.value] = true;
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
        // Reads 9 and 10 have n1 = 10, n2 = 12.
        skewedWeights[9] = 100;
        skewedWeights[10] = 100;

        for(const RunLengthEstimator estimator:
            {RunLengthEstimator::Median, RunLengthEstimator::Mode}) {
            vector< pair<Base, uint64_t> > skewedConsensus;
            vector<ExtendedBase> ignore1;
            vector<uint64_t> ignore2;
            extendedConsensus(extendedAlignment, alignedRunLengthsAllRows,
                skewedWeights, estimator, skewedConsensus, ignore1, ignore2);
            string skewedString;
            for(const auto& [base, coverage]: skewedConsensus) {
                skewedString.push_back(base.character());
            }

            // Both estimators must now report the run lengths of the heavy
            // reads, 10 and 12, in each of the two tandem copies.
            const string expected =
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT" +
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT";
            SHASTA2_ASSERT(skewedString == expected);
            SHASTA2_ASSERT(skewedString != consensusString);
        }
    }



    // A real substitution must survive. Change one base in some rows to make a
    // genuine two allele column, and check the consensus takes the majority
    // allele rather than absorbing or dropping it.
    {
        vector< vector<ExtendedBase> > snpAlignment = extendedAlignment;
        vector< vector<uint64_t> > snpRunLengths = alignedRunLengthsAllRows;

        // Find a column holding a plain symbol in every row.
        uint64_t column = invalid<uint64_t>;
        for(uint64_t j=0; j<snpAlignment.front().size(); j++) {
            if(not snpAlignment.front()[j].isPoly()) {
                column = j;
                break;
            }
        }
        SHASTA2_ASSERT(column != invalid<uint64_t>);
        const ExtendedBase original = snpAlignment.front()[column];

        // Flip it in 8 of the 19 rows, so the original stays in the majority.
        const ExtendedBase mutated = ExtendedBase::fromBase(original.base().complement());
        SHASTA2_ASSERT(mutated != original);
        for(uint64_t i=0; i<8; i++) {
            snpAlignment[i][column] = mutated;
        }

        vector< pair<Base, uint64_t> > snpConsensus;
        vector<ExtendedBase> snpAlignedConsensus;
        vector<uint64_t> snpConsensusRunLengths;
        extendedConsensus(snpAlignment, snpRunLengths, weights,
            RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus, snpConsensusRunLengths);

        // The majority allele wins, and the column is not turned into a gap.
        SHASTA2_ASSERT(snpAlignedConsensus[column] == original);
        SHASTA2_ASSERT(snpConsensusRunLengths[column] == 1);

        // Its coverage reflects the split: 11 of 19, not all 19.
        uint64_t positionInConsensus = 0;
        for(uint64_t j=0; j<column; j++) {
            positionInConsensus += snpConsensusRunLengths[j];
        }
        SHASTA2_ASSERT(snpConsensus[positionInConsensus].second == 11);
        cout << "With a substitution in 8 of 19 rows, the majority allele is kept "
            "with coverage " << snpConsensus[positionInConsensus].second << "." << endl;

        // Now flip it in 12 rows, making the new allele the majority, and check
        // the consensus follows.
        for(uint64_t i=0; i<12; i++) {
            snpAlignment[i][column] = mutated;
        }
        extendedConsensus(snpAlignment, snpRunLengths, weights,
            RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus, snpConsensusRunLengths);
        SHASTA2_ASSERT(snpAlignedConsensus[column] == mutated);
    }



    // A column that is a gap in a majority of rows must become a gap, and
    // contribute nothing to the consensus.
    {
        vector< vector<ExtendedBase> > gapAlignment = extendedAlignment;
        vector< vector<uint64_t> > gapRunLengths = alignedRunLengthsAllRows;
        for(uint64_t i=0; i<gapAlignment.size(); i++) {
            gapAlignment[i].insert(gapAlignment[i].begin(), ExtendedBase::gap());
            gapRunLengths[i].insert(gapRunLengths[i].begin(), 0);
        }
        // One row has a real base in the new column, the other 18 have gaps.
        gapAlignment[0][0] = ExtendedBase::fromCharacter('C');
        gapRunLengths[0][0] = 1;

        vector< pair<Base, uint64_t> > gapConsensus;
        vector<ExtendedBase> gapAlignedConsensus;
        vector<uint64_t> gapConsensusRunLengths;
        extendedConsensus(gapAlignment, gapRunLengths, weights,
            RunLengthEstimator::Mode,
            gapConsensus, gapAlignedConsensus, gapConsensusRunLengths);
        SHASTA2_ASSERT(gapAlignedConsensus[0].isGap());
        SHASTA2_ASSERT(gapConsensusRunLengths[0] == 0);

        // The consensus is unchanged by the added column.
        string gapString;
        for(const auto& [base, coverage]: gapConsensus) {
            gapString.push_back(base.character());
        }
        SHASTA2_ASSERT(gapString == consensusString);
    }



    // No regression when nothing is long enough to become a poly symbol.
    // With a threshold above every run, the encoding is the identity, so this
    // must reduce to plain weighted column majority.
    {
        const vector<string> rows = {
            "ACGTACGTAC",
            "ACGTACGTAC",
            "ACGTTCGTAC"
        };
        const vector<uint64_t> rowWeights = {1, 1, 1};
        vector< vector<ExtendedBase> > a;
        vector< vector<uint64_t> > r;
        for(const string& row: rows) {
            a.push_back(vectorOfExtendedBasesFromString(row));
            r.push_back(vector<uint64_t>(row.size(), 1));
        }
        for(const vector<ExtendedBase>& row: a) {
            for(const ExtendedBase e: row) {
                SHASTA2_ASSERT(not e.isPoly());
            }
        }
        vector< pair<Base, uint64_t> > c;
        vector<ExtendedBase> ac;
        vector<uint64_t> crl;
        extendedConsensus(a, r, rowWeights, RunLengthEstimator::Mode, c, ac, crl);
        string s;
        for(const auto& [base, coverage]: c) {
            s.push_back(base.character());
        }
        // Plain column majority of the three rows.
        SHASTA2_ASSERT(s == "ACGTACGTAC");
        // The column with the substitution has coverage 2, the rest 3.
        SHASTA2_ASSERT(c[4].second == 2);
        SHASTA2_ASSERT(c[0].second == 3);
    }



    // The encoding must survive a round trip through a string, because that is
    // how it reaches theseus and comes back: theseusExtended passes
    // toString(encoded) in and parses the MSA rows with
    // vectorOfExtendedBasesFromString. If ExtendedBase::fromCharacter ever
    // folded 'a' into 'A', as AlignedBase::fromCharacter does, every poly symbol
    // would silently become a plain one and the whole scheme would collapse
    // without any visible error.
    {
        for(const vector<ExtendedBase>& row: extendedAlignment) {
            const string s = toString(row);
            SHASTA2_ASSERT(vectorOfExtendedBasesFromString(s) == row);
        }

        // Including gaps, which theseus emits in its MSA rows.
        const string withGaps = "CTC-a-G-a--GTT";
        const vector<ExtendedBase> parsed = vectorOfExtendedBasesFromString(withGaps);
        SHASTA2_ASSERT(toString(parsed) == withGaps);

        // And the case distinction really is preserved.
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a') != ExtendedBase::fromCharacter('A'));
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').isPoly());
        SHASTA2_ASSERT(not ExtendedBase::fromCharacter('A').isPoly());
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').base() ==
            ExtendedBase::fromCharacter('A').base());

        // AlignedBase, by contrast, would lose it. This documents exactly why
        // theseusExtended must not use vectorOfAlignedBasesFromString.
        SHASTA2_ASSERT(AlignedBase::fromCharacter('a') == AlignedBase::fromCharacter('A'));
    }



    // alignedRunLengths lines run lengths up with the columns of a row.
    {
        const vector<ExtendedBase> encoded = vectorOfExtendedBasesFromString("CTCaGa");
        const vector<uint64_t> runLengths = {1, 1, 1, 12, 1, 11};
        const vector<ExtendedBase> row = vectorOfExtendedBasesFromString("C-TC-aG-a");
        vector<uint64_t> aligned;
        alignedRunLengths(row, encoded, runLengths, Anchoring::BothSides, aligned);
        const vector<uint64_t> expected = {1, 0, 1, 1, 0, 12, 1, 0, 11};
        SHASTA2_ASSERT(aligned == expected);
    }



    // A row can cover only part of its encoding.
    //
    // Theseus trims the overhang of a sequence fixed on one side only: a left
    // fixed sequence loses symbols off its right end, a right fixed one off its
    // left. So the aligned row is a contiguous window of the encoding and the run
    // lengths have to be read from the matching offset. This was found by running
    // msa1 on theseus's own test data, where a left fixed sequence ending in a
    // poly-T came back with the poly-T removed. Reading the run lengths from
    // offset 0 regardless would pair each symbol with another symbol's length.
    {
        // A left fixed sequence: the trailing poly symbol is trimmed.
        const vector<ExtendedBase> encoded = vectorOfExtendedBasesFromString("CTCaGat");
        const vector<uint64_t> runLengths = {1, 1, 1, 12, 1, 11, 7};
        const vector<ExtendedBase> row = vectorOfExtendedBasesFromString("CTC-aGa");
        vector<uint64_t> aligned;
        alignedRunLengths(row, encoded, runLengths, Anchoring::LeftOnly, aligned);
        const vector<uint64_t> expected = {1, 1, 1, 0, 12, 1, 11};
        SHASTA2_ASSERT(aligned == expected);
    }
    {
        // A right fixed sequence: the leading poly symbol is trimmed, so the run
        // lengths must be read from offset 1, not 0.
        const vector<ExtendedBase> encoded = vectorOfExtendedBasesFromString("gCTCaGa");
        const vector<uint64_t> runLengths = {9, 1, 1, 1, 12, 1, 11};
        const vector<ExtendedBase> row = vectorOfExtendedBasesFromString("CTCaG-a");
        vector<uint64_t> aligned;
        alignedRunLengths(row, encoded, runLengths, Anchoring::RightOnly, aligned);
        const vector<uint64_t> expected = {1, 1, 1, 12, 1, 0, 11};
        SHASTA2_ASSERT(aligned == expected);
        // The trimmed poly symbol's length of 9 must not appear anywhere.
        for(const uint64_t length: aligned) {
            SHASTA2_ASSERT(length != 9);
        }
    }



    // The window must come from the anchoring, not from searching the encoding
    // for the row's symbols.
    //
    // Here the encoding is a tandem repeat, A a C A a C, whose two poly runs
    // have different lengths, 10 and 20. The sequence is fixed on the right, so
    // theseus trimmed the first copy and the row is the second one, whose run is
    // 20. Searching for the row's symbols finds the FIRST copy and would pair the
    // row with 10 instead. Tandem repeats are exactly what this code is for, and
    // the 19 read locus it was developed on is itself a 2x tandem duplication, so
    // this is not a remote case.
    {
        const vector<ExtendedBase> encoded = vectorOfExtendedBasesFromString("AaCAaC");
        const vector<uint64_t> runLengths = {1, 10, 1, 1, 20, 1};
        const vector<ExtendedBase> row = vectorOfExtendedBasesFromString("AaC");

        vector<uint64_t> aligned;
        alignedRunLengths(row, encoded, runLengths, Anchoring::RightOnly, aligned);
        const vector<uint64_t> expected = {1, 20, 1};
        SHASTA2_ASSERT(aligned == expected);

        // And the mirror image: fixed on the left, so the row is the first copy
        // and its run is 10.
        alignedRunLengths(row, encoded, runLengths, Anchoring::LeftOnly, aligned);
        const vector<uint64_t> expectedLeft = {1, 10, 1};
        SHASTA2_ASSERT(aligned == expectedLeft);
    }



    // Degenerate inputs.
    {
        vector< pair<Base, uint64_t> > c;
        vector<ExtendedBase> ac;
        vector<uint64_t> crl;

        // Empty alignment.
        extendedConsensus({}, {}, {}, RunLengthEstimator::Mode, c, ac, crl);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.empty());

        // A single sequence: the consensus is that sequence.
        const vector< vector<ExtendedBase> > one = {vectorOfExtendedBasesFromString("CTCaGa")};
        const vector< vector<uint64_t> > oneRunLengths = {{1, 1, 1, 12, 1, 11}};
        extendedConsensus(one, oneRunLengths, {7}, RunLengthEstimator::Mode, c, ac, crl);
        string s;
        for(const auto& [base, coverage]: c) {
            s.push_back(base.character());
            SHASTA2_ASSERT(coverage == 7);
        }
        SHASTA2_ASSERT(s == "CTC" + string(12, 'A') + "G" + string(11, 'A'));

        // Every row a gap at every column.
        const vector< vector<ExtendedBase> > allGaps = {
            vectorOfExtendedBasesFromString("--"),
            vectorOfExtendedBasesFromString("--")};
        const vector< vector<uint64_t> > allGapRunLengths = {{0, 0}, {0, 0}};
        extendedConsensus(allGaps, allGapRunLengths, {1, 1},
            RunLengthEstimator::Mode, c, ac, crl);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.size() == 2);
        SHASTA2_ASSERT(ac[0].isGap() and ac[1].isGap());
    }

    cout << "testMsa1Consensus passed." << endl;
}
