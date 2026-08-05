// Shasta2.
#include "msa1.hpp"
#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"
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



// See msa1.hpp for comments.
bool shasta2::msa1PatternPresent(
    const vector<Base>& sequence,
    uint64_t threshold)
{
    // Walk the runs, keeping only the two previous run lengths. The pattern is a
    // run longer than the threshold, then a run of exactly one base, then
    // another run longer than the threshold. Nothing is allocated and the walk
    // stops at the first hit, which is what lets this be called on every local
    // assembly before deciding whether an alignment is even needed.
    uint64_t lengthBeforeLast = 0;
    uint64_t lastLength = 0;

    const uint64_t n = sequence.size();
    for(uint64_t runBegin=0; runBegin<n; /* incremented below */) {
        const Base base = sequence[runBegin];
        uint64_t runEnd = runBegin + 1;
        while((runEnd < n) and (sequence[runEnd] == base)) {
            ++runEnd;
        }
        const uint64_t runLength = runEnd - runBegin;

        if( (lengthBeforeLast > threshold) and
            (lastLength == 1) and
            (runLength > threshold)) {
            return true;
        }

        lengthBeforeLast = lastLength;
        lastLength = runLength;
        runBegin = runEnd;
    }
    return false;
}



// See msa1.hpp for comments.
bool shasta2::msa1PatternPresent(
    const vector< vector<Base> >& sequences,
    uint64_t threshold)
{
    for(const vector<Base>& sequence: sequences) {
        if(msa1PatternPresent(sequence, threshold)) {
            return true;
        }
    }
    return false;
}



// See msa1.hpp for comments.
bool shasta2::msa1PatternPresent(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    uint64_t threshold)
{
    for(const auto& [sequence, weight]: sequences) {
        if(msa1PatternPresent(sequence, threshold)) {
            return true;
        }
    }
    return false;
}



namespace shasta2 {

    // Return true if a region boundary at this column falls inside a homopolymer
    // run of any row, that is, if the last base of some row before the boundary
    // is the same as its next base at or after it.
    //
    // A boundary that cuts a run is a problem for everything downstream: the
    // pattern test would see a run shorter than it is and might not recognise
    // the pattern at all, and the encoding would collapse a partial run.
    static bool msa1BoundaryCutsRun(
        const vector< vector<AlignedBase> >& alignment,
        uint64_t column)
    {
        for(const vector<AlignedBase>& row: alignment) {

            // The last base before the boundary.
            AlignedBase before = AlignedBase::gap();
            for(uint64_t k=column; k>0; k--) {
                if(not row[k - 1].isGap()) {
                    before = row[k - 1];
                    break;
                }
            }
            if(before.isGap()) {
                continue;
            }

            // The first base at or after the boundary.
            AlignedBase after = AlignedBase::gap();
            for(uint64_t k=column; k<row.size(); k++) {
                if(not row[k].isGap()) {
                    after = row[k];
                    break;
                }
            }
            if(after.isGap()) {
                continue;
            }

            if(before == after) {
                return true;
            }
        }
        return false;
    }
}



// See msa1.hpp for comments.
void shasta2::msa1FindBadRegions(
    const vector< vector<AlignedBase> >& alignment,
    const vector<AlignedBase>& alignedConsensus,
    uint64_t threshold,
    uint64_t flank,
    uint64_t mergeDistance,
    vector<Msa1Region>& regions)
{
    regions.clear();
    if(alignment.empty()) {
        return;
    }
    const uint64_t alignmentLength = alignedConsensus.size();
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignmentLength);
    }

    // A column is discordant if any row differs from the aligned consensus
    // there. This is what LocalAssembly7::writeAlignment paints pink.
    vector<bool> isDiscordant(alignmentLength, false);
    for(uint64_t j=0; j<alignmentLength; j++) {
        for(const vector<AlignedBase>& row: alignment) {
            if(row[j] != alignedConsensus[j]) {
                isDiscordant[j] = true;
                break;
            }
        }
    }

    // Cluster the discordant columns, merging clusters that are close together,
    // then widen each cluster by the flank.
    for(uint64_t j=0; j<alignmentLength; /* incremented below */) {
        if(not isDiscordant[j]) {
            ++j;
            continue;
        }

        uint64_t begin = j;
        uint64_t end = j + 1;
        for(uint64_t k=end; k<alignmentLength; k++) {
            if(isDiscordant[k]) {
                end = k + 1;
            } else if(k >= end + mergeDistance) {
                break;
            }
        }
        j = end;

        // Widen by the flank, clamped to the alignment.
        Msa1Region region;
        region.begin = (begin > flank) ? (begin - flank) : 0;
        region.end = min(alignmentLength, end + flank);

        // Then widen further, until neither boundary falls inside a homopolymer
        // run of any row.
        //
        // Without this the flank has to be guessed large enough to clear the
        // runs, and a window that clips one silently does nothing: the pattern
        // test sees a run shorter than it is and does not recognise the pattern.
        // On the alignment this was developed against, a flank of 4 found the
        // regions but repaired neither, and a flank of 0 found none at all,
        // while 6 and above worked. That is not a property to leave to a
        // constant, since a locus with longer runs would need a larger one.
        // Snapping the boundaries to run boundaries removes the dependence: the
        // flank is now only a minimum amount of context.
        while((region.begin > 0) and msa1BoundaryCutsRun(alignment, region.begin)) {
            --region.begin;
        }
        while((region.end < alignmentLength) and msa1BoundaryCutsRun(alignment, region.end)) {
            ++region.end;
        }

        // Merge with the previous region if they now touch or overlap.
        if((not regions.empty()) and (region.begin <= regions.back().end)) {
            regions.back().end = max(regions.back().end, region.end);
        } else {
            regions.push_back(region);
        }
    }

    // Keep only the regions whose sequences actually show the problem pattern.
    // This is the confirmation step: a discordant cluster is common, the pattern
    // is not.
    vector<Msa1Region> confirmed;
    vector<Base> windowSequence;
    for(const Msa1Region& region: regions) {
        bool found = false;
        for(const vector<AlignedBase>& row: alignment) {
            windowSequence.clear();
            for(uint64_t j=region.begin; j<region.end; j++) {
                if(not row[j].isGap()) {
                    windowSequence.push_back(Base(row[j]));
                }
            }
            if(msa1PatternPresent(windowSequence, threshold)) {
                found = true;
                break;
            }
        }
        if(found) {
            confirmed.push_back(region);
        }
    }
    regions = confirmed;
}



namespace shasta2 {

    // Recompute one region of an alignment using the extended alphabet, and
    // return the replacement columns.
    //
    // The rows of the region are ungapped and encoded. They are repaired only if
    // their symbol strings all agree, in which case they stack without any
    // alignment and the usual two votes run, symbols first and run lengths
    // second, before expanding back to bases.
    //
    // Returns false and leaves the outputs alone if the region cannot be
    // improved, in which case the caller must leave it as it was.
    static bool msa1RepairRegion(
        const vector< vector<AlignedBase> >& alignment,
        const Msa1Region& region,
        const vector<uint64_t>& weights,
        uint64_t threshold,
        RunLengthEstimator estimator,
        vector< vector<AlignedBase> >& newRows,
        vector<AlignedBase>& newAlignedConsensus,
        vector< pair<Base, uint64_t> >& newConsensus)
    {
        const uint64_t n = alignment.size();

        // Ungap each row inside the region and encode it.
        vector< pair<ExtendedSequence, uint64_t> > encoded;
        vector<uint64_t> rowWeights;
        vector<uint64_t> nonEmptyRows;
        for(uint64_t i=0; i<n; i++) {
            vector<Base> windowSequence;
            for(uint64_t j=region.begin; j<region.end; j++) {
                if(not alignment[i][j].isGap()) {
                    windowSequence.push_back(Base(alignment[i][j]));
                }
            }
            ExtendedSequence e;
            encodeExtended(windowSequence, threshold, e);
            encoded.push_back(make_pair(e, weights[i]));
            rowWeights.push_back(weights[i]);
            if(not windowSequence.empty()) {
                nonEmptyRows.push_back(i);
            }
        }

        // A region in which some row contributes nothing is left alone. Such a
        // row cannot be given a symbol string to align, and guessing one would
        // be a change to a row that has no evidence in this region.
        if(nonEmptyRows.size() != n) {
            return false;
        }

        // All rows must encode to the same symbol string.
        //
        // When they do, the homopolymer lengths are the only thing that differs
        // between the reads in this window, the rows stack without any
        // alignment at all, and the repair is exact: there is no aligner left to
        // make a placement decision, so the defect being repaired cannot recur.
        // This is the case the extended alphabet is designed to produce and it is
        // what the locus this was developed on does.
        //
        // When they do not, the window has real variation in it as well as the
        // homopolymers, and repairing it would need the rows realigned in symbol
        // space. That is deliberately not done. Theseus is the only aligner here
        // that could do it, since abpoa is fixed at four letters, and theseus is
        // the worse aligner of the two; replacing abpoa's alignment of a window
        // with a worse one to fix a homopolymer would trade one error for
        // another. Such a window is left exactly as abpoa produced it.
        //
        // The cost of that choice is that a bad homopolymer sitting next to a
        // substitution is not repaired. Widening this is possible later by
        // clipping every run at the threshold and realigning with abpoa, which
        // gets the same effect as poly symbols in an alphabet abpoa accepts, but
        // it needs measuring against real data first.
        for(uint64_t i=1; i<n; i++) {
            if(toString(encoded[i].first) != toString(encoded[0].first)) {
                return false;
            }
        }

        vector<AlignedExtendedSequence> extendedAlignment(n);
        for(uint64_t i=0; i<n; i++) {
            extendedAlignment[i].clear();
            extendedAlignment[i].reserve(encoded[i].first.size());
            for(const auto& [symbol, runLength]: encoded[i].first) {
                extendedAlignment[i].push_back(
                    make_pair(AlignedExtendedBase(symbol), runLength));
            }
        }

        // Vote, then expand.
        AlignedExtendedSequence alignedExtendedConsensus;
        extendedConsensus(extendedAlignment, rowWeights, estimator,
            newConsensus, alignedExtendedConsensus);
        expandExtendedAlignment(extendedAlignment, alignedExtendedConsensus,
            newRows, newAlignedConsensus);

        // The bases of the region consensus must be the non-gap entries of the
        // aligned region consensus, in order. The caller splices both, so if
        // they disagreed the two would drift apart.
        uint64_t k = 0;
        for(const AlignedBase b: newAlignedConsensus) {
            if(not b.isGap()) {
                SHASTA2_ASSERT(k < newConsensus.size());
                SHASTA2_ASSERT(newConsensus[k].first == Base(b));
                ++k;
            }
        }
        SHASTA2_ASSERT(k == newConsensus.size());

        return true;
    }
}



// See msa1.hpp for comments.
uint64_t shasta2::msa1(
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus,
    vector< pair<Base, uint64_t> >& consensus,
    const vector<uint64_t>& weights,
    uint64_t threshold,
    RunLengthEstimator estimator,
    uint64_t flank,
    uint64_t mergeDistance)
{
    const uint64_t n = alignment.size();
    if(n == 0) {
        return 0;
    }
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignedConsensus.size());
    }

    // An empty weights argument means every row weighs 1.
    vector<uint64_t> rowWeights = weights;
    if(rowWeights.empty()) {
        rowWeights.assign(n, 1);
    }
    SHASTA2_ASSERT(rowWeights.size() == n);

    // Find the regions worth repairing. Usually there are none, and then nothing
    // below runs and nothing is modified.
    vector<Msa1Region> regions;
    msa1FindBadRegions(alignment, alignedConsensus, threshold, flank,
        mergeDistance, regions);
    if(regions.empty()) {
        return 0;
    }

    // The consensus must agree with the aligned consensus before anything is
    // spliced, because the base positions of a region are found by counting
    // non-gap entries of the aligned consensus.
    {
        uint64_t k = 0;
        for(const AlignedBase b: alignedConsensus) {
            if(not b.isGap()) {
                ++k;
            }
        }
        SHASTA2_ASSERT(k == consensus.size());
    }

    // Repair the regions from right to left, so that the column indices of the
    // regions still to be done, and the base positions of the consensus before
    // them, are not disturbed by the splices already made.
    uint64_t repairedCount = 0;
    for(uint64_t k=regions.size(); k>0; k--) {
        const Msa1Region& region = regions[k - 1];

        vector< vector<AlignedBase> > newRows;
        vector<AlignedBase> newAlignedConsensus;
        vector< pair<Base, uint64_t> > newConsensus;
        if(not msa1RepairRegion(alignment, region, rowWeights, threshold,
            estimator, newRows, newAlignedConsensus, newConsensus)) {
            continue;
        }
        SHASTA2_ASSERT(newRows.size() == n);

        // Find the range of consensus BASES this region covers, by counting the
        // non-gap entries of the aligned consensus. This has to be done before
        // the splice below, which changes the aligned consensus.
        uint64_t baseBegin = 0;
        for(uint64_t j=0; j<region.begin; j++) {
            if(not alignedConsensus[j].isGap()) {
                ++baseBegin;
            }
        }
        uint64_t baseEnd = baseBegin;
        for(uint64_t j=region.begin; j<region.end; j++) {
            if(not alignedConsensus[j].isGap()) {
                ++baseEnd;
            }
        }
        SHASTA2_ASSERT(baseEnd <= consensus.size());

        // Splice the alignment rows and the aligned consensus, in column space.
        // These must all be spliced together or they stop having equal length.
        for(uint64_t i=0; i<n; i++) {
            vector<AlignedBase>& row = alignment[i];
            row.erase(row.begin() + int64_t(region.begin),
                row.begin() + int64_t(region.end));
            row.insert(row.begin() + int64_t(region.begin),
                newRows[i].begin(), newRows[i].end());
        }
        alignedConsensus.erase(
            alignedConsensus.begin() + int64_t(region.begin),
            alignedConsensus.begin() + int64_t(region.end));
        alignedConsensus.insert(
            alignedConsensus.begin() + int64_t(region.begin),
            newAlignedConsensus.begin(), newAlignedConsensus.end());

        // Splice the consensus, in base space. The same erase and insert, just
        // indexed by base rather than by column.
        //
        // The consensus is spliced rather than recomputed from the repaired
        // alignment. Recomputing it would rewrite the coverage of every position
        // in the sequence, including all the positions outside any repaired
        // region, because the coverage abpoa reports is computed its own way and
        // is not necessarily what a recount of the rows would give. Rewriting
        // them would silently discard abpoa's numbers everywhere, which is
        // exactly what repairing locally is supposed to avoid.
        consensus.erase(
            consensus.begin() + int64_t(baseBegin),
            consensus.begin() + int64_t(baseEnd));
        consensus.insert(
            consensus.begin() + int64_t(baseBegin),
            newConsensus.begin(), newConsensus.end());

        ++repairedCount;
    }

    if(repairedCount == 0) {
        return 0;
    }

    // All rows, the aligned consensus and the consensus must still agree.
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignedConsensus.size());
    }
    {
        uint64_t k = 0;
        for(const AlignedBase b: alignedConsensus) {
            if(not b.isGap()) {
                SHASTA2_ASSERT(k < consensus.size());
                SHASTA2_ASSERT(consensus[k].first == Base(b));
                ++k;
            }
        }
        SHASTA2_ASSERT(k == consensus.size());
    }

    return repairedCount;
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

    static string msa1ToString(const vector<AlignedBase>& v)
    {
        string s;
        for(const AlignedBase b: v) {
            s.push_back(b.character());
        }
        return s;
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


// Test the local repair.
//
// The property that matters most here is the first one: on input with no bad
// region, nothing is modified at all. The whole point of repairing locally is
// that sequence accuracy in the rest of the assembly is left alone, so a test
// that the repair is a no-op where it should be is more important than a test
// that it works where it should.
void shasta2::testMsa1Repair()
{
    const uint64_t threshold = defaultHomopolymerThreshold;

    // Build the real, defective abpoa alignment.
    vector< vector<AlignedBase> > alignment;
    vector<string> reads;
    for(const string& row: msa1BadAlignmentRows) {
        alignment.push_back(vectorOfAlignedBasesFromString(row));
        string ungapped;
        for(const char c: row) {
            if(c != '-') {
                ungapped.push_back(c);
            }
        }
        reads.push_back(ungapped);
    }
    const vector<uint64_t> weights(alignment.size(), 1);

    // Its column-wise consensus, which is what abpoa reports.
    const auto columnConsensus = [&](
        const vector< vector<AlignedBase> >& a,
        vector<AlignedBase>& ac,
        vector< pair<Base, uint64_t> >& c)
    {
        ac.clear();
        c.clear();
        for(uint64_t j=0; j<a.front().size(); j++) {
            array<uint64_t, 5> w;
            fill(w.begin(), w.end(), 0);
            for(uint64_t i=0; i<a.size(); i++) {
                w[a[i][j].value] += weights[i];
            }
            const auto it = std::ranges::max_element(w);
            const AlignedBase b = AlignedBase::fromInteger(uint64_t(it - w.begin()));
            ac.push_back(b);
            if(not b.isGap()) {
                c.push_back(make_pair(Base(b), *it));
            }
        }
    };
    vector<AlignedBase> alignedConsensus;
    vector< pair<Base, uint64_t> > consensus;
    columnConsensus(alignment, alignedConsensus, consensus);



    // The prescreen must fire on these reads, which contain two long A runs
    // separated by a single G.
    {
        vector< vector<Base> > sequences;
        for(const string& r: reads) {
            sequences.push_back(vectorOfBasesFromString(r));
        }
        SHASTA2_ASSERT(msa1PatternPresent(sequences, threshold));

        // And it must not fire on ordinary sequence with no long runs. This is
        // the common case in real reads, where 97% of runs are 4 bases or
        // shorter, and it is what keeps the repair from running at all.
        const vector< vector<Base> > ordinary = {
            vectorOfBasesFromString("ACGTACGTAACCGGTTACGTACGT"),
            vectorOfBasesFromString("ACGTACGTAACCGGTTACGTACGT")};
        SHASTA2_ASSERT(not msa1PatternPresent(ordinary, threshold));

        // Nor on a single long run with no separating base.
        const vector< vector<Base> > oneRun = {
            vectorOfBasesFromString("ACGT" + string(20, 'A') + "CGTA")};
        SHASTA2_ASSERT(not msa1PatternPresent(oneRun, threshold));

        // Nor when the two long runs are separated by more than one base.
        const vector< vector<Base> > twoBases = {
            vectorOfBasesFromString(string(12, 'A') + "GT" + string(11, 'A'))};
        SHASTA2_ASSERT(not msa1PatternPresent(twoBases, threshold));

        // But yes when they are separated by exactly one.
        const vector< vector<Base> > oneBase = {
            vectorOfBasesFromString(string(12, 'A') + "G" + string(11, 'A'))};
        SHASTA2_ASSERT(msa1PatternPresent(oneBase, threshold));

        // The exact boundaries of the pattern, on the single sequence primitive.
        // A run must be strictly longer than the threshold on BOTH sides, and
        // the separator must be exactly one base.
        const auto pattern = [&](uint64_t left, uint64_t middle, uint64_t right) {
            return msa1PatternPresent(vectorOfBasesFromString(
                string(left, 'A') + string(middle, 'G') + string(right, 'A')),
                threshold);
        };
        SHASTA2_ASSERT(pattern(threshold + 1, 1, threshold + 1));
        SHASTA2_ASSERT(not pattern(threshold, 1, threshold + 1));  // left too short
        SHASTA2_ASSERT(not pattern(threshold + 1, 1, threshold));  // right too short
        SHASTA2_ASSERT(not pattern(threshold + 1, 2, threshold + 1)); // separator too long
        SHASTA2_ASSERT(pattern(100, 1, 100));

        // The pattern is found wherever it sits, including at the very start and
        // the very end, which a scan that keeps only two previous run lengths
        // could get wrong.
        SHASTA2_ASSERT(msa1PatternPresent(vectorOfBasesFromString(
            string(12, 'A') + "G" + string(11, 'A') + "CGTACGT"), threshold));
        SHASTA2_ASSERT(msa1PatternPresent(vectorOfBasesFromString(
            "CGTACGT" + string(12, 'A') + "G" + string(11, 'A')), threshold));

        // Degenerate inputs.
        SHASTA2_ASSERT(not msa1PatternPresent(vector<Base>(), threshold));
        SHASTA2_ASSERT(not msa1PatternPresent(
            vectorOfBasesFromString("A"), threshold));
        SHASTA2_ASSERT(not msa1PatternPresent(
            vector<Base>(1000, Base::fromCharacter('A')), threshold));
    }



    // THE IMPORTANT TEST. On an alignment with no bad region, nothing changes.
    {
        const vector<string> cleanRows = {
            "ACGTACGTAACCGGTTACGTACGT",
            "ACGTACGTAACCGGTTACGTACGT",
            "ACGTACGTAACCTGTTACGTACGT"};
        vector< vector<AlignedBase> > a;
        for(const string& r: cleanRows) {
            a.push_back(vectorOfAlignedBasesFromString(r));
        }
        vector<AlignedBase> ac;
        vector< pair<Base, uint64_t> > c;
        const vector<uint64_t> w(a.size(), 1);
        {
            // Local column consensus for this small case.
            for(uint64_t j=0; j<a.front().size(); j++) {
                array<uint64_t, 5> t;
                fill(t.begin(), t.end(), 0);
                for(uint64_t i=0; i<a.size(); i++) {
                    t[a[i][j].value] += w[i];
                }
                const auto it = std::ranges::max_element(t);
                const AlignedBase b = AlignedBase::fromInteger(uint64_t(it - t.begin()));
                ac.push_back(b);
                if(not b.isGap()) {
                    c.push_back(make_pair(Base(b), *it));
                }
            }
        }

        // Note this alignment DOES have a discordant column, at the T for C
        // substitution, so the region finder has something to cluster. It is the
        // pattern confirmation that must reject it.
        const auto aBefore = a;
        const auto acBefore = ac;
        const auto cBefore = c;

        const uint64_t repaired = msa1(a, ac, c, w, threshold);
        SHASTA2_ASSERT(repaired == 0);
        SHASTA2_ASSERT(a == aBefore);
        SHASTA2_ASSERT(ac == acBefore);
        SHASTA2_ASSERT(c == cBefore);
        cout << "Repair left a clean alignment untouched." << endl;
    }



    // On the real defective alignment, the repair fires and fixes it.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;

        // Before: 4 columns hold more than one base.
        {
            vector<string> rows;
            for(const auto& row: a) {
                rows.push_back(msa1ToString(row));
            }
            SHASTA2_ASSERT(msa1ImpureColumnCount(rows) == 4);
        }

        vector<Msa1Region> regions;
        msa1FindBadRegions(a, ac, threshold, 10, 20, regions);
        cout << "Found " << regions.size() << " bad region(s) in the real alignment." << endl;
        SHASTA2_ASSERT(not regions.empty());

        const uint64_t repaired = msa1(a, ac, c, weights, threshold);
        cout << "Repaired " << repaired << " region(s)." << endl;
        SHASTA2_ASSERT(repaired > 0);

        // After: no column holds more than one base.
        {
            vector<string> rows;
            for(const auto& row: a) {
                rows.push_back(msa1ToString(row));
            }
            const uint64_t impure = msa1ImpureColumnCount(rows);
            cout << "After repair the alignment has " << impure <<
                " columns containing more than one base." << endl;
            SHASTA2_ASSERT(impure == 0);
        }

        // Every row still ungaps to its original read: the repair rearranges
        // columns, it does not alter the reads.
        for(uint64_t i=0; i<a.size(); i++) {
            string ungapped;
            for(const AlignedBase b: a[i]) {
                if(not b.isGap()) {
                    ungapped.push_back(b.character());
                }
            }
            SHASTA2_ASSERT(ungapped == reads[i]);
        }

        // All rows still have the length of the aligned consensus.
        for(const auto& row: a) {
            SHASTA2_ASSERT(row.size() == ac.size());
        }

        // consensus is alignedConsensus with the gaps removed.
        {
            string fromAligned;
            for(const AlignedBase b: ac) {
                if(not b.isGap()) {
                    fromAligned.push_back(b.character());
                }
            }
            SHASTA2_ASSERT(fromAligned == msa1ToString(c));
        }

        // And the consensus is the known true one: the two A runs are 12 and 11.
        const string s = msa1ToString(c);
        const string expected =
            "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
            "AAACTATAAAGTAAATTCCTCCCATAGTT"
            "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
            "AAACTATAAAGTAAATTCCTCCCATAGTT";
        cout << "Consensus after repair has length " << s.size() << "." << endl;
        SHASTA2_ASSERT(s == expected);
    }



    // Only the bad regions change. Everything outside them is untouched, which
    // is the whole reason for repairing locally.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;

        vector<Msa1Region> regions;
        msa1FindBadRegions(a, ac, threshold, 10, 20, regions);
        SHASTA2_ASSERT(not regions.empty());

        // Record the columns before the first region and after the last one.
        const uint64_t firstBegin = regions.front().begin;
        const uint64_t lastEnd = regions.back().end;
        const uint64_t tailLength = ac.size() - lastEnd;
        vector<string> prefixBefore;
        vector<string> suffixBefore;
        for(const auto& row: a) {
            prefixBefore.push_back(msa1ToString(row).substr(0, firstBegin));
            suffixBefore.push_back(msa1ToString(row).substr(lastEnd));
        }

        msa1(a, ac, c, weights, threshold);

        // The prefix is at the same columns and unchanged.
        for(uint64_t i=0; i<a.size(); i++) {
            SHASTA2_ASSERT(msa1ToString(a[i]).substr(0, firstBegin) == prefixBefore[i]);
        }
        // The suffix is unchanged, though it has moved if the repair changed the
        // number of columns.
        for(uint64_t i=0; i<a.size(); i++) {
            const string row = msa1ToString(a[i]);
            SHASTA2_ASSERT(row.substr(row.size() - tailLength) == suffixBefore[i]);
        }
        cout << "Columns outside the repaired regions are unchanged." << endl;
    }



    // The consensus outside a repaired region is untouched, coverage included.
    //
    // This is the reason the consensus is spliced rather than recomputed from
    // the repaired alignment. Recomputing it rewrote the coverage of every
    // position in the sequence, because the coverage abpoa reports is computed
    // its own way and is not what a recount of the rows gives. That went
    // unnoticed at first because the test built the consensus with the same
    // formula the recomputation used, so the two agreed by construction. Here
    // the coverage is deliberately set to a value no recount could produce.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;
        const uint64_t sentinel = 777;
        for(auto& [base, coverage]: c) {
            coverage = sentinel;
        }
        const auto before = c;

        // Which consensus bases lie inside a region, worked out before the
        // repair moves anything.
        vector<Msa1Region> regions;
        msa1FindBadRegions(a, ac, threshold, 10, 20, regions);
        SHASTA2_ASSERT(not regions.empty());
        vector<bool> isInsideRegion(c.size(), false);
        {
            uint64_t base = 0;
            for(uint64_t j=0; j<ac.size(); j++) {
                if(ac[j].isGap()) {
                    continue;
                }
                for(const Msa1Region& region: regions) {
                    if((j >= region.begin) and (j < region.end)) {
                        isInsideRegion[base] = true;
                    }
                }
                ++base;
            }
        }
        uint64_t insideCount = 0;
        for(const bool b: isInsideRegion) {
            if(b) {
                ++insideCount;
            }
        }
        SHASTA2_ASSERT(insideCount > 0);
        SHASTA2_ASSERT(insideCount < c.size());

        SHASTA2_ASSERT(msa1(a, ac, c, weights, threshold) > 0);

        // The bases before the first repaired region and after the last must be
        // unchanged, coverage and all.
        uint64_t firstInside = 0;
        while((firstInside < isInsideRegion.size()) and (not isInsideRegion[firstInside])) {
            ++firstInside;
        }
        uint64_t lastInside = 0;
        for(uint64_t i=0; i<isInsideRegion.size(); i++) {
            if(isInsideRegion[i]) {
                lastInside = i;
            }
        }
        const uint64_t tailLength = before.size() - lastInside - 1;

        for(uint64_t i=0; i<firstInside; i++) {
            SHASTA2_ASSERT(c[i] == before[i]);
        }
        for(uint64_t i=0; i<tailLength; i++) {
            SHASTA2_ASSERT(c[c.size() - 1 - i] == before[before.size() - 1 - i]);
        }

        // And the coverage inside the regions really was replaced, or the test
        // above would pass for the wrong reason.
        uint64_t stillSentinel = 0;
        for(const auto& [base, coverage]: c) {
            if(coverage == sentinel) {
                ++stillSentinel;
            }
        }
        SHASTA2_ASSERT(stillSentinel == c.size() - insideCount);
        cout << "Consensus outside the repaired regions is unchanged: " <<
            stillSentinel << " of " << c.size() <<
            " positions keep their original coverage." << endl;
    }



    // A region where some row contributes no bases is left alone rather than
    // guessed at.
    {
        const vector<string> rows = {
            "AAAAAAAAAAAAGAAAAAAAAAAA",
            "AAAAAAAAAAAAGAAAAAAAAAAA",
            "------------------------"};
        vector< vector<AlignedBase> > a;
        for(const string& r: rows) {
            a.push_back(vectorOfAlignedBasesFromString(r));
        }
        vector<AlignedBase> ac = vectorOfAlignedBasesFromString(rows[0]);
        vector< pair<Base, uint64_t> > c;
        for(const char ch: rows[0]) {
            c.push_back(make_pair(Base::fromCharacter(ch), 2UL));
        }
        const vector<uint64_t> w(3, 1);
        const auto aBefore = a;
        const uint64_t repaired = msa1(a, ac, c, w, threshold);
        SHASTA2_ASSERT(repaired == 0);
        SHASTA2_ASSERT(a == aBefore);
    }



    // Adversarial cases, and then many random ones. What is checked in both is
    // the set of invariants that must hold whatever the input:
    //  - every row still has the length of the aligned consensus
    //  - no read was altered, only the columns it occupies
    //  - the consensus is still the aligned consensus with the gaps removed
    //  - the alignment was not made worse
    //  - if nothing was repaired, nothing at all changed
    {
        // Build an alignment from rows, give it a column majority consensus with
        // recognisable coverage, repair it, and check the invariants.
        const auto checkOne = [&](const vector<string>& rows) -> bool
        {
            vector< vector<AlignedBase> > a;
            for(const string& r: rows) {
                a.push_back(vectorOfAlignedBasesFromString(r));
            }
            const vector<uint64_t> w(a.size(), 1);

            vector<AlignedBase> ac;
            vector< pair<Base, uint64_t> > c;
            for(uint64_t j=0; j<a.front().size(); j++) {
                array<uint64_t, 5> t;
                fill(t.begin(), t.end(), 0);
                for(uint64_t i=0; i<a.size(); i++) {
                    t[a[i][j].value] += w[i];
                }
                const auto it = std::ranges::max_element(t);
                const AlignedBase b = AlignedBase::fromInteger(uint64_t(it - t.begin()));
                ac.push_back(b);
                if(not b.isGap()) {
                    c.push_back(make_pair(Base(b), 500 + c.size()));
                }
            }

            const auto aBefore = a;
            const auto cBefore = c;
            vector<string> readsBefore;
            for(const auto& row: a) {
                string s;
                for(const AlignedBase b: row) {
                    if(not b.isGap()) {
                        s.push_back(b.character());
                    }
                }
                readsBefore.push_back(s);
            }
            const uint64_t impureBefore = msa1ImpureColumnCount(
                [&]{ vector<string> v; for(const auto& r: a) v.push_back(msa1ToString(r)); return v; }());

            const uint64_t repaired = msa1(a, ac, c, w, threshold);

            for(const auto& row: a) {
                SHASTA2_ASSERT(row.size() == ac.size());
            }
            for(uint64_t i=0; i<a.size(); i++) {
                string s;
                for(const AlignedBase b: a[i]) {
                    if(not b.isGap()) {
                        s.push_back(b.character());
                    }
                }
                SHASTA2_ASSERT(s == readsBefore[i]);
            }
            {
                string fromAligned;
                for(const AlignedBase b: ac) {
                    if(not b.isGap()) {
                        fromAligned.push_back(b.character());
                    }
                }
                SHASTA2_ASSERT(fromAligned == msa1ToString(c));
            }
            const uint64_t impureAfter = msa1ImpureColumnCount(
                [&]{ vector<string> v; for(const auto& r: a) v.push_back(msa1ToString(r)); return v; }());
            SHASTA2_ASSERT(impureAfter <= impureBefore);
            if(repaired == 0) {
                SHASTA2_ASSERT(a == aBefore);
                SHASTA2_ASSERT(c == cBefore);
            }
            return repaired > 0;
        };

        const string A12(12, 'A'), A11(11, 'A'), A10(10, 'A'), A6(6, 'A');
        const vector< vector<string> > cases = {
            // The pattern flush against each end of the alignment, and filling it.
            {A12 + "G" + A11 + "CGT", A10 + "--G" + A11 + "CGT", A12 + "G" + A11 + "CGT"},
            {"CGT" + A12 + "G" + A11, "CGT" + A10 + "--G" + A11, "CGT" + A12 + "G" + A11},
            {A12 + "G" + A11, A10 + "--G" + A11, A12 + "G" + A11},
            // Two patterns, close enough to merge and far enough to stay apart.
            {A12 + "G" + A11 + "C" + A12 + "G" + A11,
             A10 + "--G" + A11 + "C" + A10 + "--G" + A11,
             A12 + "G" + A11 + "C" + A12 + "G" + A11},
            {A12 + "G" + A11 + string(40, 'C') + A12 + "G" + A11,
             A10 + "--G" + A11 + string(40, 'C') + A10 + "--G" + A11,
             A12 + "G" + A11 + string(40, 'C') + A12 + "G" + A11},
            // Rows that contribute nothing, or much less than the others.
            {A12 + "G" + A11, A10 + "--G" + A11, string(24, '-')},
            {A12 + "G" + A11, A12 + "G" + A11, "-----" + A6 + "G" + A6 + "------"},
            // Very long runs.
            {string(100, 'A') + "G" + string(100, 'A'),
             string(98, 'A') + "--G" + string(100, 'A'),
             string(100, 'A') + "G" + string(100, 'A')},
            // A real substitution where the separating base is, which is not a
            // homopolymer problem and must be left alone.
            {A12 + "G" + A11, A12 + "T" + A11, A12 + "G" + A11},
            // Three long runs in a row, and two patterns sharing a run.
            {A12 + "G" + string(11, 'C') + "T" + string(12, 'G'),
             A10 + "--G" + string(11, 'C') + "T" + string(12, 'G'),
             A12 + "G" + string(11, 'C') + "T" + string(12, 'G')},
            {A12 + "G" + A11 + "G" + A12, A10 + "--G" + A11 + "G" + A12, A12 + "G" + A11 + "G" + A12},
            // Small and degenerate shapes.
            {A12 + "G" + A11, A10 + "--G" + A11},
            {A12 + "G" + A11},
            {A12 + "G" + A11, A12 + "G" + A11, A12 + "G" + A11},
            // Runs sitting on the threshold, a gap column inside a run, a long
            // gap stretch, and every row a different run length.
            {"CC" + string(4, 'A') + "G" + string(4, 'A') + "CC",
             "CC" + string(5, 'A') + "G" + string(4, 'A') + "C",
             "CC" + string(4, 'A') + "G" + string(4, 'A') + "CC"},
            {A6 + "-" + A6 + "G" + A11, A6 + "-" + A6 + "G" + A11, A6 + "A" + A6 + "G" + A11},
            {A12 + "G" + string(10, '-') + A11, A10 + "--G" + string(10, '-') + A11,
             A12 + "G" + string(10, '-') + A11},
            {string(9, 'A') + "---G" + A11, string(10, 'A') + "--G" + A11,
             string(11, 'A') + "-G" + A11, string(12, 'A') + "G" + A11}
        };
        uint64_t repairedCases = 0;
        for(const vector<string>& rows: cases) {
            if(checkOne(rows)) {
                ++repairedCases;
            }
        }
        cout << "Checked " << cases.size() << " adversarial alignments, " <<
            repairedCases << " of them repaired." << endl;
        SHASTA2_ASSERT(repairedCases > 0);



        // Random alignments. The generator makes rows that differ only in the
        // lengths of their long runs, which is the case the repair targets, then
        // sprinkles gaps to make an alignment.
        uint64_t random = 1;
        const auto next = [&random]() {
            random = random * 6364136223846793005ULL + 1442695040888963407ULL;
            return uint64_t(random >> 33);
        };
        uint64_t randomRepaired = 0;
        const uint64_t trialCount = 3000;
        for(uint64_t trial=0; trial<trialCount; trial++) {
            random = 1000003ULL * (trial + 1);

            const uint64_t rowCount = 1 + next() % 6;
            const uint64_t blockCount = 1 + next() % 5;
            vector<string> ungapped(rowCount);
            uint64_t previousBase = 4;
            for(uint64_t b=0; b<blockCount; b++) {
                uint64_t baseIndex = next() % 4;
                if(baseIndex == previousBase) {
                    baseIndex = (baseIndex + 1) % 4;
                }
                previousBase = baseIndex;
                const char base = "ACGT"[baseIndex];
                const bool isLong = (next() % 2) == 0;
                for(uint64_t i=0; i<rowCount; i++) {
                    const uint64_t runLength = isLong ? (5 + next() % 10) : (1 + next() % 2);
                    for(uint64_t k=0; k<runLength; k++) {
                        ungapped[i].push_back(base);
                    }
                }
                if(b + 1 < blockCount) {
                    uint64_t separatorIndex = next() % 4;
                    if(separatorIndex == previousBase) {
                        separatorIndex = (separatorIndex + 1) % 4;
                    }
                    for(uint64_t i=0; i<rowCount; i++) {
                        ungapped[i].push_back("ACGT"[separatorIndex]);
                    }
                    previousBase = separatorIndex;
                }
            }

            uint64_t width = 0;
            for(const string& s: ungapped) {
                width = max(width, uint64_t(s.size()));
            }
            width += next() % 5;
            vector<string> rows;
            for(uint64_t i=0; i<rowCount; i++) {
                string padded = ungapped[i];
                while(padded.size() < width) {
                    const uint64_t p = next() % (padded.size() + 1);
                    padded.insert(padded.begin() + int64_t(p), '-');
                }
                rows.push_back(padded);
            }
            if(checkOne(rows)) {
                ++randomRepaired;
            }
        }
        cout << "Checked " << trialCount << " random alignments, " <<
            randomRepaired << " of them repaired." << endl;
        SHASTA2_ASSERT(randomRepaired > 0);
    }



    // Degenerate inputs.
    {
        vector< vector<AlignedBase> > a;
        vector<AlignedBase> ac;
        vector< pair<Base, uint64_t> > c;
        SHASTA2_ASSERT(msa1(a, ac, c, {}, threshold) == 0);

        // A single row cannot disagree with the consensus, so nothing is found.
        a.push_back(vectorOfAlignedBasesFromString(
            string(12, 'A') + "G" + string(11, 'A')));
        ac = a.front();
        c.clear();
        for(const AlignedBase b: ac) {
            c.push_back(make_pair(Base(b), 1UL));
        }
        const auto aBefore = a;
        SHASTA2_ASSERT(msa1(a, ac, c, {}, threshold) == 0);
        SHASTA2_ASSERT(a == aBefore);
    }

    cout << "testMsa1Repair passed." << endl;
}
