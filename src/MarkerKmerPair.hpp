#pragma once

// Shasta.
#include "invalid.hpp"
#include "Kmer.hpp"
#include "MarkerInfo.hpp"
#include "Markers.hpp"
#include "Reads.hpp"

// Standard library.
#include <map>
#include <vector.hpp>

namespace shasta {
    class MarkerKmerPair;

    class MarkerKmers;
    class Markers;
}



// A MarkerKmerPair between two MarkerKmers kmer0 and kmer1 is defined using
// OrientedReadIds that appear in kmer0 at ordinal0 and in kmer1 at ordinal1 and such that:
// - ordinal0 < ordinal1
// - position1 - position0 <= maxPositionOffset
// - The ReadId appears exactly once in both kmer0 and kmer1.
//   This implies that the OrientedReadId also appears exactly once in both kmer0 and kmer1.
class shasta::MarkerKmerPair {
public:

    MarkerKmerPair(
        const MarkerKmers&,
        const Kmer& kmer0,
        const Kmer& kmer1,
        uint64_t maxPositionOffset);

    // The left and right Kmer of the pair.
    Kmer kmer0;
    Kmer kmer1;

    // The MarkerInfos of kmer0 and kmer1, excluding ReadIds that appear
    // more than once in kmer0 or kmer1.
    vector<MarkerInfo> markerInfos0;
    vector<MarkerInfo> markerInfos1;
    void getMarkerInfos(const MarkerKmers&);



    // The distinct sequences of the common oriented reads stored
    // below between the midpoint of the two Kmers.
    class SequenceInfo {
    public:

        // The indexes of the oriented reads in the commonOrientedReads vector
        // that have this sequence.
        vector<uint64_t> orientedReadIndexes;
        uint64_t coverage() const
        {
            return orientedReadIndexes.size();
        }

        // The rank of this sequence in order of decreasing coverage.
        // The sequence with the most coverage has rank 0.
        uint64_t rank = invalid<uint64_t>;

    };
    using SequenceMap = std::map<vector<Base>, SequenceInfo>;
    SequenceMap sequenceMap;
    vector<SequenceMap::iterator> sequencesByRank;
    void gatherSequences(const Reads&);
    void rankSequences();

    // Compute the alignment.
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    void align();

    // Return the edit distance between two sequences with given rank.
    uint64_t editDistance(uint64_t rank0, uint64_t rank1) const;


    // Common oriented reads that have kmer0 at ordinal0
    // and kmer1 at ordinal1, with ordinal0 < ordinal1.
    class CommonOrientedRead {
    public:
        OrientedReadId orientedReadId;

        uint32_t ordinal0;
        uint32_t ordinal1;
        uint32_t ordinalOffset() const
        {
            return ordinal1 - ordinal0;
        }

        // Positions of the mid point of these markers.
        uint32_t position0;
        uint32_t position1;
        uint32_t positionOffset() const
        {
            return position1 - position0;
        }

        void getSequence(const Reads&, vector<Base>&) const;

        // An iterator pointing to the sequenceMap entry
        // corresponding to the sequence of this oriented read
        // between the midpoints of the two markers.
        SequenceMap::const_iterator sequenceMapIterator;

        bool operator<(const CommonOrientedRead& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<CommonOrientedRead> commonOrientedReads;
    void gatherCommonOrientedReads(const Markers&, uint64_t maxPositionOffset);

    void writeSummary(ostream& html, uint64_t k) const;
    void writeSequences(ostream& html) const;
    void writeCommonOrientedReads(ostream& html) const;
    void writeAlignment(ostream& html) const;
    void writePairAlignmentDistances(ostream& html) const;
};
