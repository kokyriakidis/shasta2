#pragma once

// Shasta.
#include "Anchor.hpp"
#include "MarkerInfo.hpp"

namespace shasta {
    class AnchorPair;
}



// An AnchorPair is a set of OrientedReadIds that visit AnchorIdA
// and then later in their Journey, AnchorIdB.
class shasta::AnchorPair {
public:

    AnchorPair(
        const Anchors&,
        AnchorId,
        AnchorId);

    AnchorId anchorIdA;
    AnchorId anchorIdB;

    vector<MarkerInterval> markerIntervals;

    // Get a vector of pairs(positionA, positionB), each corresponding
    // to one of the MarkerIntervals.
    // The positions returned are the midpoint of the markers
    // corresponding to anchorIdA and anchorIdB.
    void getPositions(
        const Markers&,
        vector< pair<uint32_t, uint32_t> >& positions) const;

    // Get a vector of base sequences, each corresponding
    // to one of the MarkerIntervals.
    // The sequeences returned are between the midpoint of the markers
    // corresponding to anchorIdA and anchorIdB.
    void getSequences(
        const Markers&,
        vector< vector<Base> >&) const;

    // Combine the two above in a single call.
    void getPositionsAndSequences(
        const Markers&,
        vector< pair<uint32_t, uint32_t> >& positions,
        vector< vector<Base> >&) const;
};
