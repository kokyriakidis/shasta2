#pragma once

// Shasta.
#include "MarkerInfo.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include <cstdint.hpp>
#include <utility.hpp>
#include <vector.hpp>

namespace shasta {
    class AnchorPair;
    class Anchors;
    class Base;
}



// An AnchorPair is a set of OrientedReadIds that visit AnchorIdA
// and then later in their Journey, AnchorIdB.
// That is, each of the OrientedReadIds appear both anchorIdA and anchorIdB,
// and the position at anchorIdB is greater than the position at anchorIdA.
class shasta::AnchorPair {
public:

    // The constructor creates an AnchorPair that includes common
    // oriented reads between anchorIdA and anchorIdB that visit
    // anchorIdB immediately after anchorIdA.
    AnchorPair(
        const Anchors&,
        AnchorId,
        AnchorId);

    AnchorId anchorIdA;
    AnchorId anchorIdB;

    vector<OrientedReadId> orientedReadIds;
    uint64_t size() const
    {
        return orientedReadIds.size();
    }

    // Get positions in journey, ordinals, base positions
    // for each of the two reads and for each of the two anchors.
    // The positions returned are the midpoint of the markers
    // corresponding to anchorIdA and anchorIdB.
    class Positions {
    public:
        uint32_t positionInJourney;
        uint32_t ordinal;
        uint32_t basePosition;
        Positions(
            uint32_t positionInJourney,
            uint32_t ordinal,
            uint32_t basePosition) :
            positionInJourney(positionInJourney),
            ordinal(ordinal),
            basePosition(basePosition)
        {}
    };
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions) const;

    // Same as the above, but also returns compute the sequences.
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions,
        vector< vector<Base> >&) const;
};
