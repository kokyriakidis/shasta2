#pragma once

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include <cstdint.hpp>
#include <utility.hpp>
#include <vector.hpp>

namespace shasta {
    class AnchorPair;
    class Anchors;
    class Base;
    class Journeys;
}



// An AnchorPair is a set of OrientedReadIds that visit AnchorIdA
// and then later in their Journey, AnchorIdB.
// That is, each of the OrientedReadIds appear both anchorIdA and anchorIdB,
// and the position at anchorIdB is greater than the position at anchorIdA.
class shasta::AnchorPair {
public:

    // The constructor creates an AnchorPair such that:
    // - anchorIdA and anchorIdB are as specified.
    // - All oriented reads have a journey offset equal to 1.
    AnchorPair(
        const Anchors&,
        AnchorId anchorIdA,
        AnchorId anchorIdB);

    AnchorPair() {}

    // This finds AnchorPairs as follows:
    // - anchorIdA is as specified.
    // - Coverage is at least minCoverage.
    // - All oriented reads have a journey offset equal to 1.
    static void createChildren(
        const Anchors&,
        const Journeys&,
        AnchorId anchorIdA,
        uint64_t minCoverage,
        vector<AnchorPair>&
        );


    AnchorId anchorIdA = invalid<AnchorId>;
    AnchorId anchorIdB = invalid<AnchorId>;

    vector<OrientedReadId> orientedReadIds;
    uint64_t size() const
    {
        return orientedReadIds.size();
    }

    void getOffsetStatistics(
        const Anchors&,
        uint32_t& averageBaseOffset,
        uint32_t& minBaseOffset,
        uint32_t& maxBaseOffset) const;

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

        template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
        {
            ar & positionInJourney;
            ar & ordinal;
            ar & basePosition;
        }
    };
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions) const;

    // Same as the above, but also returns compute the sequences.
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions,
        vector< vector<Base> >&) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIds;
    }
};
