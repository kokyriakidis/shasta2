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



// An AnchorPair is a set of OrientedReadIds that visit anchorIdA
// and then, later in their Journey, anchorIdB.
// That is, each of the OrientedReadIds appear both in anchorIdA and anchorIdB,
// and the position at anchorIdB is greater than the position at anchorIdA.
// The AnchorPair can use a subset of all possible OrientedReadIds
// that satisfy the above.
class shasta::AnchorPair {
public:

    // The constructor creates an AnchorPair between anchorIdA and anchorIdB.
    // - If adjacentInJourney is false, it includes all OrientedReadIds
    //  that visit anchorIdA and then, later in their Journey, anchorIdB.
    // - If adjacentInJourney is true, only OrientedReadIds that visit anchorIdB
    //   immediately after visiting anchorIdA are included.
    //   That is, the journey offset between anchorIdA and anchorIdB
    //   for OrientedReadIds that are included must be exactly 1.
    AnchorPair(
        const Anchors&,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        bool adjacentInJourney);

    AnchorPair() {}

    AnchorPair(const AnchorPair& that) :
        anchorIdA(that.anchorIdA),
        anchorIdB(that.anchorIdB),
        orientedReadIds(that.orientedReadIds)
    {}

    // Copy from another AnchorPair, but excluding some OrientedReadIds.
    AnchorPair(
        const AnchorPair&,
        const vector<OrientedReadId>& excludedOrientedReadIds);

    // "Join" constructor from two AnchorPairs.
    // This constructs a new AnchorPair as follows:
    // - anchorIdA is the same as anchorPair0.anchorIdA.
    // - anchorIdB is the same as anchorPair1.anchorIdB.
    // - orientedReadIds are the intersection of
    //   anchorPair0.orientedReadIds and anchorPair1.orientedReadIds,
    //   with the additional requirement that the new journey offset
    //   is positive. That is, each OrientedReadId of the new AnchorPair
    //   visits anchorIdB after visiting the new anchorIdA.
    // This constructor is used in detangling.
    AnchorPair(
        const Anchors&,
        const AnchorPair& anchorPair0,
        const AnchorPair& anchorPair1);

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

    uint32_t getAverageOffset(const Anchors&) const;
    void getOffsets(
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

    // Just return the ordinals.
    void getOrdinals(const Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Split the AnchorPair into one or more AnchorPairs with consistent offsets.
    // In the resulting AnchorPairs, if the position offsets are sorted in
    // increasing order, any two adjacent offsets D0 and D1
    // will satisfy D1 - D0 <= aDrift + bDrift * (D0 + D1) / 2.
    void split(
        const Anchors&,
        double aDrift,
        double bDrift,
        vector<AnchorPair>&
        ) const;

    // This returns true if a call to split with the same arguments would not split this Anchor.
    // The second and third areguments are work vectors added as arguments for performancew,
    bool isConsistent(
        const Anchors&,
        double aDrift,
        double bDrift,
        vector< pair<Positions, Positions> >&,
        vector<uint64_t>&) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIds;
    }
};
