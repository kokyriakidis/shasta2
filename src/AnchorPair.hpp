#pragma once

// Shasta.
#include "shastaTypes.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Standard library.
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {
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
class shasta2::AnchorPair {
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

    // Just return the positions in journeys.
    void getPositionsInJourneys(const Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Count OrientedReadIds in common with another AnchorPair.
    uint64_t countCommon(const AnchorPair&) const;

    // Remove from the AnchorPair OrientedReadIds that have negative offsets.
    void removeNegativeOffsets(const Anchors&);

    bool contains(OrientedReadId) const;

    // Return the url for the exploreAnchorPair2 page for this AnchorPair.
    string url() const;

    // Html output.
    void writeAllHtml(ostream&, const Anchors&, const Journeys&) const;
    void writeSummaryHtml(ostream&, const Anchors&) const;
    void writeOrientedReadIdsHtml(ostream&, const Anchors&) const;
    void writeJourneysHtml(
        ostream&,
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys   // As computed by getPositionsInJourneys.
        ) const;



    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIds;
    }
};
