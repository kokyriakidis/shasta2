// Shasta.
#include "AnchorPair.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/iterator/function_output_iterator.hpp>

// Standard library.
#include "stdexcept.hpp"



AnchorPair::AnchorPair(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    bool adjacentInJourney) :
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    // Loop over common oriented reads between these two anchors.
    // If adjacentInJourney is false, the journey offset is required to be positive.
    // If adjacentInJourney is true, the journey offset is required to be exactly 1.

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        if(adjacentInJourney) {
            if(itB->positionInJourney == itA->positionInJourney + 1) {
                orientedReadIds.push_back(orientedReadId);
            }
        } else {
            if(itB->positionInJourney > itA->positionInJourney) {
                orientedReadIds.push_back(orientedReadId);
            }
        }

        ++itA;
        ++itB;
    }
}



// Copy from another AnchorPair, but excluding some OrientedReadIds.
AnchorPair::AnchorPair(
    const AnchorPair& that,
    const vector<OrientedReadId>& excludedOrientedReadIds) :
    anchorIdA(that.anchorIdA),
    anchorIdB(that.anchorIdB)
{
    std::set_difference(
        that.orientedReadIds.begin(), that.orientedReadIds.end(),
        excludedOrientedReadIds.begin(), excludedOrientedReadIds.end(),
        back_inserter(orientedReadIds));
}


// Get positions in journey, ordinals, and base positions
// for each of the two reads and for each of the two anchors.
// The positions returned are the midpoint of the markers
// corresponding to anchorIdA and anchorIdB.
void AnchorPair::get(
    const Anchors& anchors,
    vector< pair<Positions, Positions> >& positions) const
{

    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);
    positions.clear();

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;
            SHASTA_ASSERT(positionInJourneyB >= positionInJourneyA);    // Allow degenerate AnchorPair witn anchorIdA==anchorIdB
            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;

            positions.push_back(make_pair(
                Positions(positionInJourneyA, ordinalA, positionA),
                Positions(positionInJourneyB, ordinalB, positionB)
                ));
        }

        ++itA;
        ++itB;
    }
}



// Remove from the AnchorPair OrientedReadIds that have negative offsets.
void AnchorPair::removeNegativeOffsets(const Anchors& anchors)
{
    vector< pair<uint32_t, uint32_t> > ordinals;
    getOrdinals(anchors, ordinals);
    SHASTA_ASSERT(ordinals.size() == orientedReadIds.size());

    vector<OrientedReadId> newOrientedReadIds;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const auto& p = ordinals[i];
        if(p.second >= p.first) {
            newOrientedReadIds.push_back(orientedReadIds[i]);
        }
    }

    orientedReadIds.swap(newOrientedReadIds);
}




// Same as the above, but also returns compute the sequences.
void AnchorPair::get(
    const Anchors& anchors,
    vector< pair<Positions, Positions> >& positions,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);
    const Reads& reads = anchors.markers.reads;

    get(anchors, positions);

    sequences.clear();
    sequences.resize(orientedReadIds.size());

    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];
        vector<Base>& sequence = sequences[i];

        const uint32_t positionA = positionsAB.first.basePosition + kHalf;
        const uint32_t positionB = positionsAB.second.basePosition + kHalf;

        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}



// Same as the above, but only compute the ordinals.
void AnchorPair::getOrdinals(
    const Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& ordinals) const
{
    ordinals.clear();

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;

            ordinals.push_back(make_pair(ordinalA, ordinalB));
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());
}



// This finds AnchorPairs as follows:
// - anchorIdA is as specified.
// - Coverage is at least minCoverage.
// - All oriented reads have a journey offset equal to 1.
void AnchorPair::createChildren(
    const Anchors& anchors,
    const Journeys& journeys,
    AnchorId anchorIdA,
    uint64_t minCoverage,
    vector<AnchorPair>& anchorPairs
    )
{
    // Find possible choices for anchorIdB.
    vector<AnchorId> anchorIdsB;
    vector<uint64_t> coverage;
    anchors.findChildren(journeys, anchorIdA, anchorIdsB, coverage, minCoverage);

    anchorPairs.clear();
    for(const AnchorId anchorIdB: anchorIdsB) {
        anchorPairs.emplace_back(anchors, anchorIdA, anchorIdB, true);
    }
}


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
AnchorPair::AnchorPair(
    const Anchors& anchors,
    const AnchorPair& anchorPair0,
    const AnchorPair& anchorPair1) :
        anchorIdA(anchorPair0.anchorIdA),
        anchorIdB(anchorPair1.anchorIdB)
{
    vector< pair<Positions, Positions> > positions0;
    vector< pair<Positions, Positions> > positions1;
    anchorPair0.get(anchors, positions0);
    anchorPair1.get(anchors, positions1);

    const uint64_t n0 = anchorPair0.size();
    const uint64_t n1 = anchorPair1.size();
    SHASTA_ASSERT(positions0.size() == n0);
    SHASTA_ASSERT(positions1.size() == n1);

    // Joint loop over OrientedReadIds of the two AnchorPairs.
    uint64_t i0 = 0;
    uint64_t i1 = 0;
    while((i0 != n0) and (i1 != n1)) {
        const OrientedReadId orientedReadId0 = anchorPair0.orientedReadIds[i0];
        const OrientedReadId orientedReadId1 = anchorPair1.orientedReadIds[i1];

        if(orientedReadId0 < orientedReadId1) {
            ++i0;
            continue;
        }

        if(orientedReadId1 < orientedReadId0) {
            ++i1;
            continue;
        }

        SHASTA_ASSERT(orientedReadId0 == orientedReadId1);
        const OrientedReadId orientedReadId = orientedReadId0;

        // We found an OrientedReadId common between anchorPair0 and anchorPair1.
        // We also have to check the journey offsets.
        if(positions0[i0].first.positionInJourney < positions1[i1].second.positionInJourney) {
            orientedReadIds.push_back(orientedReadId);
        }

        ++i0;
        ++i1;
    }
}



uint32_t AnchorPair::getAverageOffset(const Anchors& anchors) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);

    uint64_t sumBaseOffset = 0;

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB <= ordinalA) {
                throw runtime_error(
                    "Order violation at anchor pair " +
                    anchorIdToString(anchorIdA) + " " +
                    anchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            SHASTA_ASSERT(positionB > positionA);

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());

    return uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



void AnchorPair::getOffsets(
    const Anchors& anchors,
    uint32_t& averageBaseOffset,
    uint32_t& minBaseOffset,
    uint32_t& maxBaseOffset) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);

    uint64_t sumBaseOffset = 0;
    minBaseOffset = std::numeric_limits<uint32_t>::max();
    maxBaseOffset = 0;

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB <= ordinalA) {
                throw runtime_error(
                    "Order violation at anchor pair " +
                    anchorIdToString(anchorIdA) + " " +
                    anchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            SHASTA_ASSERT(positionB > positionA);

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
            minBaseOffset = min(minBaseOffset, offset);
            maxBaseOffset = max(maxBaseOffset, offset);
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());

    averageBaseOffset = uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



// Split the AnchorPair into one or more AnchorPairs with consistent offsets.
// In the resulting AnchorPairs, if the position offsets are sorted in
// increasing order, any two adjacent offsets D0 and D1
// will satisfy D1 - D0 <= aDrift + bDrift * (D0 + D1) / 2.
void AnchorPair::split(
    const Anchors& anchors,
    double aDrift,
    double bDrift,
    vector<AnchorPair>& newAnchorPairs) const
{
    vector< pair<Positions, Positions> > positions;
    get(anchors, positions);
    const uint64_t n = orientedReadIds.size();
    SHASTA_ASSERT(positions.size() == n);

    // Gather pairs(index, offset) where index is the index
    // in the OrientedReadIds, vector.
    vector< pair<uint64_t, uint64_t> > offsets;
    for(uint64_t i=0; i<n; i++) {
        const uint32_t positionA = positions[i].first.basePosition;
        const uint32_t positionB = positions[i].second.basePosition;
        SHASTA_ASSERT(positionB > positionA);
        const uint64_t offset = positionB - positionA;
        offsets.push_back(make_pair(i, offset));
    }
    sort(offsets.begin(), offsets.end(), OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Find places where we have to split.
    vector<uint64_t> splitPoints;
    splitPoints.push_back(0);
    for(uint64_t i1=1; i1<n; i1++) {
        const uint64_t i0 = i1 - 1;
        const double offset0 = double(offsets[i0].second);
        const double offset1 = double(offsets[i1].second);
        if(offset1 - offset0 > aDrift + .5 * bDrift  * (offset1 + offset0)) {
            splitPoints.push_back(i1);
        }
    }
    splitPoints.push_back(n);

    // Each interval between split points generates a new AnchorPair.
    newAnchorPairs.clear();
    newAnchorPairs.resize(splitPoints.size() - 1);
    for(uint64_t i=0; i<splitPoints.size() -1 ; i++) {
        const uint64_t j0 = splitPoints[i];
        const uint64_t j1 = splitPoints[i + 1];

        newAnchorPairs[i].anchorIdA = anchorIdA;
        newAnchorPairs[i].anchorIdB = anchorIdB;
        for(uint64_t j=j0; j!=j1; j++) {
            newAnchorPairs[i].orientedReadIds.push_back(orientedReadIds[offsets[j].first]);
        }
        sort(newAnchorPairs[i].orientedReadIds.begin(), newAnchorPairs[i].orientedReadIds.end());
    }


    // Sort them by decreasing coverage.
    class SortHelper {
    public:
        bool operator() (const AnchorPair& x, const AnchorPair& y) const
        {
            return x.orientedReadIds.size() > y.orientedReadIds.size();
        }
    };
    sort(newAnchorPairs.begin(), newAnchorPairs.end(), SortHelper());

}



// This returns true if a call to split with the same arguments would not split this Anchor.
// The second and third areguments are work vectors added as arguments for performancew,
bool AnchorPair::isConsistent(
    const Anchors& anchors,
    double aDrift,
    double bDrift,
    vector< pair<Positions, Positions> >& positions,
    vector<uint64_t>& offsets) const
{

    get(anchors, positions);
    const uint64_t n = orientedReadIds.size();
    SHASTA_ASSERT(positions.size() == n);

    // Gather offsets.
    offsets.clear();
    offsets.resize(n);
    for(uint64_t i=0; i<n; i++) {
        const uint32_t positionA = positions[i].first.basePosition;
        const uint32_t positionB = positions[i].second.basePosition;
        SHASTA_ASSERT(positionB > positionA);
        const uint64_t offset = positionB - positionA;
        offsets[i] = offset;
    }
    sort(offsets.begin(), offsets.end());

    for(uint64_t i1=1; i1<n; i1++) {
        const uint64_t i0 = i1 - 1;
        const double offset0 = double(offsets[i0]);
        const double offset1 = double(offsets[i1]);
        if(offset1 - offset0 > aDrift + .5 * bDrift  * (offset1 + offset0)) {
            return false;
        }
    }

    return true;
}



// Count OrientedReadIds in common with another AnchorPair.
uint64_t AnchorPair::countCommon(const AnchorPair& y) const
{
    const AnchorPair& x = *this;

    uint64_t n = 0;
    auto counter = [&n](auto){++n;};

    std::set_intersection(
        x.orientedReadIds.begin(),
        x.orientedReadIds.end(),
        y.orientedReadIds.begin(),
        y.orientedReadIds.end(),
        boost::make_function_output_iterator(counter)
    );

    return n;

}
