// Shasta.
#include "AnchorPair.hpp"
#include "Anchor.hpp"
#include "HttpServer.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/iterator/function_output_iterator.hpp>

// Standard library.
#include <cmath>



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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;
            SHASTA2_ASSERT(positionInJourneyB >= positionInJourneyA);    // Allow degenerate AnchorPair witn anchorIdA==anchorIdB
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
    SHASTA2_ASSERT(ordinals.size() == orientedReadIds.size());

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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

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

    SHASTA2_ASSERT(it == orientedReadIds.end());
}



// Just return the journey positions.
void AnchorPair::getPositionsInJourneys(
    const Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& positionsInJourneys) const
{
    positionsInJourneys.clear();

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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;

            positionsInJourneys.push_back(make_pair(positionInJourneyA, positionInJourneyB));
        }

        ++itA;
        ++itB;
    }

    SHASTA2_ASSERT(it == orientedReadIds.end());

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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {       // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
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
            SHASTA2_ASSERT(positionB >= positionA);      // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
        }

        ++itA;
        ++itB;
    }

    SHASTA2_ASSERT(it == orientedReadIds.end());

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
        SHASTA2_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {          // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
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
            SHASTA2_ASSERT(positionB > positionA);

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
            minBaseOffset = min(minBaseOffset, offset);
            maxBaseOffset = max(maxBaseOffset, offset);
        }

        ++itA;
        ++itB;
    }

    SHASTA2_ASSERT(it == orientedReadIds.end());

    averageBaseOffset = uint32_t(std::round(double(sumBaseOffset) / double(size())));
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



bool AnchorPair::contains(OrientedReadId orientedReadId) const
{
    const auto it = std::lower_bound(orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
    return it != orientedReadIds.end() and (*it == orientedReadId);
}



// Return the url for the exploreAnchorPair2 page for this AnchorPair.
string AnchorPair::url() const
{
    string s =
        "exploreAnchorPair2?"
        "anchorIdAString=" + HttpServer::urlEncode(anchorIdToString(anchorIdA)) +
        "&anchorIdBString=" + HttpServer::urlEncode(anchorIdToString(anchorIdB)) +
        "&orientedReadIdsString=";

    for(const OrientedReadId orientedReadId: orientedReadIds) {
        s += orientedReadId.getString();
        s += ",";
    }

    return s;
}



void AnchorPair::writeAllHtml(
    ostream& html,
    const Anchors& anchors,
    const Journeys& journeys) const
{
    // Write the summary and oriented reads.
    writeSummaryHtml(html, anchors);
    writeOrientedReadIdsHtml(html, anchors);

    // Get the positions of anchorIdA and anchorIdB in the
    // journeys of all OrientedReadids.
    vector< pair<uint32_t, uint32_t> > positionsInJourneys;
    getPositionsInJourneys(anchors, positionsInJourneys);

    // Write the journey portions between anchorIdA and anchorIdB.
    writeJourneysHtml(html, journeys, positionsInJourneys);

}



void AnchorPair::writeSummaryHtml(ostream& html, const Anchors& anchors) const
{
    const uint64_t offset = getAverageOffset(anchors);
    html <<
        "<table>"
        "<tr><th>Anchor A<td class=centered>" << anchorIdToString(anchorIdA) <<
        "<tr><th>Anchor B<td class=centered>" << anchorIdToString(anchorIdB) <<
        "<tr><th>Coverage<td class=centered>" << size() <<
        "<tr><th>Average offset<td class=centered>" << offset <<
        "</table>";
}



void AnchorPair::writeOrientedReadIdsHtml(ostream& html, const Anchors& anchors) const
{
    vector< pair<Positions, Positions> > positions;
    get(anchors, positions);

    html << "<h3>Oriented reads</h3>";
    html <<
        "<p>"
        "<table>"
        "<tr><th>Oriented<br>read id"
        "<th>Position<br>in journey<br>A<th>Position<br>in journey<br>B<th>Journey<br>offset"
        "<th>OrdinalA<th>OrdinalB<th>Ordinal<br>offset"
        "<th>A middle<br>position"
        "<th>B middle<br>position"
        "<th>Sequence<br>length";

    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];

        const auto& positionsA = positionsAB.first;
        const auto& positionsB = positionsAB.second;

        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << positionsA.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney - positionsA.positionInJourney <<
            "<td class=centered>" << positionsA.ordinal <<
            "<td class=centered>" << positionsB.ordinal <<
            "<td class=centered>" << positionsB.ordinal - positionsA.ordinal <<
            "<td class=centered>" << positionsA.basePosition <<
            "<td class=centered>" << positionsB.basePosition <<
            "<td class=centered>" << positionsB.basePosition - positionsA.basePosition;
    }
    html << "</table>";

}



void AnchorPair::writeJourneysHtml(
    ostream& html,
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys   // As computed by getPositionsInJourneys.
    ) const
{
    html << "<h3>Journey portions within this anchor pair</h3>";
    html << "<table>";

    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey& journey = journeys[orientedReadId];

        const auto& positionInJourney = positionsInJourneys[i];
        const auto positionInJourneyA = positionInJourney.first;
        const auto positionInJourneyB = positionInJourney.second;

        if(html) {
            html << "<tr><th class=centered>" << orientedReadId;
        }
        for(auto position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const AnchorId anchorId = journey[position];
            html << "<td class=centered>" << anchorIdToString(anchorId);
        }
    }

    html << "</table>";

}
