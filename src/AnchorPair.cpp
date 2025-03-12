#include "AnchorPair.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta;



AnchorPair::AnchorPair(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB) :
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    // Loop over common oriented reads between these two anchors,
    // for which the positions in journeys differ by 1.

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

        if(itB->positionInJourney == itA->positionInJourney + 1) {
            markerIntervals.push_back(
                MarkerInterval(orientedReadId, itA->ordinal, itB->ordinal));
        }

        ++itA;
        ++itB;
    }
}



// Get a vector of pairs(positionA, positionB), each corresponding
// to one of the MarkerIntervals.
// The positions returned are the midpoint of the markers
// correspoding to anchorIdA and anchorIdB.
void AnchorPair::getPositions(
    const Markers& markers,
    vector< pair<uint32_t, uint32_t> >& positions) const
{
    const uint32_t kHalf = uint32_t(markers.k / 2);
    positions.clear();

    for(const MarkerInterval& markerInterval: markerIntervals) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t positionA = orientedReadMarkers[markerInterval.ordinalA].position + kHalf;
        const uint32_t positionB = orientedReadMarkers[markerInterval.ordinalB].position + kHalf;
        positions.push_back(make_pair(positionA, positionB));
    }
}



// Get a vector of base sequences, each corresponding
// to one of the MarkerIntervals.
// The sequences returned are between the midpoint of the markers
// corresponding to anchorIdA and anchorIdB.
void AnchorPair::getSequences(
    const Markers& markers,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(markers.k / 2);
    sequences.clear();

    for(const MarkerInterval& markerInterval: markerIntervals) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t positionA = orientedReadMarkers[markerInterval.ordinalA].position + kHalf;
        const uint32_t positionB = orientedReadMarkers[markerInterval.ordinalB].position + kHalf;

        sequences.resize(sequences.size() + 1);
        vector<Base>& sequence = sequences.back();
        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(markers.reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}


// Combine the two above in a single call.
void AnchorPair::getPositionsAndSequences(
    const Markers& markers,
    vector< pair<uint32_t, uint32_t> >& positions,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(markers.k / 2);
    positions.clear();
    sequences.clear();

    for(const MarkerInterval& markerInterval: markerIntervals) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t positionA = orientedReadMarkers[markerInterval.ordinalA].position + kHalf;
        const uint32_t positionB = orientedReadMarkers[markerInterval.ordinalB].position + kHalf;
        positions.push_back(make_pair(positionA, positionB));

        sequences.resize(sequences.size() + 1);
        vector<Base>& sequence = sequences.back();
        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(markers.reads.getOrientedReadBase(orientedReadId, position));
        }
    }

}

