#include "Transition.hpp"
using namespace shasta;


Transition::Transition(
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
