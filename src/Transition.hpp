#pragma once

// Shasta.
#include "Anchor.hpp"
#include "MarkerInfo.hpp"

namespace shasta {
    class Transition;
}



// A transition is the set of oriented reads that visit
// anchorIdB immediately after visiting anchorIdA.
class shasta::Transition {
public:

    Transition(
        const Anchors&,
        AnchorId,
        AnchorId);

    AnchorId anchorIdA;
    AnchorId anchorIdB;

    vector<MarkerInterval> markerIntervals;
};
