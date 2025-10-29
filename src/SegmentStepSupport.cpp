// Shasta.
#include "SegmentStepSupport.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
using namespace shasta;

// Standard library.
#include "algorithm.hpp"



// Get a vector of SegmentStepSupport for a given range of steps
// of an AssemblyGraphEdge. This discard any previous contents of the vector.
void SegmentStepSupport::get(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepBegin,
    uint32_t stepEnd,
    vector<SegmentStepSupport>& v
    )
{
    v.clear();
    for(uint32_t stepId=stepBegin; stepId<stepEnd; stepId++) {
        append(assemblyGraph, e, stepId, v);
    }
}



// Get a vector of SegmentStepSupport for the first stepCount steps
// of an AssemblyGraphEdge. This discard any previous contents of the vector.
void SegmentStepSupport::getInitial(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepCount,
    vector<SegmentStepSupport>& v
    )
{
    const uint32_t totalStepCount = uint32_t(assemblyGraph[e].size());
    const uint32_t begin = 0;
    const uint32_t end = min(stepCount, totalStepCount);

    get(assemblyGraph, e, begin, end, v);
}



// Get a vector of SegmentStepSupport for the last stepCount steps
// of an AssemblyGraphEdge. This discard any previous contents of the vector.
void SegmentStepSupport::getFinal(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepCount,
    vector<SegmentStepSupport>& v
    )
{
    const uint32_t totalStepCount = uint32_t(assemblyGraph[e].size());
    const uint32_t begin =
        ((totalStepCount >= stepCount) ? (totalStepCount - stepCount) : 0);
    const uint32_t end = totalStepCount;

    get(assemblyGraph, e, begin, end, v);

}



// Append to a vector the SegmentStepSupport for a given step
// of an AssemblyGraphEdge.
void SegmentStepSupport::append(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepId,
    vector<SegmentStepSupport>& v
    )
{
    const uint32_t kHalf = uint32_t(assemblyGraph.anchors.k / 2);

    const AssemblyGraphEdge& edge = assemblyGraph[e];
    SHASTA_ASSERT(stepId < uint32_t(edge.size()));
    const AssemblyGraphEdgeStep& step = edge[stepId];
    const AnchorPair& anchorPair = step.anchorPair;

    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        const auto orientedReadMarkers = assemblyGraph.anchors.markers[orientedReadId.getValue()];

        v.emplace_back();
        SegmentStepSupport& stepSupport = v.back();

        stepSupport.e = e;
        stepSupport.stepId = stepId;
        stepSupport.orientedReadId = orientedReadId;

        // Get position information at the left AnchorId of this step.
        const AnchorId anchorIdA = anchorPair.anchorIdA;
        const AnchorMarkerInfo& anchorMarkerInfoA = assemblyGraph.anchors.getAnchorMarkerInfo(anchorIdA, orientedReadId);
        stepSupport.positionInJourneyA = anchorMarkerInfoA.positionInJourney;
        stepSupport.ordinalA = anchorMarkerInfoA.ordinal;
        stepSupport.positionA = orientedReadMarkers[anchorMarkerInfoA.ordinal].position + kHalf;

        // Get position information at the right AnchorId of this step.
        const AnchorId anchorIdB = anchorPair.anchorIdB;
        const AnchorMarkerInfo& anchorMarkerInfoB = assemblyGraph.anchors.getAnchorMarkerInfo(anchorIdB, orientedReadId);
        stepSupport.positionInJourneyB = anchorMarkerInfoB.positionInJourney;
        stepSupport.ordinalB = anchorMarkerInfoB.ordinal;
        stepSupport.positionB = orientedReadMarkers[anchorMarkerInfoB.ordinal].position + kHalf;

    }
}
