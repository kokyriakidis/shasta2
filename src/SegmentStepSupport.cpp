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



void SegmentStepSupport::getInitialFirst(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepCount,
    vector<SegmentStepSupport>& v
    )
{
    getInitial(assemblyGraph, e, stepCount, v);
    keepFirst(v);
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



void SegmentStepSupport::getFinalLast(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint32_t stepCount,
    vector<SegmentStepSupport>& v
    )
{
    getFinal(assemblyGraph, e, stepCount, v);
    keepLast(v);
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



// Given a vector of SegmentStepSupport, for each OrientedReadId keep only the one
// with the largest stepId.
void SegmentStepSupport::keepLast(vector<SegmentStepSupport>& v)
{
    // Sort by OrientedReadId, then by stepId.
    std::ranges::sort(v, std::ranges::less(),
        [](const SegmentStepSupport& s) {return std::tie(s.orientedReadId, s.stepId);}
        );

    // For each streak with the same OrientedReadId, keep the last.
    auto it = v.begin();
    const auto end = v.end();
    auto out = v.begin();
    while(it != end) {
        const OrientedReadId orientedReadId = it->orientedReadId;

        // Find the streak with this OrientedReadId.
        auto streakBegin = it;
        auto streakEnd = streakBegin + 1;
        while((streakEnd != end) and (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }

        // Store the last one in the streak.
        *out++ = *(streakEnd - 1);

        // Prepare to process the next streak;
        it = streakEnd;
    }

    // Only keep the portion we filled in.
    v.resize(out - v.begin());
}



// Given a vector of SegmentStepSupport, for each OrientedReadId keep only the one
// with the smallest stepId.
void SegmentStepSupport::keepFirst(vector<SegmentStepSupport>& v)
{
    // Sort by OrientedReadId, then by stepId.
    std::ranges::sort(v, std::ranges::less(),
        [](const SegmentStepSupport& s) {return std::tie(s.orientedReadId, s.stepId);}
        );

    // For each streak with the same OrientedReadId, keep the last.
    auto it = v.begin();
    const auto end = v.end();
    auto out = v.begin();
    while(it != end) {
        const OrientedReadId orientedReadId = it->orientedReadId;

        // Find the streak with this OrientedReadId.
        auto streakBegin = it;
        auto streakEnd = streakBegin + 1;
        while((streakEnd != end) and (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }

        // Store the last one in the streak.
        *out++ = *streakBegin;

        // Prepare to process the next streak;
        it = streakEnd;
    }

    // Only keep the portion we filled in.
    v.resize(out - v.begin());
}



// Output a vector of SegmentStepSupport to a html table.
void SegmentStepSupport::writeHtml(
    ostream& html,
    const AssemblyGraph& assemblyGraph,
    const vector<SegmentStepSupport>& v)
{
    html << "<table><tr>"
        "<th>Segment"
        "<th>Step"
        "<th>Left<br>AnchorId"
        "<th>Right<br>AnchorId"
        "<th>OrientedReadId"
        "<th>Left<br>position<br>in journey"
        "<th>Right<br>position<br>in journey"
        "<th>Position<br>in journey<br>offset"
        "<th>Left<br>ordinal"
        "<th>Right<br>ordinal"
        "<th>Ordinal<br>offset"
        "<th>Left<br>position"
        "<th>Right<br>position"
        "<th>Position<br>offset";

    for(const SegmentStepSupport& stepSupport: v) {
        const edge_descriptor e = stepSupport.e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const AssemblyGraphEdgeStep& step = edge[stepSupport.stepId];
        html <<
            "<tr>" <<
            "<td class=centered>" << edge.id <<
            "<td class=centered>" << stepSupport.stepId <<
            "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
            "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<

            "<td class=centered>" << stepSupport.orientedReadId <<

            "<td class=centered>" << stepSupport.positionInJourneyA <<
            "<td class=centered>" << stepSupport.positionInJourneyB <<
            "<td class=centered>" << stepSupport.positionInJourneyOffset() <<

            "<td class=centered>" << stepSupport.ordinalA <<
            "<td class=centered>" << stepSupport.ordinalB <<
            "<td class=centered>" << stepSupport.ordinalOffset() <<

            "<td class=centered>" << stepSupport.positionA <<
            "<td class=centered>" << stepSupport.positionB <<
            "<td class=centered>" << stepSupport.positionOffset();
    }
    html << "</table>";

}


