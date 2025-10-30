// Shasta.
#include "SegmentStepSupport.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta;

// Standard library.
#include "algorithm.hpp"
#include <iomanip>



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
        "<th>Read<br>length"
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

            "<td class=centered>" << assemblyGraph.anchors.reads.getReadSequenceLength(stepSupport.orientedReadId.getReadId()) <<
            "<td class=centered>" << stepSupport.positionA <<
            "<td class=centered>" << stepSupport.positionB <<
            "<td class=centered>" << stepSupport.positionOffset();
    }
    html << "</table>";

}






// Analyze SegmentStepSupport by comparing
// the final representative region of e0
// with the initial representative region of e1.
SegmentPairInformation SegmentStepSupport::analyzeSegmentPair(
    ostream& html,
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e0,
    edge_descriptor e1,
    uint32_t representativeRegionStepCount
    )
{
    const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
    const AssemblyGraphEdge& edge1 = assemblyGraph[e1];

    // Get SegmentStepSupport for the final representative region of e0,
    // then for each OrientedReadId keep only the one with the largest stepId.
    vector<SegmentStepSupport> support0;
    SegmentStepSupport::getFinalLast(assemblyGraph, e0, representativeRegionStepCount, support0);

    // Get SegmentStepSupport for the initial representative region of e1,
    // then for each OrientedReadId keep only the one with the largest stepId.
    vector<SegmentStepSupport> support1;
    SegmentStepSupport::getInitialFirst(assemblyGraph, e1, representativeRegionStepCount, support1);

    if(html) {
        html << "<h3>Final support for " << edge0.id << "</h3>";
        SegmentStepSupport::writeHtml(html, assemblyGraph, support0);

        html << "<h3>Initial support for " << edge1.id << "</h3>";
        SegmentStepSupport::writeHtml(html, assemblyGraph, support1);
    }



    // Use the common OrientedReadIds to estimate offset.
    auto it0 = support0.begin();
    const auto end0 = support0.end();
    auto it1 = support1.begin();
    const auto end1 = support1.end();
    uint32_t commonCount = 0;
    int64_t offsetSum = 0;
    if(html) {
        html << "<h3>Common support</h3>"
            "<table><tr>"
            "<th>OrientedReadId<th>Estimated<br>offset<br>between<br>segments";
    }
    while((it0 != end0) and (it1 != end1)) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
            continue;
        }
        ++commonCount;
        const int32_t offset = estimateOffset(assemblyGraph, *it0, *it1);
        offsetSum += offset;

        if(html) {
            html << "<tr><td class=centered>" << it0->orientedReadId <<
                "<td class=centered>" << offset;
        }

        ++it0;
        ++it1;
    }
    if(html) {
        html << "</table>";
    }

    if(commonCount == 0) {
        if(html) {
            html << "<br>No common coverage is available to estimate segment offset.";
        }
        return SegmentPairInformation();
    }

    const int32_t segmentOffset = int32_t(std::round(double(offsetSum) / double(commonCount)));
    if(html) {
        html << "<br>Common count " << commonCount <<
            "<br>Estimated segment offset " << segmentOffset;
    }

    const uint64_t intersectionCount = commonCount;
    const uint64_t unionCount = support0.size() + support1.size() - intersectionCount;
    const double jaccard = double(intersectionCount) / double(unionCount);
    if(html) {
        html << "<br>Jaccard similarity " << std::fixed << std::setprecision(2) << jaccard;
    }



    // Support present in e0 but not e1.
    uint32_t missing0 = 0;
    if(html) {
        html << "<h3>Support present in " << edge0.id <<
            " but not in " << edge1.id <<
            "</h3>"
            "<table><tr>"
            "<th>OrientedReadId"
            "<th>Status";
    }
    it0 = support0.begin();
    it1 = support1.begin();
    while(it0 != end0) {
        while((it1 != end1) and (it1->orientedReadId < it0->orientedReadId)) {
            ++it1;
        }
        if((it1 != end1) and (it1->orientedReadId == it0->orientedReadId)) {
            ++it0;
            continue;
        }

        // Figure out if this read is too short to appear in the first step of e1.
        int32_t position = it0->positionB;
        for(uint32_t stepId=it0->stepId+1; stepId<edge0.size(); stepId++) {
            const uint64_t stepOffset = (edge0.wasAssembled ? edge0[stepId].sequence.size() : edge0[stepId].offset);
            position += int32_t(stepOffset);
        }
        position += segmentOffset;
        const uint64_t step1Offset = (edge1.wasAssembled ? edge1.front().sequence.size() : edge1.front().offset);
        position += int32_t(step1Offset);
        const int32_t readLength = int32_t(assemblyGraph.anchors.reads.getReadSequenceLength(it0->orientedReadId.getReadId()));
        const bool isShort = readLength < position;

        if(not isShort) {
            ++missing0;
        }


        if(html) {
            html << "<tr><td class=centered>" << it0->orientedReadId <<
                "<td class=centered>";
            html << (isShort ? "Short" : "Missing");
        }
        ++it0;
    }
    if(html) {
        html << "</table>";
        html << "<br>" << missing0 << " missing.";
    }



    // Support present in e1 but not e0.
    uint32_t missing1 = 0;
    if(html) {
        html << "<h3>Support present in " << edge1.id <<
            " but not in " << edge0.id <<
            "</h3>"
            "<table><tr>"
            "<th>OrientedReadId"
            "<th>Status";
    }
    it0 = support0.begin();
    it1 = support1.begin();
    while(it1 != end1) {
        while((it0 != end0) and (it0->orientedReadId < it1->orientedReadId)) {
            ++it0;
        }
        if((it0 != end0) and (it0->orientedReadId == it1->orientedReadId)) {
            ++it1;
            continue;
        }

        // Figure out if this read is too short to appear in the last step of e0.
        int32_t position = it1->positionA;
        for(uint32_t stepId=0; stepId<it1->stepId; stepId++) {
            const uint64_t stepOffset = (edge1.wasAssembled ? edge1[stepId].sequence.size() : edge1[stepId].offset);
            position -= int32_t(stepOffset);
        }
        position -= segmentOffset;
        const uint64_t step0Offset = (edge0.wasAssembled ? edge0.back().sequence.size() : edge0.back().offset);
        position -= int32_t(step0Offset);
        const bool isShort = position < 0;

        if(not isShort) {
            ++missing1;
        }

        if(html) {
            html << "<tr><td class=centered>" << it1->orientedReadId <<
                "<td class=centered>";
            html << (isShort ? "Short" : "Missing");
        }

        ++it1;
    }
    if(html) {
        html << "</table>";
        html << "<br>" << missing1 << " missing.";
    }

    const uint64_t correctedUnionSize = missing0 + missing1 + commonCount;
    const double correctedJaccard = double(intersectionCount) / double(correctedUnionSize);
    if(html) {
        html << "<br>Corrected jaccard similarity " << std::fixed << std::setprecision(2) << correctedJaccard;
    }


    SegmentPairInformation info;
    info.commonCount = commonCount;
    info.missing0 = missing0;
    info.missing1 = missing1;
    info.segmentOffset = segmentOffset;
    info.jaccard = jaccard;
    info.correctedJaccard = correctedJaccard;
    return info;
}



// Estimate the offset between two segments using a SegmentStepSupport
// on the first segment and a SegmentStepSupport on the second segment,
// both for the same OrientedReadId.
int32_t SegmentStepSupport::estimateOffset(
    const AssemblyGraph& assemblyGraph,
    SegmentStepSupport& s0,
    SegmentStepSupport& s1)
{
    SHASTA_ASSERT(s0.orientedReadId == s1.orientedReadId);

    const edge_descriptor e0 = s0.e;
    const edge_descriptor e1 = s1.e;

    const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
    const AssemblyGraphEdge& edge1 = assemblyGraph[e1];

    const uint32_t position0 = s0.positionB;
    const uint32_t position1 = s1.positionA;

    // We can't simply compute position1 - position0.
    // We have to subtract the offset between position0 and the end of e0,
    // and the offset between the beginning of e1 and position1.

    // Compute the offset between position0 and the end of e0.
    uint64_t offset0 = 0;
    for(uint32_t stepId=s0.stepId+1; stepId<edge0.size(); stepId++) {
        if(edge0.wasAssembled) {
            offset0 += edge0[stepId].sequence.size();
        } else {
            offset0 += edge0[stepId].offset;
        }
    }

    // Compute the offset between the beginning of e1 and position1.
    uint64_t offset1 = 0;
    for(uint32_t stepId=0; stepId<s1.stepId; stepId++) {
        if(edge1.wasAssembled) {
            offset1 += edge1[stepId].sequence.size();
        } else {
            offset1 += edge1[stepId].offset;
        }
    }

    return int32_t(position1) - int32_t(position0) - int32_t(offset0) - int32_t(offset1);
}
