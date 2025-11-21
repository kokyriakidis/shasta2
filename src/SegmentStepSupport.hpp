#pragma once

// A SegmentStepSupport object contains detailed information
// about one OrientedReadId that appears in an AssemblyGraphEdgeStep.

#include "AssemblyGraph.hpp"

namespace shasta2 {
    class SegmentStepSupport;
    class SegmentPairInformation;
}



class shasta2::SegmentPairInformation {
public:
    uint64_t commonCount = 0;
    uint64_t missing0 = 0;
    uint64_t missing1 = 0;
    int32_t segmentOffset = invalid<uint32_t>;
    double jaccard = 0.;
    double correctedJaccard = 0.;
};



class shasta2::SegmentStepSupport {
public:
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    edge_descriptor e;
    uint32_t stepId;

    OrientedReadId orientedReadId;

    // Position information at the left AnchorId of this step.
    uint32_t positionInJourneyA;
    uint32_t ordinalA;
    uint32_t positionA;

    // Position information at the right AnchorId of this step.
    uint32_t positionInJourneyB;
    uint32_t ordinalB;
    uint32_t positionB;

    uint32_t positionInJourneyOffset() const
    {
        return positionInJourneyB - positionInJourneyA;
    }

    uint32_t ordinalOffset() const
    {
        return ordinalB - ordinalA;
    }

    uint32_t positionOffset() const
    {
        return positionB - positionA;
    }

    // Get a vector of SegmentStepSupport for a given range of steps
    // of an AssemblyGraphEdge. This discard any previous contents of the vector.
    static void get(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepBegin,
        uint32_t stepEnd,
        vector<SegmentStepSupport>&
        );

    // Get a vector of SegmentStepSupport for the first stepCount steps
    // of an AssemblyGraphEdge. This discard any previous contents of the vector.
    static void getInitial(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<SegmentStepSupport>&
        );
    // Same, followed by a call to keepFirst.
    static void getInitialFirst(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<SegmentStepSupport>&
        );

    // Get a vector of SegmentStepSupport for the last stepCount steps
    // of an AssemblyGraphEdge. This discard any previous contents of the vector.
    static void getFinal(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<SegmentStepSupport>&
        );
    // Same, followed by a call to keepLast.
    static void getFinalLast(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<SegmentStepSupport>&
        );

    // Append to a vector the SegmentStepSupport for a given step
    // of an AssemblyGraphEdge.
    static void append(
        const AssemblyGraph&,
        edge_descriptor,
        uint32_t stepId,
        vector<SegmentStepSupport>&
        );

    // Given a vector of SegmentStepSupport, for each OrientedReadId keep only the one
    // with the largest stepId.
    static void keepLast(vector<SegmentStepSupport>&);

    // Given a vector of SegmentStepSupport, for each OrientedReadId keep only the one
    // with the smallest stepId.
    static void keepFirst(vector<SegmentStepSupport>&);

    // Output a vector of SegmentStepSupport to a html table.
    static void writeHtml(ostream& html, const AssemblyGraph&, const vector<SegmentStepSupport>&);

    // Analyze SegmentStepSupport by comparing
    // the final representative region of e0
    // with the initial representative region of e1.
    static SegmentPairInformation analyzeSegmentPair(
        ostream& html,
        const AssemblyGraph&,
        edge_descriptor e0,
        edge_descriptor e1,
        uint32_t representativeRegionStepCount
        );

    // Estimate the offset between two segments using a SegmentStepSupport
    // on the first segment and a SegmentStepSupport on the second segment,
    // both for the same OrientedReadId.
    static int32_t estimateOffset(
        const AssemblyGraph&,
        SegmentStepSupport&,
        SegmentStepSupport&);
};
