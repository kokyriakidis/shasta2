#pragma once

// A SegmentStepSupport object contains detailed information
// about one OrientedReadId that appears in an AssemblyGraphEdgeStep.

#include "AssemblyGraph.hpp"

namespace shasta {
    class SegmentStepSupport;
}


class shasta::SegmentStepSupport {
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

    // Get a vector of SegmentStepSupport for the last stepCount steps
    // of an AssemblyGraphEdge. This discard any previous contents of the vector.
    static void getFinal(
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

};
