#pragma once

// In the AnchorGraph, each vertex corresponds to an AnchorId
// and each edge corresponds to an AnchorPair.
// It uses boost::vecS as its second template argument,
// and as a result its vertex descriptors are AnchorIds.
// The AnchorGraph is not stored. It is used to create
// the initial Assembly graph and then discarded.

// Shasta.
#include "AnchorPair.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {

    class AnchorGraph;
    using AnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        boost::no_property,
        AnchorPair>;

    class Anchors;
    class Journeys;
}



class shasta::AnchorGraph : public AnchorGraphBaseClass {
public:
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        uint64_t minEdgeCoverage);

    // Constructor that splits edges that have an AnchorPair
    // with inconsistent offsets.
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        uint64_t minEdgeCoverage,
        double aDrift,
        double bDrift);

    // Compute the edge journeys.
    // The edge journey of an OrientedReadId is the sequence of
    // AnchorGraph edges visited by the OrientedReadId.
    // An OrientedReadId visits an AnchorGraph edges if it
    // appears in the AnchorPairt for the edge.
    // Edge journeys are indexed by OrientedReadId::getValue().
    void computeEdgeJourneys(
        const Anchors&,
        vector< vector<edge_descriptor> >&
    ) const;
};

