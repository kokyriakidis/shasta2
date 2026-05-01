#pragma once

// The AnchorSimilarityGraph is a directed graph in which
// each vertex represents an anchor. An edge anchorId0->anchorId1
// is created if the oriented read compositions of
// anchorId0 and anchorId1 are sufficiently similar.
// The connectivity of this graph is very high,
// so this works better when --min-anchor-coverage is
// low and most anchors correspond to a single copy of sequence.

// Shasta2.
#include "Anchor.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"



namespace shasta2 {

    class AnchorSimilarityGraphEdge;
    class AnchorSimilarityGraph;
    using AnchorSimilarityGraphBaseClass = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::bidirectionalS,
        boost::no_property,
        AnchorSimilarityGraphEdge>;

    class AnchorGraph;
}



class shasta2::AnchorSimilarityGraphEdge {
public:

    // The weight is computed by analyzing the AnchorPair for this edge.
    // It is low if the two anchors have very similar read composition.
    double weight = invalid<double>;

    uint64_t baseOffset = invalid<uint64_t>;

    AnchorSimilarityGraphEdge(
        double weight = invalid<double>,
        uint64_t baseOffset = invalid<uint64_t>) :
        weight(weight),
        baseOffset(baseOffset)
    {}

    template<class Archive> void serialize(
        Archive& ar,
        [[maybe_unused]] unsigned int version)
    {
        ar & weight;
        ar & baseOffset;
    }
};



class shasta2::AnchorSimilarityGraph :
    public AnchorSimilarityGraphBaseClass,
    public MappedMemoryOwner {
public:

    // Construct the AnchorSimilarityGraph from the completeAnchorGraph.
    // Only include edges with at least the specified minCoverage.
    AnchorSimilarityGraph(
        const Anchors&,
        const AnchorGraph& completeAnchorGraph);

    // Constructor from binary data.
    AnchorSimilarityGraph(const MappedMemoryOwner&, const string& name);

    // Compute a shortest path tree starting at the given AnchorId.
    void shortestPaths(AnchorId) const;

    // More efficient version with work areas.
    // The 3 work areas (predecessorMap, distanceMap, colorMap)
    // are vectors of size anchors.size().
    // On input, they must be set as follows.
    // On exit, they are returned in exactly the same state.
    // For performance, these conditions are not checked.
    // - predecessorMap[anchorId] == anchorId for all anchorIds.
    // - distanceMap[anchorId] == std::numeric_limits<double>::max().
    // - colorMap[anchorId] == boost::default_color_type::white_color.
    // This allows using dijkstra_shortest_paths_no_init.
    // On exit, the accessibleVertices contains the AnchorIds
    // that were seen and for which the predecessorMap and distanceMap
    // contains a valid value.
    void shortestPathsFast(
        AnchorId,
        vector<AnchorId>& predecessorMap,
        vector<double>& distanceMap,
        vector<boost::default_color_type>& colorMap,
        vector<AnchorId>& accessibleVertices
        ) const;
    void shortestPathsFast(AnchorId, const Anchors&) const;

    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorSimilarityGraphBaseClass>(*this);
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);

private:

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = 3;
    const double a = 3.;
    const double b = 10.;
    const double minLogP = 0.;


    void createVertices(const Anchors&);

    // Create all the edges.
    void createEdges(const Anchors&, const AnchorGraph& completeAnchorGraph);

    // Create the edges with source anchorIdA.
    void createEdges(
        const Anchors&,
        const AnchorGraph& completeAnchorGraph,
        AnchorId anchorIdA,
        vector<uint8_t>& color);
};

