#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

// Standard library.
#include <map>



namespace shasta {

    class ReadFollowing;
    class ReadFollowingVertex;
    class ReadFollowingEdge;

    using ReadFollowingBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        ReadFollowingVertex,
        ReadFollowingEdge>;
}



class shasta::ReadFollowingVertex {
public:
    AssemblyGraph::edge_descriptor ae;
    uint64_t length;
    bool isLong;
    ReadFollowingVertex(
        const AssemblyGraph&,
        AssemblyGraph::edge_descriptor,
        uint64_t longLengthThreshold);
};



class shasta::ReadFollowingEdge {
public:
    uint64_t coverage = invalid<uint64_t>;
    ReadFollowingEdge(uint64_t coverage) : coverage(coverage) {}
};



class shasta::ReadFollowing : public ReadFollowingBaseClass {
public:
    ReadFollowing(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    using AVertex = AssemblyGraph::vertex_descriptor;
    using AEdge = AssemblyGraph::edge_descriptor;
    using AssemblyGraphVertexPair = pair<AVertex, AVertex>;
    using AssemblyGraphEdgePair = pair<AEdge, AEdge>;
    using AVertexPair = AssemblyGraphVertexPair;
    using AEdgePair = AssemblyGraphEdgePair;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t representativeRegionLength = 10;
    const uint64_t minCoverage = 3;
    const double minCoverageFraction = 0.8;
    const uint64_t maxAppearanceCount = 25;
    const double minJaccard = 0.05;
    const uint64_t shortLengthThreshold = 1000;
    const uint64_t longLengthThreshold = 500000;

    // Find appearances of OrientedReadIds in the initial/final
    // representative regions of each edge.
    // Store them by OrientedReadId, indexed by OrientedReadId::getValue().
    // For the initial representative region, store the last appearance
    // (in journey order) of each OrientedReadId.
    // For the final representative region, store the first appearance
    // (in journey order) of each OrientedReadId.
    class Appearance {
    public:
        AEdge e;
        uint32_t positionInJourney;
        Appearance(AEdge e, uint32_t positionInJourney) : e(e), positionInJourney(positionInJourney) {}
    };
    vector< vector<Appearance> > initialAppearances;
    vector< vector<Appearance> > finalAppearances;
    void findAppearances();


    // The number of initial/final Appearances in each AssemblyGraph edge.
    std::map<AEdge, uint64_t> initialAppearancesCount;
    std::map<AEdge, uint64_t> finalAppearancesCount;
    void countAppearances();
    uint64_t getInitialAppearancesCount(AEdge) const;
    uint64_t getFinalAppearancesCount(AEdge) const;



    // Pairs of edges (e0, e1) such that one or more OrientedReadIds
    // appear in the final representative region of e0 and
    // in the initial representative region of e1,
    // with the journey position in e0 less than the journey position in e1.
    // The map is keyed by (e0, e1) and contains the number of
    // such OrientedReadIds.
    std::map<AssemblyGraphEdgePair, uint64_t> edgePairs;
    void findEdgePairs();

    // Some of the edgePairs are classified as "strong" - see createGraph for details.
    // The strong edgePairs define a graph in which vertex corresponds to
    // an AssemblyGraph edge, and a directed edge is added for each strong edgePair.
    std::map<AEdge, vertex_descriptor> vertexMap;
    void createGraph();
    void writeGraph() const;

    // Jaccard similarity for an EdgePairsGraph edge.
    // Computed using the finalAppearancesCount of the source vertex
    // and the initialAppearancesCount of the target vertex.
    double jaccard(AEdge, AEdge, uint64_t coverage) const;
    double jaccard(edge_descriptor) const;

    // In the graph, find a path that starts at a given AEdge
    // and moves forward/backward. At each step we choose the child vertex
    // corresponding to best Jaccard similarity.
public:
    void findPath(AEdge, uint64_t direction) const;
    void findForwardPath(AEdge) const;
    void findBackwardPath(AEdge) const;
};
