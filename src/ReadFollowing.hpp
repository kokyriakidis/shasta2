#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>



namespace shasta {
    class ReadFollowing;
}

class shasta::ReadFollowing {
public:
    ReadFollowing(
        const AssemblyGraph&,
        uint64_t representativeRegionLength,
        uint64_t minCoverage);
    void followForward(AssemblyGraph::edge_descriptor) const;

private:
    const AssemblyGraph& assemblyGraph;
    using AVertex = AssemblyGraph::vertex_descriptor;
    using AEdge = AssemblyGraph::edge_descriptor;
    using AssemblyGraphVertexPair = pair<AVertex, AVertex>;
    using AssemblyGraphEdgePair = pair<AEdge, AEdge>;


    // A "line graph" representation of the AssemblyGraph.
    // Names or vertices and edges of the line graph are prefixed with "l".
    // Each vertex of the LineGraph corresponds to an edge of the AssemblyGraph.
    // Given LineGraph vertices lv0 and lv1 corresponding to AssemblyGraph
    // edges e0 and e1, a LineGraph edge lv0->lv1 is created if
    // target(e0, assemblyGraph) == source(e1, assemblyGraph).
    // The LineGraph is similar to the Bandage representation of the AssemblyGraph.
    using LineGraph = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::bidirectionalS,
        AEdge>;
    LineGraph lineGraph;
    std::map<AEdge, LineGraph::vertex_descriptor> lineGraphVertexMap;
    void createLineGraph();


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
    void findAppearances(uint64_t representativeRegionLength);



    // Pairs of edges (e0, e1) such that one or more OrientedReadIds
    // appear in the final representative region of e0 and
    // in the initial representative region of e1.
    // with the journey position in e0 less than the journey position in e1.
    // The map is keyed by (e0, e1) and contains the number of
    // such OrientedReadIds.
    std::map<AssemblyGraphEdgePair, uint64_t> edgePairs;
    void findEdgePairs();

    // The edge pairs with coverage >= minCoverage define a directed graph
    // in which each vertex corresponds to an AssemblyGraph edge.
    using EdgePairsGraph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AEdge,
        uint64_t>;
    EdgePairsGraph edgePairsGraph;
    std::map<AEdge, EdgePairsGraph::vertex_descriptor> edgePairsVertexMap;
    void createEdgePairsGraph(uint64_t minCoverage);
    void writeEdgePairsGraph();
};
