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
