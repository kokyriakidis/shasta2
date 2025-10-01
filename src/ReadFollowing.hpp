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
}

class shasta::ReadFollowing {
public:
    ReadFollowing(
        const AssemblyGraph&,
        uint64_t representativeRegionLength);

private:
    const AssemblyGraph& assemblyGraph;
    using AVertex = AssemblyGraph::vertex_descriptor;
    using AEdge = AssemblyGraph::edge_descriptor;
    using AssemblyGraphVertexPair = pair<AVertex, AVertex>;
    using AssemblyGraphEdgePair = pair<AEdge, AEdge>;
    using AVertexPair = AssemblyGraphVertexPair;
    using AEdgePair = AssemblyGraphEdgePair;

#if 0
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



    // A filtered version of the LineGraph that contains
    // only the vertices corresponding to some AssemblyGraph edges.
    // This is not permanently stored. It is only created as needed.
    class FilteringPredicate {
    public:
        FilteringPredicate(
            const LineGraph* lineGraph = 0,
            const std::set<AEdge>* assemblyGraphEdges = 0) :
            lineGraph(lineGraph),
            assemblyGraphEdges(assemblyGraphEdges)
        {}
        const LineGraph* lineGraph;
        const std::set<AEdge>* assemblyGraphEdges;

        bool operator()(const LineGraph::vertex_descriptor& lv) const
        {
            const AEdge ae = (*lineGraph)[lv];
            return assemblyGraphEdges->contains(ae);
        }

        bool operator()(const LineGraph::edge_descriptor& le) const
        {
            const LineGraph::vertex_descriptor lv0 = source(le, *lineGraph);
            const LineGraph::vertex_descriptor lv1 = target(le, *lineGraph);
            const AEdge ae0 = (*lineGraph)[lv0];
            const AEdge ae1 = (*lineGraph)[lv1];
            return assemblyGraphEdges->contains(ae0) and assemblyGraphEdges->contains(ae1);
        }
    };
    using FilteredLineGraph = boost::filtered_graph<LineGraph, FilteringPredicate, FilteringPredicate>;
#endif



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


    // The number of initial/final Appearances in each AssemblyGraph edge.
    std::map<AEdge, uint64_t> initialAppearancesCount;
    std::map<AEdge, uint64_t> finalAppearancesCount;
    void countAppearances();



    // Pairs of edges (e0, e1) such that one or more OrientedReadIds
    // appear in the final representative region of e0 and
    // in the initial representative region of e1,
    // with the journey position in e0 less than the journey position in e1.
    // The map is keyed by (e0, e1) and contains the number of
    // such OrientedReadIds.
    std::map<AssemblyGraphEdgePair, uint64_t> edgePairs;
    void findEdgePairs();

    // Some of the edgePairs are classified as "strong" - see findStrongEdgePairs for details.
    // The strong edgePairs define a graph in which vertex corresponds to
    // an AssemblyGraph edge, and a directed edge is added for each strong edgePair.
    using EdgePairsGraph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AEdge,
        uint64_t>;
    EdgePairsGraph edgePairsGraph;
    std::map<AEdge, EdgePairsGraph::vertex_descriptor> edgePairsVertexMap;
    void createEdgePairsGraph();
    void writeEdgePairsGraph();
};
