#pragma once

// Given a directed graph, this computes a "partial" transitive reduction
// as follows:
// - It computes the CondensedGraph, which is acyclic.
// - it calls shasta2:;transitiveReduction on the CondensedGraph.
// - It then transports the result of the transitive reduction to the
//   original Graph.

#include "CondensedGraph.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <map>
#include <set>
#include "vector.hpp"


namespace shasta2 {
    template<class Graph> void partialTransitiveReduction(Graph&);
}


template<class Graph> void shasta2::partialTransitiveReduction(Graph& graph)
{
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // Create the CondensedGraph.
    class CondensedGraphVertex {
    public:
        vector<vertex_descriptor> vertices;
    };
    class CondensedGraphEdge {
    public:
        edge_descriptor e = {0, 0, 0};
    };
    using CondensedGraph = adjacency_list<listS, vecS, bidirectionalS,
        CondensedGraphVertex, CondensedGraphEdge>;
    std::map<vertex_descriptor, typename CondensedGraph::vertex_descriptor> vertexMap;
    CondensedGraph condensedGraph = createCondensedGraph<Graph, CondensedGraph>(graph, vertexMap);

    // Do transitive reduction on the CondensedGraph.
    transitiveReduction(condensedGraph);

    // Gather the edges that survived the transitive reduction.
    std::set<edge_descriptor> survivingEdges;
    BGL_FORALL_EDGES_T(ce, condensedGraph, CondensedGraph) {
        survivingEdges.insert(condensedGraph[ce].e);
    }

    // Remove all remaining edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        if(not survivingEdges.contains(e)) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}
