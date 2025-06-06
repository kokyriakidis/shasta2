#pragma once

// Shasta.
#include "invalid.hpp"
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/rational.hpp>


// Standard library.
#include <cstdint.hpp>
#include <map>
#include <queue>
#include <utility.hpp>
#include <vector.hpp>

// REMOVE THIS WHEN DONE DEBUGGING.
// #include <iostream.hpp>


// Give a vertex vA of a directed graph, this finds the first vertex,
// (first in topological ordering) such that:
// - Distance(vA, vB) <= maxDistance.
// - All paths that start at vA go through vB.
// If such a vertex is not found or a loop was encountered,
// this returns Graph::null_vertex().

namespace shasta {
   template<class Graph> typename Graph::vertex_descriptor
       findConvergingVertex(const Graph&, typename Graph::vertex_descriptor, uint64_t maxDistance);

   // Test function for the above.
   void testFindConvergingVertex();
}



template<class Graph> typename Graph::vertex_descriptor shasta::findConvergingVertex(
    const Graph& graph,
    typename Graph::vertex_descriptor vA,
    uint64_t maxDistance)
{
    using VertexDescriptor = typename Graph::vertex_descriptor;
    using Rational = boost::rational<uint64_t>;

    // const bool debug = false;

    // The LocalGraph is constructed using a BFS starting at vA up to distance maxDistance.
    // Each vertex contains a vertex_descriptor of the Graph.
    class LocalVertex {
    public:
        VertexDescriptor v = Graph::null_vertex();
        uint64_t distance = invalid<uint64_t>;
        Rational flow = 0;
        LocalVertex(VertexDescriptor v, uint64_t distance) : v(v), distance(distance) {}
        LocalVertex() {}
    };

    using LocalGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, LocalVertex>;
    using LocalVertexDescriptor = typename LocalGraph::vertex_descriptor;
    LocalGraph localGraph;

    // Initialize the BFS.
    std::queue<VertexDescriptor> q;
    q.push(vA);

    // Map Graph vertices to LocalGraph vertices.
    std::map<VertexDescriptor, LocalVertexDescriptor> m;
    const LocalVertexDescriptor lvA = add_vertex(LocalVertex(vA, 0), localGraph);
    m.insert({vA, lvA});



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const VertexDescriptor v0 = q.front();
        q.pop();
        auto it0 = m.find(v0);
        SHASTA_ASSERT(it0 != m.end());
        const LocalVertexDescriptor& lv0 = it0->second;
        const LocalVertex& localVertex0 = localGraph[lv0];
        const uint64_t distance0 = localVertex0.distance;
        const uint64_t distance1 = distance0 + 1;

#if 0
        if(debug) {
            cout << "Dequeued " << graph[v0] << " at distance " << distance0 << endl;
        }
#endif

        // Loop over its out-edges.
        BGL_FORALL_OUTEDGES_T(v0, e, graph, Graph) {
            const VertexDescriptor v1 = target(e, graph);

#if 0
            if(debug) {
                cout << "Encountered " << graph[v1] << endl;
            }
#endif

            // Check if we already encountered v1.
            auto it1 = m.find(v1);
            if(it1 == m.end()) {
#if 0
                if(debug) {
                    cout << "Adding " << graph[v1] << " to local graph at distance " << distance1 << endl;
                }
#endif
                const LocalVertexDescriptor lv1 = add_vertex(LocalVertex(v1, distance1), localGraph);
                m.insert({v1, lv1});
                add_edge(lv0, lv1, localGraph);
#if 0
                if(debug) {
                    cout << "Adding " << graph[v0] << "->" << graph[v1] << endl;
                }
#endif
                if(distance1 < maxDistance) {
                   q.push(v1);
                }
            } else {
                LocalVertexDescriptor lv1 = it1->second;
                add_edge(lv0, lv1, localGraph);
#if 0
                if(debug) {
                    cout << "Adding " << graph[v0] << "->" << graph[v1] << endl;
                }
#endif
            }

        }
    }

#if 0
    if(debug) {
        cout << "The BFS found " << num_vertices(localGraph) <<
            " vertices and " << num_edges(localGraph) << " edges." << endl;
    }
#endif

    if(num_vertices(localGraph) == 1) {
        return Graph::null_vertex();
    }

#if 0
    // Write out the LocalGraph in Graphviz format.
    if(debug) {
        cout << "digraph LocalGraph {\n";
        BGL_FORALL_VERTICES_T(lv, localGraph, LocalGraph) {
            cout << graph[localGraph[lv].v] << ";\n";
        }
        BGL_FORALL_EDGES_T(le, localGraph, LocalGraph) {
            const LocalVertexDescriptor lv0 = source(le, localGraph);
            const LocalVertexDescriptor lv1 = target(le, localGraph);
            cout << graph[localGraph[lv0].v] << "->";
            cout << graph[localGraph[lv1].v] << ";\n";
        }
        cout << "}\n";
    }
#endif

    // Do a topological sort of the LocalGraph.
    vector<LocalVertexDescriptor> topologicalOrder;
    try {
        boost::topological_sort(localGraph, back_inserter(topologicalOrder));
    } catch(const boost::not_a_dag&) {
#if 0
        if(debug) {
            cout << "Found a loop." << endl;
        }
#endif
        return Graph::null_vertex();
    }
    reverse(topologicalOrder.begin(), topologicalOrder.end());

#if 0
    if(debug) {
        cout << "Topological order:";
        for(const LocalVertexDescriptor lv: topologicalOrder) {
            cout << " " << graph[localGraph[lv].v];
        }
        cout << endl;
    }
#endif

    // Sanity check on the topological ordering.
    SHASTA_ASSERT(topologicalOrder.size() == num_vertices(localGraph));
    SHASTA_ASSERT(localGraph[topologicalOrder.front()].v == vA);


    // Compute flow in topological order.
    localGraph[topologicalOrder.front()].flow = 1;
    for(uint64_t i=1; i<topologicalOrder.size(); i++) {
        const LocalVertexDescriptor lv1 = topologicalOrder[i];
        LocalVertex& localVertex1 = localGraph[lv1];
        localVertex1.flow = 0;

        BGL_FORALL_INEDGES_T(lv1, e, localGraph, LocalGraph){
            const LocalVertexDescriptor lv0 = source(e, localGraph);
            localVertex1.flow += localGraph[lv0].flow / out_degree(lv0, localGraph);
        }
    }

#if 0
    if(debug) {
        BGL_FORALL_VERTICES_T(lv, localGraph, LocalGraph) {
            const LocalVertex& localVertex = localGraph[lv];
            cout << graph[localVertex.v] << " has flow " << localVertex.flow << endl;
        }
    }
#endif

    // Find the first vertex, in topological order, with flow equal to 1.
    for(uint64_t i=1; i<topologicalOrder.size(); i++) {
        const LocalVertexDescriptor lv = topologicalOrder[i];
        LocalVertex& localVertex = localGraph[lv];
        if(localVertex.flow == 1) {
            return localVertex.v;
        }
    }

    return Graph::null_vertex();
}
