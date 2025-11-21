#pragma once

// Given a directed graph, this does a forward BFS from a given start vertex
// up to a specified distance to construct a local subgraph.

// This is intrusive: the LocalSubgraph vertex is required
// to have the following two fields:
// - Graph::vertex_descriptor v
// - uint64_t distance
// The second overload also requires the LocalSubgraph edge
// to have a field Graph::edge_descriptor e.

#include "SHASTA2_ASSERT.hpp"

#include "cstdint.hpp"
#include <map>
#include <queue>
#include "vector.hpp"


namespace shasta {

    // Forward only, single start vertex.
    template<class Graph, class LocalSubgraph> LocalSubgraph createLocalSubgraph(
        const Graph&,
        typename Graph::vertex_descriptor startVertex,
        uint64_t maxDistance);

   // Forward and/or backward, multiple start vertices.
    template<class Graph, class LocalSubgraph> LocalSubgraph createLocalSubgraph(
        const Graph&,
        const vector<typename Graph::vertex_descriptor>& startVertices,
        bool forward,
        bool backward,
        uint64_t maxDistance);
}



// Forward only, single start vertex.
template<class Graph, class LocalSubgraph> LocalSubgraph shasta::createLocalSubgraph(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxDistance)
{
    using VertexDescriptor = typename Graph::vertex_descriptor;
    using EdgeDescriptor = typename Graph::edge_descriptor;
    using LocalVertexDescriptor = typename LocalSubgraph::vertex_descriptor;
    using OutEdgeIterator = typename Graph::out_edge_iterator;

    LocalSubgraph localSubgraph;

    // Create the local start vertex.
    const LocalVertexDescriptor lvStart = add_vertex(localSubgraph);
    auto& localVertex = localSubgraph[lvStart];
    localVertex.v = vStart;
    localVertex.distance = 0;

    // Map Graph vertices to LocalGraph vertices.
    std::map<VertexDescriptor, LocalVertexDescriptor> m;
    m.insert({vStart, lvStart});

    // Initialize the BFS.
    std::queue<VertexDescriptor> q;
    q.push(vStart);



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a global vertex.
        const VertexDescriptor v0 = q.front();
        q.pop();
        auto it0 = m.find(v0);
        SHASTA2_ASSERT(it0 != m.end());

        // Get the corresponding local vertex.
        const LocalVertexDescriptor& lv0 = it0->second;
        const auto& localVertex0 = localSubgraph[lv0];
        const uint64_t distance0 = localVertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over its out-edges.
        OutEdgeIterator itOut, itOutEnd;
        tie(itOut, itOutEnd) = out_edges(v0, graph);
        for(; itOut!=itOutEnd; ++itOut) {
            const EdgeDescriptor e = *itOut;
            const VertexDescriptor v1 = target(e, graph);

            // Check if we already encountered v1.
            auto it1 = m.find(v1);

            if(it1 == m.end()) {

                // We have not yet encountered v1.
                const LocalVertexDescriptor lv1 = add_vertex(localSubgraph);
                auto& localVertex1 = localSubgraph[lv1];
                localVertex1.v = v1;
                localVertex1.distance = distance1;
                m.insert({v1, lv1});
                add_edge(lv0, lv1, localSubgraph);
                if(distance1 < maxDistance) {
                   q.push(v1);
                }
            } else {

                // We already encountered v1. Just add this edge.
                LocalVertexDescriptor lv1 = it1->second;
                add_edge(lv0, lv1, localSubgraph);
            }

        }
    }

    return localSubgraph;
}



// Forward and/or backward, multiple start vertices.
// We want to make sure to get edges between vertices at maximum distance,
// so we first find all vertices, then all edges.
 template<class Graph, class LocalSubgraph> LocalSubgraph shasta::createLocalSubgraph(
     const Graph& graph,
     const vector<typename Graph::vertex_descriptor>& startVertices,
     bool forward,
     bool backward,
     uint64_t maxDistance)
 {
     using VertexDescriptor = typename Graph::vertex_descriptor;
     using EdgeDescriptor = typename Graph::edge_descriptor;
     using LocalVertexDescriptor = typename LocalSubgraph::vertex_descriptor;
     using LocalEdgeDescriptor = typename LocalSubgraph::edge_descriptor;
     using OutEdgeIterator = typename Graph::out_edge_iterator;
     using InEdgeIterator = typename Graph::in_edge_iterator;

     // Map Graph vertices to LocalGraph vertices.
     std::map<VertexDescriptor, LocalVertexDescriptor> m;

     // BFS queue.
     std::queue<VertexDescriptor> q;

     LocalSubgraph localSubgraph;

     // Create the local start vertices.
     for(const VertexDescriptor v: startVertices) {
         q.push(v);
         const LocalVertexDescriptor lv = add_vertex(localSubgraph);
         auto& localVertex = localSubgraph[lv];
         localVertex.v = v;
         localVertex.distance = 0;
         m.insert({v, lv});
     }

     if(maxDistance == 0) {
         return localSubgraph;
     }



     // Main BFS loop to find the vertices.
     while(not q.empty()) {

         // Dequeue a global vertex.
         const VertexDescriptor v0 = q.front();
         q.pop();
         auto it0 = m.find(v0);
         SHASTA2_ASSERT(it0 != m.end());

         // Get the corresponding local vertex.
         const LocalVertexDescriptor& lv0 = it0->second;
         const auto& localVertex0 = localSubgraph[lv0];
         const uint64_t distance0 = localVertex0.distance;
         const uint64_t distance1 = distance0 + 1;

         // Loop over its out-edges.
         if(forward) {
             OutEdgeIterator itOut, itOutEnd;
             tie(itOut, itOutEnd) = out_edges(v0, graph);
             for(; itOut!=itOutEnd; ++itOut) {
                 const EdgeDescriptor e = *itOut;
                 const VertexDescriptor v1 = target(e, graph);

                 // Check if we already encountered v1.
                 auto it1 = m.find(v1);

                 if(it1 == m.end()) {

                     // We have not yet encountered v1.
                     const LocalVertexDescriptor lv1 = add_vertex(localSubgraph);
                     auto& localVertex1 = localSubgraph[lv1];
                     localVertex1.v = v1;
                     localVertex1.distance = distance1;
                     m.insert({v1, lv1});
                     if(distance1 < maxDistance) {
                        q.push(v1);
                     }
                 }
             }
         }

         // Loop over its in-edges.
         if(backward) {
             InEdgeIterator itOut, itOutEnd;
             tie(itOut, itOutEnd) = in_edges(v0, graph);
             for(; itOut!=itOutEnd; ++itOut) {
                 const EdgeDescriptor e = *itOut;
                 const VertexDescriptor v1 = source(e, graph);

                 // Check if we already encountered v1.
                 auto it1 = m.find(v1);

                 if(it1 == m.end()) {

                     // We have not yet encountered v1.
                     const LocalVertexDescriptor lv1 = add_vertex(localSubgraph);
                     auto& localVertex1 = localSubgraph[lv1];
                     localVertex1.v = v1;
                     localVertex1.distance = distance1;
                     m.insert({v1, lv1});
                     if(distance1 < maxDistance) {
                        q.push(v1);
                     }
                 }
             }
         }
     }



     // Now add all the edges. This way we are guaranteed to get edges
     // between vertices at maximum distance.
     BGL_FORALL_VERTICES_T(lv0, localSubgraph, LocalSubgraph) {
         const VertexDescriptor v0 = localSubgraph[lv0].v;
         BGL_FORALL_OUTEDGES_T(v0, e, graph, Graph) {
             const VertexDescriptor v1 = target(e, graph);
             const auto it1 = m.find(v1);
             if(it1 != m.end()) {
                 const LocalVertexDescriptor lv1 = it1->second;
                 LocalEdgeDescriptor le;
                 bool edgeWasAdded = false;
                 tie(le, edgeWasAdded) = add_edge(lv0, lv1, localSubgraph);
                 SHASTA2_ASSERT(edgeWasAdded);
                 localSubgraph[le].e = e;
             }
         }
     }




     return localSubgraph;
 }
