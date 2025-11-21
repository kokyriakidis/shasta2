#pragma once

/*******************************************************************************

A class that can be used to avoid creating cycles in a Boost directed graph.

The Graph MUST not call boost::add_vertex and boost::add_edge
to create vertices and edges. Instead, it must call
CycleAvoider::addVertex and CycleAvoider::addEdge.

The Graph MUST NOT remove vertices or edges.

The vertex type MUST be derived from CycleAvoiderVertex
and must be default constructible.

CycleAvoider::addEdge returns true if the edge was added, in which case
no cycles were created. If adding the edge would create cycles,
or boost::add_edge fails, addEdge returns false and the edge is not added.

This uses the PK algorithm for incremental topological sort described in:
David. J. Pearce and Paul H. J. Kelly,
A Dynamic Topological Sort Algorithm for Directed Acyclic graphs,
ACM Journal of Experimental Algorithmics 11, 1-24 (2006).

Note that a C++ implementation of this algorithm
is available on GitHub:
https://github.com/DavePearce/DynamicTopologicalSort
However, this code was written from scratch using the algorithm
description in the paper.

********************************************************************************/

// Shasta2.
#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "cstdint.hpp"
#include <stack>
#include "utility.hpp"
#include "vector.hpp"



namespace shasta2 {
    template<class Graph> class CycleAvoider;
    class CycleAvoiderVertex;

    void testCycleAvoider();
}



class shasta2::CycleAvoiderVertex {
public:
    uint64_t rank = invalid<uint64_t>;
    uint64_t color = 0;
};



template<class Graph> class shasta2::CycleAvoider {
public:

    CycleAvoider(Graph& graph) :
        graph(graph)
    {

    }

    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    uint64_t nextRank = 0;
    vertex_descriptor addVertex()
    {
        vertex_descriptor v = boost::add_vertex(graph);
        graph[v].rank = nextRank++;
        return v;
    }



    pair<edge_descriptor, bool> addEdge(vertex_descriptor vX, vertex_descriptor vY)
    {
        const CycleAvoiderVertex& vertexX = graph[vX];
        const CycleAvoiderVertex& vertexY = graph[vY];

        // If this edge is consistent with the existing ranks, add it.
        // The topological sort is not affected.
        if(vertexX.rank < vertexY.rank) {
            return boost::add_edge(vX, vY, graph);
        }

        deltaRanks.clear();

        // Use a forward DFS starting at vY, and limited to the affected region,
        // to compute deltaF (see definition 2.6 in the paper).
        deltaF.clear();
        SHASTA2_ASSERT(vertexStack.empty());
        vertexStack.push(vY);
        deltaF.push_back(make_pair(vertexY.rank, vY));
        deltaRanks.push_back(vertexY.rank);
        graph[vY].color = 1;
        while(!vertexStack.empty()) {
            const vertex_descriptor v0 = vertexStack.top();
            vertexStack.pop();
            BGL_FORALL_OUTEDGES_T(v0, e01, graph, Graph) {
                const vertex_descriptor v1 = target(e01, graph);
                auto& vertex1 = graph[v1];
                if(vertex1.rank > vertexX.rank) {
                    // Outside the affected region.
                    continue;
                }
                if(vertex1.color == 1) {
                    // We already have this vertex.
                    continue;
                }
                vertex1.color = 1;
                vertexStack.push(v1);
                deltaF.push_back(make_pair(vertex1.rank, v1));
                deltaRanks.push_back(vertex1.rank);
            }
        }
        sort(deltaF.begin(), deltaF.end());
        for(const auto& p: deltaF) {
            graph[p.second].color = 0;
        }



        // Use a backward DFS starting at vX, and limited to the affected region,
        // to compute deltaB (see definition 2.6 in the paper).
        deltaB.clear();
        SHASTA2_ASSERT(vertexStack.empty());
        vertexStack.push(vX);
        deltaB.push_back(make_pair(vertexX.rank, vX));
        deltaRanks.push_back(vertexX.rank);
        graph[vX].color = 1;
        while(!vertexStack.empty()) {
            const vertex_descriptor v0 = vertexStack.top();
            vertexStack.pop();
            BGL_FORALL_INEDGES_T(v0, e01, graph, Graph) {
                const vertex_descriptor v1 = source(e01, graph);
                auto& vertex1 = graph[v1];
                if(vertex1.rank < vertexY.rank) {
                    // Outside the affected region.
                    continue;
                }
                if(vertex1.color == 1) {
                    // We already have this vertex.
                    continue;
                }
                vertex1.color = 1;
                vertexStack.push(v1);
                deltaB.push_back(make_pair(vertex1.rank, v1));
                deltaRanks.push_back(vertex1.rank);
            }
        }
        sort(deltaB.begin(), deltaB.end());
        for(const auto& p: deltaB) {
            graph[p.second].color = 0;
        }


        // Sort the ranks. If we find any duplicates, deltaF and deltaB
        // intersect, which means this edge introduces a cycle.
        // So we don't add it.
        sort(deltaRanks.begin(), deltaRanks.end());
        for(size_t i=1; i<deltaRanks.size(); i++) {
            if(deltaRanks[i-1] == deltaRanks[i]) {
                return make_pair(edge_descriptor(), false);
            }
        }

        // Adding this edge will not create a cycle.
        // redistribute the deltaRanks to vertices in deltaB and deltaF.
        size_t i = 0;
        for(const auto& p: deltaB) {
           graph[p.second].rank = deltaRanks[i++];
        }
        for(const auto& p: deltaF) {
           graph[p.second].rank = deltaRanks[i++];
        }

        // Now we can create the edge.
        return boost::add_edge(vX, vY, graph);
    }



private:
    Graph& graph;

    // Vectors to hold deltaF and deltaB
    // (see definition 2.6 in the paper).
    // Stored as pairs (rank, vertex descriptor).
    vector< pair<size_t, vertex_descriptor> > deltaF;
    vector< pair<size_t, vertex_descriptor> > deltaB;

    // Vector to hold the ranks of vertices in deltaF and deltaB.
    vector<size_t> deltaRanks;

    // Stack used for DFS's.
    std::stack<vertex_descriptor> vertexStack;
};
