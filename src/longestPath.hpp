#pragma once

// Find the longest path in a directed acyclic graph
// stored as a boost::adjacency_list.
// If the graph has cycles, this throws boost::not_a_dag.
// This uses the algorithm described here:
// https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include <map>
#include <vector.hpp>


namespace shasta {
    template<class Graph> void longestPath(
        const Graph &graph,
        vector<typename Graph::edge_descriptor>& longestPath);
}



template<class Graph> void shasta::longestPath(
    const Graph &graph,
    vector<typename Graph::edge_descriptor>& longestPath)
{
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        vertexIndexMap.insert({v, vertexIndex++});
    }

    // Use boost topological_sort to get a vector of vertex descriptors
    // in topological order. The output from the boost call is in
    // reverse topological order.
    // If the graph has cycles, this throws boost::not_a_dag.
    vector<vertex_descriptor> sortedVertices;
    topological_sort(
        graph,
        back_inserter(sortedVertices),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));
    std::reverse(sortedVertices.begin(), sortedVertices.end());


    // Map to contain the length of the longest path ending at each vertex.
    std::map<vertex_descriptor, uint64_t> lengthMap;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        lengthMap.insert(make_pair(v, 0));
    }

    // Compute the maximum length of a path ending at each vertex.
    for(const vertex_descriptor v: sortedVertices) {
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            maximumLength = max(maximumLength, lengthMap[source(e, graph)]);
        }
        lengthMap[v] = maximumLength + 1;
    }

    // Find the vertex with the longest length.
    // This will be the end of the longest path.
    vertex_descriptor v = Graph::null_vertex();
    uint64_t maximumLength = 0;
    for(const auto& p: lengthMap) {
        if(p.second > maximumLength) {
            v = p.first;
            maximumLength = p.second;
        }
    }


    // Construct the path, moving backward from here.
    longestPath.clear();
    while(true) {
        vertex_descriptor vPrevious = Graph::null_vertex();
        edge_descriptor ePrevious;
        uint64_t maximumLength = 0;
        BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const uint64_t length = lengthMap[v0];
            if(length > maximumLength) {
                vPrevious = v0;
                ePrevious = e;
                maximumLength = length;
            }
        }
        if(vPrevious == Graph::null_vertex()) {
            break;
        }
        v = vPrevious;
        longestPath.push_back(ePrevious);

    }
    std::reverse(longestPath.begin(), longestPath.end());
}




