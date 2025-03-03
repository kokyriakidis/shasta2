#ifndef SHASTA_FIND_LINEAR_CHAINS_HPP
#define SHASTA_FIND_LINEAR_CHAINS_HPP

// Find linear chains in a directed graph

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <list>
#include "vector.hpp"

namespace shasta {

    // Find linear chains of edges (paths).
    template<class Graph> void findLinearChains(
        const Graph&,
        uint64_t minimumLength,
        vector< std::list<typename Graph::edge_descriptor> >&);
    template<class Graph> void findLinearChains(
        const Graph&,
        uint64_t minimumLength,
        vector< vector<typename Graph::edge_descriptor> >&);


    // Find linear chains of vertices.
    template<class Graph> void findLinearVertexChains(
        const Graph&,
        vector< std::list<typename Graph::vertex_descriptor> >&);
    template<class Graph> void findLinearVertexChains(
        const Graph&,
        vector< vector<typename Graph::vertex_descriptor> >&);
}



// Version that uses vectors.
template<class Graph> void shasta::findLinearChains(
    const Graph& graph,
    uint64_t minimumLength,
    vector< vector<typename Graph::edge_descriptor> >& chainVectors)
{
    using edge_descriptor = typename Graph::edge_descriptor;

    // Call the version that uses lists.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(graph, minimumLength, chains);

    // Copy the lists to vectors.
    chainVectors.clear();
    for(const auto& chain: chains) {
        chainVectors.push_back(vector<edge_descriptor>(chain.begin(), chain.end()));
    }
}



// Find linear chains of edges (paths).
template<class Graph> inline void shasta::findLinearChains(
    const Graph& graph,
    uint64_t minimumLength,
    vector< std::list<typename Graph::edge_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // The edges we have already encountered.
    std::set<edge_descriptor> edgesFound;


    chains.clear();

    // Consider all possible start edges for the chain.
    BGL_FORALL_EDGES_T(eStart, graph, Graph) {

        // If we already assigned this edge to a chain, skip it.
        if(edgesFound.find(eStart) != edgesFound.end()) {
            continue;
        }

        // Add a new chain consisting of the start edge.
        chains.resize(chains.size() + 1);
        std::list<edge_descriptor>& chain = chains.back();
        chain.push_back(eStart);
        edgesFound.insert(eStart);

        // Extend forward.
        bool isCircular = false;
        edge_descriptor e = eStart;
        while(true) {
            const vertex_descriptor v = target(e, graph);
            if(in_degree(v, graph) != 1) {
                break;
            }
            if(out_degree(v, graph) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, eNext, graph, Graph) {
                e = eNext;
                break;
            }
            if(e == eStart) {
                isCircular = true;
                break;
            }
            chain.push_back(e);
            SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
            edgesFound.insert(e);
        }


        // Extend backward.
        if(not isCircular) {
            edge_descriptor e = eStart;
            while(true) {
                const vertex_descriptor v = source(e, graph);
                if(in_degree(v, graph) != 1) {
                    break;
                }
                if(out_degree(v, graph) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, ePrevious, graph, Graph) {
                    e = ePrevious;
                    break;
                }
                if(e == eStart) {
                    isCircular = true;
                    break;
                }
                chain.push_front(e);
                SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
                edgesFound.insert(e);
            }
        }

        // If the chain is too short, get rid of it.
        if(chain.size() < minimumLength) {
            chains.resize(chains.size() - 1);
        }

    }

    // Check that all edges were found.
    // Just using num_edges does not work if the graph is a filtered_graph.
    // SHASTA_ASSERT(edgesFound.size() == num_edges(graph));
    uint64_t edgeCount = 0;
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        ++edgeCount;
    }
    SHASTA_ASSERT(edgesFound.size() == edgeCount);
}



// Find linear chains of vertices.
template<class Graph> void shasta::findLinearVertexChains(
    const Graph& graph,
    vector< std::list<typename Graph::vertex_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // The vertices we have already encountered.
    std::set<vertex_descriptor> verticesFound;

    chains.clear();

    // Consider all possible start vertices for the chain.
    BGL_FORALL_VERTICES_T(vStart, graph, Graph) {

        // If we already assigned this vertex to a chain, skip it.
        if(verticesFound.find(vStart) != verticesFound.end()) {
            continue;
        }

        // Add a new chain consisting of the start vertex.
        chains.resize(chains.size() + 1);
        std::list<vertex_descriptor>& chain = chains.back();
        chain.push_back(vStart);
        verticesFound.insert(vStart);

        // Extend forward.
        bool isCircular = false;
        vertex_descriptor v = vStart;
        while(true) {
            if(out_degree(v, graph) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
                v = target(e, graph);
                break;
            }
            if(v == vStart) {
                isCircular = true;
                break;
            }
            if(in_degree(v, graph) != 1) {
                break;
            }
            chain.push_back(v);
            SHASTA_ASSERT(verticesFound.find(v) == verticesFound.end());
            verticesFound.insert(v);
        }

        // Extend backward.
        if(not isCircular) {
            vertex_descriptor v = vStart;
            while(true) {
                if(in_degree(v, graph) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
                    v = source(e, graph);
                    break;
                }
                if(out_degree(v, graph) != 1) {
                    break;
                }
                chain.push_front(v);
                SHASTA_ASSERT(verticesFound.find(v) == verticesFound.end());
                verticesFound.insert(v);
            }
        }
    }



    // Check that all vertices were found.
    // Just using num_vertices does not work if the graph is a filtered_graph.
    // SHASTA_ASSERT(verticesFound.size() == num_vertices(graph));
    uint64_t vertexCount = 0;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        if(v != Graph::null_vertex()) { // Just to avoid compiler warning.
            ++vertexCount;
        }
    }
    SHASTA_ASSERT(verticesFound.size() == vertexCount);
}



template<class Graph> void shasta::findLinearVertexChains(
    const Graph& graph,
    vector< vector<typename Graph::vertex_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Find the chains.
    vector< std::list<typename Graph::vertex_descriptor> > chainLists;
    findLinearVertexChains(graph, chainLists);

    // Copy lists to vectors.
    chains.clear();
    chains.reserve(chainLists.size());
    for(const auto& chain: chainLists) {
        chains.push_back(vector<vertex_descriptor>());
        copy(chain.begin(), chain.end(), back_inserter(chains.back()));
    }
}


#endif

