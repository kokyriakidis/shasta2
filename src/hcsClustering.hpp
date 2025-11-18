#pragma once

// Suppress some warnings in Boost code.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

// HCS clustering for an undirected graph.
// https://en.wikipedia.org/wiki/HCS_clustering_algorithm
// The min-cut steps uses stoer_wagner_min_cut
// https://www.boost.org/doc/libs/latest/libs/graph/doc/stoer_wagner_min_cut.html

#include "SHASTA2_ASSERT.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#include <boost/graph/stoer_wagner_min_cut.hpp>
#pragma GCC diagnostic pop

#include <cstdint.hpp>
#include <map>
#include <tuple.hpp>
#include <vector.hpp>



namespace shasta {
    template<class InputGraph> void hcsClustering(
        const InputGraph&,
        vector< vector<typename InputGraph::vertex_descriptor> >& clusters
        );

    template <class InputGraph> class WorkGraph;
}



// We recursively use a WorkGraph in which each vertex contains
// a vertex_descriptor of the InputGraph.
template <class InputGraph> class shasta::WorkGraph :
    public boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    typename InputGraph::vertex_descriptor> {
public:

    // Create a WorkGraph from a connected component of the InputGraph.
    WorkGraph(
        const InputGraph& inputGraph,
        const vector<typename InputGraph::vertex_descriptor>& component)
    {
        // Add vertices.
        for(const typename InputGraph::vertex_descriptor v: component) {
            const typename WorkGraph::vertex_descriptor w = boost::add_vertex(v, *this);
            vertexMap.insert({v, w});
        }

        // Add edges.
        for(const typename InputGraph::vertex_descriptor v0: component) {
            const auto it0 = vertexMap.find(v0);
            SHASTA2_ASSERT(it0 != vertexMap.end());
            const typename WorkGraph::vertex_descriptor w0 = it0->second;
            typename InputGraph::out_edge_iterator it1, it1End;
            tie(it1, it1End) = boost::out_edges(v0, inputGraph);
            for(; it1!=it1End; ++it1) {
                const typename InputGraph::edge_descriptor e = *it1;
                const typename InputGraph::vertex_descriptor v1 = boost::target(e, inputGraph);
                const auto it1 = vertexMap.find(v1);
                SHASTA2_ASSERT(it1 != vertexMap.end());
                const typename WorkGraph::vertex_descriptor w1 = it1->second;

                // Don't add the edge twice.
                if(w0 < w1) {
                    boost::add_edge(w0, w1, *this);
                }
            }
        }
    }



    // Create the WorkGraph from one side of the min cut of another WorkGraph.
    WorkGraph(
        const WorkGraph& that,
        std::map<typename WorkGraph::vertex_descriptor, bool>& parityMap,
        bool side,
        const InputGraph& inputGraph)
    {
        // Add vertices.
        typename WorkGraph::vertex_iterator it, itEnd;
        tie(it, itEnd) = boost::vertices(that);
        for(; it!=itEnd; ++it) {
            const typename WorkGraph::vertex_descriptor wOld = *it;
            const auto it = parityMap.find(wOld);
            SHASTA2_ASSERT(it != parityMap.end());
            if(it->second == side) {
                const typename InputGraph::vertex_descriptor v = that[wOld];
                const typename WorkGraph::vertex_descriptor wNew = boost::add_vertex(v, *this);
                vertexMap.insert({v, wNew});
            }
        }

        // Add edges.
        tie(it, itEnd) = boost::vertices(*this);
        for(; it!=itEnd; ++it) {
            const typename WorkGraph::vertex_descriptor w0 = *it;
            const typename InputGraph::vertex_descriptor v = (*this)[w0];

            // Loop over the out-edges of v in the InputGraph.
            typename InputGraph::out_edge_iterator it1, it1End;
            tie(it1, it1End) = boost::out_edges(v, inputGraph);
            for(; it1!=it1End; ++it1) {
                const typename InputGraph::edge_descriptor e = *it1;
                const typename InputGraph::vertex_descriptor v1 = boost::target(e, inputGraph);
                const auto it2 = vertexMap.find(v1);
                if(it2 != vertexMap.end()) {
                    const typename WorkGraph::vertex_descriptor w1 = it2->second;
                    if(w0 < w1) {
                        boost::add_edge(w0, w1, *this);
                    }
                }
            }
        }

    }



    // Map vertex descriptors of the InputGraph to vertex descriptors of the WorkGraph.
    std::map<typename InputGraph::vertex_descriptor, typename WorkGraph::vertex_descriptor> vertexMap;



    void hcsClustering(
        const InputGraph& inputGraph,
        vector< vector<typename InputGraph::vertex_descriptor> >& clusters) const
    {
        // If just one vertex, generate a single cluster consisting of the entire WorkGraph.
        if(boost::num_vertices(*this) < 2) {
            addEntireGraphAsCluster(clusters);
            return;
        }

        // Compute the minimum cut.
        std::map<typename WorkGraph::vertex_descriptor, bool> parityMap;
        boost::stoer_wagner_min_cut(
            *this,
            boost::make_static_property_map<typename WorkGraph::edge_descriptor>(1),
            boost::parity_map(boost::make_assoc_property_map(parityMap)));

        // Count the edges on the min cut.
        uint64_t cutCount = 0;
        typename WorkGraph::edge_iterator it, itEnd;
        tie(it, itEnd) = boost::edges(*this);
        for(; it!=itEnd; ++it) {
            const typename WorkGraph::edge_descriptor e = *it;
            const typename WorkGraph::vertex_descriptor v0 = boost::source(e, *this);
            const typename WorkGraph::vertex_descriptor v1 = boost::target(e, *this);
            if(parityMap[v0] != parityMap[v1]) {
                 ++cutCount;
            }
        }


        // If this is a highly connected subgraph, add it as a cluster.
        if(2 * cutCount > num_vertices(*this)) {
            addEntireGraphAsCluster(clusters);
        } else {

            // This is not a highly connected subgraph.
            // Construct two new WorkGraphs, one for each side of the min cut.
            // Then do hcsClustering on both of them.
            for(uint64_t side=0; side<2; side++) {
                WorkGraph newWorkGraph(*this, parityMap, side==0, inputGraph);
                newWorkGraph.hcsClustering(inputGraph, clusters);
            }
        }

    }



    void addEntireGraphAsCluster(vector< vector<typename InputGraph::vertex_descriptor> >& clusters) const
    {
        clusters.emplace_back();
        vector<typename InputGraph::vertex_descriptor>& cluster = clusters.back();
        typename WorkGraph::vertex_iterator it, itEnd;
        tie(it, itEnd) = boost::vertices(*this);
        for(; it!=itEnd; ++it) {
            const uint64_t w = *it;
            const typename InputGraph::vertex_descriptor v = (*this)[w];
            cluster.push_back(v);
        }
    }

};



template<class InputGraph> void shasta::hcsClustering(
    const InputGraph& inputGraph,
    vector< vector<typename InputGraph::vertex_descriptor> >& clusters)
{
    using std::map;

    // Compute connected components of the InputGraph.

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<typename InputGraph::vertex_descriptor, uint64_t> vertexIndexMap;
    typename InputGraph::vertex_iterator it, itEnd;
    tie(it, itEnd) = boost::vertices(inputGraph);
    for(; it!=itEnd; ++it) {
        vertexIndexMap.insert({*it, vertexIndex++});
    }

    // Compute connected components.
    std::map<typename InputGraph::vertex_descriptor, uint64_t> componentMap;
    const uint64_t componentCount = boost::connected_components(
        inputGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));
    vector< vector<typename InputGraph::vertex_descriptor> > components(componentCount);
    for(const auto& p: componentMap) {
        const typename InputGraph::vertex_descriptor v = p.first;
        const uint64_t componentId = p.second;
        SHASTA2_ASSERT(componentId < componentCount);
        components[componentId].push_back(v);
    }



    // For each connected component construct a WorkGraph and
    // perform HCS clustering on that WorkGraph.
    clusters.clear();
    for(const vector<typename InputGraph::vertex_descriptor>& component: components) {
        const WorkGraph workGraph(inputGraph, component);
        SHASTA2_ASSERT(boost::num_vertices(workGraph) == component.size());
        workGraph.hcsClustering(inputGraph, clusters);
    }
}


#pragma GCC diagnostic pop
