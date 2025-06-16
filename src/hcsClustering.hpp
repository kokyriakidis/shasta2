#pragma once

// HCS clustering for an undirected graph.
// https://en.wikipedia.org/wiki/HCS_clustering_algorithm
// The min-cut steps uses stoer_wagner_min_cut
// https://www.boost.org/doc/libs/latest/libs/graph/doc/stoer_wagner_min_cut.html

#include "SHASTA_ASSERT.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>

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
            SHASTA_ASSERT(it0 != vertexMap.end());
            const typename WorkGraph::vertex_descriptor w0 = it0->second;
            typename InputGraph::out_edge_iterator it1, it1End;
            tie(it1, it1End) = boost::out_edges(v0, inputGraph);
            for(; it1!=it1End; ++it1) {
                const typename InputGraph::edge_descriptor e = *it1;
                const typename InputGraph::vertex_descriptor v1 = boost::target(e, inputGraph);
                const auto it1 = vertexMap.find(v1);
                SHASTA_ASSERT(it1 != vertexMap.end());
                const typename WorkGraph::vertex_descriptor w1 = it1->second;
                boost::add_edge(w0, w1, *this);
            }
        }
    }

    // Map vertex descriptors of the InputGraph to vertex descriptors of the WorkGraph.
    std::map<typename InputGraph::vertex_descriptor, typename WorkGraph::vertex_descriptor> vertexMap;

    void hcsClustering(vector< vector<typename WorkGraph::vertex_descriptor> >& componentClusters) const
    {
        // For now return a single cluster consisting of the entire WorkGraph.
        componentClusters.clear();
        componentClusters.emplace_back();
        auto& componentCluster = componentClusters.front();
        typename WorkGraph::vertex_iterator it, itEnd;
        tie(it, itEnd) = boost::vertices(*this);
        for(; it!=itEnd; ++it) {
            const uint64_t w = *it;
            const typename InputGraph::vertex_descriptor v = (*this)[w];
            componentCluster.push_back(v);
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
    std::map<typename InputGraph::vertex_descriptor, uint64_t> vertexMap;
    typename InputGraph::vertex_iterator it, itEnd;
    tie(it, itEnd) = boost::vertices(inputGraph);
    for(; it!=itEnd; ++it) {
        vertexMap.insert({*it, vertexIndex++});
    }

    // Compute connected components.
    std::map<typename InputGraph::vertex_descriptor, uint64_t> componentMap;
    const uint64_t componentCount = boost::connected_components(
        inputGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));
    vector< vector<typename InputGraph::vertex_descriptor> > components(componentCount);
    for(const auto& p: componentMap) {
        const typename InputGraph::vertex_descriptor v = p.first;
        const uint64_t componentId = p.second;
        SHASTA_ASSERT(componentId < componentCount);
        components[componentId].push_back(v);
    }



    // For each connected component construct a WorkGraph and
    // perform HCS clustering on that WorkGraph.
    clusters.clear();
    for(const vector<typename InputGraph::vertex_descriptor>& component: components) {
        const WorkGraph workGraph(inputGraph, component);
        vector< vector<typename InputGraph::vertex_descriptor> > componentClusters;
        workGraph.hcsClustering(componentClusters);

        for(const auto& componentCluster: componentClusters) {
            clusters.push_back(componentCluster);
        }

    }
}
