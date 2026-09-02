#pragma once

#include <queue>
#include "vector.hpp"

namespace shasta2 {

    template<class Graph> void findDescendants(
        typename Graph::vertex_descriptor,
        uint64_t maxDistance,
        vector<typename Graph::vertex_descriptor>& descendants
        );

}



template<class Graph> void shasta2::findDescendants(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxDistance,
    vector<typename Graph::vertex_descriptor>& descendants
    )
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    descendants.clear();

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::map<vertex_descriptor, uint64_t> distanceMap;
    distanceMap.insert({vStart, 0});

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        const uint64_t distance0 = distanceMap.at(v0);
        const uint64_t distance1 = distance0 + 1;

        BGL_FORALL_OUTEDGES_T(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);
            if(not distanceMap.contains(v1)) {
                descendants.push_back(v1);
                distanceMap.insert({v1, distance1});
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }
    }

}
