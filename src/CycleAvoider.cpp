#include "CycleAvoider.hpp"
using namespace shasta;

#include "iostream.hpp"



void shasta::testCycleAvoider()
{
    // Define the graph.
    class Vertex : public CycleAvoiderVertex {};
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex>;
    using edge_descriptor = Graph::edge_descriptor;
    Graph graph;

    // Create the CyckeAvoider.
    CycleAvoider<Graph> cycleAvoider(graph);

    // Add the vertices.
    const uint64_t n = 6;
    for(uint64_t i=0; i<n; i++) {
        cycleAvoider.addVertex();
    }

    // Try adding edges. The ones that generate cycles will not be added.
    vector< pair<uint64_t, uint64_t> > edges = {
        {3, 4},
        {4, 6},
        {6, 3},
        {1, 2}
    };

    for(const auto& edge: edges) {
        const uint64_t v0 = edge.first;
        const uint64_t v1 = edge.second;
        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = cycleAvoider.addEdge(v0, v1);

        cout << "Edge " << v0 << " " << v1 << " was ";
        if(not edgeWasAdded) {
            cout << "not ";
        }
        cout << "added." << endl;
    }
}
