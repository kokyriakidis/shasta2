// Shasta.
#include "findConvergingVertex.hpp"
using namespace shasta;

// Standard library.
#include <iostream.hpp>
#include <vector.hpp>



void shasta::testFindConvergingVertex()
{
    using Graph = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, uint64_t>;
    using vertex_descriptor = Graph::vertex_descriptor;
    Graph graph;

    // Create the vertices.
    vector<vertex_descriptor> v;
    for(uint64_t i=0; i<7; i++) {
        v.push_back(add_vertex(i, graph));
    }

    // Create the edges.
    add_edge(v[0], v[1], graph);
    add_edge(v[1], v[2], graph);
    add_edge(v[1], v[3], graph);
    add_edge(v[2], v[5], graph);
    add_edge(v[3], v[2], graph);
    add_edge(v[3], v[4], graph);
    add_edge(v[4], v[2], graph);
    add_edge(v[4], v[5], graph);
    add_edge(v[5], v[6], graph);
    // add_edge(v[5], v[1], graph);    // Loop edge

    for(const vertex_descriptor vA: v) {
        if(graph[vA] != 1) {
            // continue;
        }

        const vertex_descriptor vB = findConvergingVertex(graph, vA, 6);

        cout << "Starting at " << graph[vA] << ": ";
        if(vB == Graph::null_vertex()) {
            cout << "not found." << endl;
        } else {
            cout << "found " << graph[vB] << endl;
        }

    }

}
