// Shasta.
#include "findConvergingVertex.hpp"
using namespace shasta;

// Standard library.
#include "iostream.hpp"
#include "vector.hpp"



void shasta::testFindConvergingVertex()
{
    using Graph = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, uint64_t>;
    using vertex_descriptor = Graph::vertex_descriptor;
    Graph graph;

    // Create the vertices.
    vector<vertex_descriptor> v;
    for(uint64_t i=0; i<10; i++) {
        v.push_back(add_vertex(i, graph));
    }

    // Create the edges.
    add_edge(v[0], v[1], graph);
    add_edge(v[1], v[2], graph);
    add_edge(v[2], v[3], graph);
    add_edge(v[2], v[4], graph);
    add_edge(v[3], v[5], graph);
    add_edge(v[4], v[6], graph);
    add_edge(v[5], v[7], graph);
    add_edge(v[6], v[7], graph);
    add_edge(v[7], v[8], graph);
    add_edge(v[8], v[9], graph);
    add_edge(v[6], v[3], graph);    // Loop edge

    for(const vertex_descriptor vA: v) {
        cout << "Starting at " << graph[vA] << endl;

        const vertex_descriptor vB = findConvergingVertexGeneral(graph, vA, 10);

        if(vB == Graph::null_vertex()) {
            cout << "Not found." << endl;
        } else {
            cout << "Found " << graph[vB] << endl;
        }

    }

}
