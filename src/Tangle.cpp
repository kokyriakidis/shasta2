// Shasta.
#include "Tangle.hpp"
using namespace shasta;

// Boost libraries,
#include <boost/graph/iteration_macros.hpp>



// Constructor from a set of AssemblyGraph vertices.
Tangle::Tangle(
    const AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& unsortedTangleVertices) :
    sorter(assemblyGraph)
{
    construct(assemblyGraph, unsortedTangleVertices);
}



// Constructor from a single AssemblyGraph vertex.
Tangle::Tangle(
    const AssemblyGraph& assemblyGraph,
    vertex_descriptor v) :
    sorter(assemblyGraph)
{
    vector<vertex_descriptor> unsortedTangleVertices = {v};
    construct(assemblyGraph, unsortedTangleVertices);
}



// Constructor from a single AssemblyGraph edge.
Tangle::Tangle(
    const AssemblyGraph& assemblyGraph,
    edge_descriptor e) :
    sorter(assemblyGraph)
{
    vector<vertex_descriptor> unsortedTangleVertices =
        {source(e, assemblyGraph), target(e, assemblyGraph)};
    construct(assemblyGraph, unsortedTangleVertices);
}



void Tangle::construct(
    const AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& unsortedTangleVertices)
{
    tangleVertices = unsortedTangleVertices;
    sort(tangleVertices.begin(), tangleVertices.end(), sorter);

    // Find the entrance edges.
    vector<edge_descriptor> entranceEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = source(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                entranceEdges.push_back(e);
            }
        }
    }

    // Find the exit edges.
    vector<edge_descriptor>exitEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                exitEdges.push_back(e);
            }
        }
    }

    // Now we can compute the TangleMatrix.
    tangleMatrix.construct(assemblyGraph, entranceEdges, exitEdges);
}
