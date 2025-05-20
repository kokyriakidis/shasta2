// Shasta.
#include "Tangle3.hpp"
#include "Anchor.hpp"
#include "TangleMatrix3.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>


// Constructor from a single AssemblyGraph edge.
Tangle3::Tangle3(
    AssemblyGraph3& assemblyGraph3,
    edge_descriptor e,
    double aDrift,
    double bDrift) :
    Tangle3(assemblyGraph3,
        vector<vertex_descriptor>({source(e, assemblyGraph3), target(e, assemblyGraph3)}),
        aDrift, bDrift)
{
}



// Constructor from a set of AssemblyGraph3 vertices.
Tangle3::Tangle3(
    AssemblyGraph3& assemblyGraph3,
    const vector<vertex_descriptor>& tangleVerticesArgument,
    double aDrift,
    double bDrift) :
    assemblyGraph3(assemblyGraph3),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangleVertices so we can do binary searches in it using isTangleVertex.
    sort(tangleVertices.begin(), tangleVertices.end());

    // Find the entrance edges.
    vector<edge_descriptor> entranceEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph3, AssemblyGraph3) {
            const vertex_descriptor v1 = source(e, assemblyGraph3);
            if(not isTangleVertex(v1)) {
                entranceEdges.push_back(e);
            }
        }
    }

    // Find the exit edges.
    vector<edge_descriptor> exitEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph3, AssemblyGraph3) {
            const vertex_descriptor v1 = target(e, assemblyGraph3);
            if(not isTangleVertex(v1)) {
                exitEdges.push_back(e);
            }
        }
    }

    // Now we can create the TangleMatrix3.
    tangleMatrix = make_shared<TangleMatrix3>(
        assemblyGraph3, entranceEdges, exitEdges, aDrift, bDrift);
}
