// Shasta.
#include "Tangle2.hpp"
#include "Anchor.hpp"
#include "TangleMatrix2.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// Constructor from a single AssemblyGraph vertex.
Tangle2::Tangle2(
    AssemblyGraph2& assemblyGraph2,
    vertex_descriptor v,
    double aDrift,
    double bDrift) :
    Tangle2(assemblyGraph2, vector<vertex_descriptor>(1, v), aDrift, bDrift)
{
}



// Constructor from a set of AssemblyGraph2 vertices.
Tangle2::Tangle2(
    AssemblyGraph2& assemblyGraph2,
    const vector<vertex_descriptor>& tangleVerticesArgument,
    double aDrift,
    double bDrift) :
    assemblyGraph2(assemblyGraph2),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangleVertices so we can do binary searches in it using isTangleVertex.
    sort(tangleVertices.begin(), tangleVertices.end());

    // Find the entrance vertices.
    vector<vertex_descriptor> entranceVertices;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph2, AssemblyGraph2) {
            const vertex_descriptor v1 = source(e, assemblyGraph2);
            if(not isTangleVertex(v1)) {
                entranceVertices.push_back(v1);
            }
        }
    }

    // Find the exit vertices.
    vector<vertex_descriptor> exitVertices;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph2, AssemblyGraph2) {
            const vertex_descriptor v1 = target(e, assemblyGraph2);
            if(not isTangleVertex(v1)) {
                exitVertices.push_back(v1);
            }
        }
    }

    // Now we can create the TangleMatrix2.
    tangleMatrix = make_shared<TangleMatrix2>(
        assemblyGraph2, entranceVertices, exitVertices, aDrift, bDrift);
}



void Tangle2::connect(uint64_t iEntrance, uint64_t iExit) {
    SHASTA_ASSERT(iEntrance < tangleMatrix->entrances.size());
    SHASTA_ASSERT(iExit < tangleMatrix->exits.size());
    connectList.push_back({iEntrance, iExit});
}


