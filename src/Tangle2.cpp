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




void Tangle2::detangle()
{

    // For each entrance/exit pair in the connectList, create a new
    // vertex joining the entrance and exit.
    for(const auto& p: connectList) {
        const uint64_t iEntrance = p.first;
        const uint64_t iExit = p.second;

        const auto& entrance = tangleMatrix->entrances[iEntrance];
        const auto& exit = tangleMatrix->exits[iExit];

        const vertex_descriptor vEntrance = entrance.v;
        const vertex_descriptor vExit = exit.v;

        const AnchorPair& bridgeAnchorPair = tangleMatrix->tangleMatrix[iEntrance][iExit];
        const uint64_t bridgeOffset = bridgeAnchorPair.getAverageOffset(assemblyGraph2.anchors);

        // Create the new vertex.
        const vertex_descriptor vNew = add_vertex(AssemblyGraph2Vertex(assemblyGraph2.nextVertexId++), assemblyGraph2);
        AssemblyGraph2Vertex& vertexNew = assemblyGraph2[vNew];
        vertexNew.emplace_back(bridgeAnchorPair, bridgeOffset);

        // Connect it to the entrance and exit.
        add_edge(vEntrance, vNew, assemblyGraph2);
        add_edge(vNew, vExit, assemblyGraph2);

    }

    // Now we can remove all the tangle vertices.
    // The entrance and exit vertices are not removed.
    for(const vertex_descriptor v: tangleVertices) {
        removedVertices.push_back(v);
        clear_vertex(v, assemblyGraph2);
        remove_vertex(v, assemblyGraph2);
    }

}
