// Shasta.
#include "Tangle.hpp"
#include "Anchor.hpp"
#include "TangleMatrix3.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <tuple.hpp>



// Constructor from a single AssemblyGraph vertex.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph3,
    vertex_descriptor v,
    double aDrift,
    double bDrift) :
    Tangle(assemblyGraph3,
        vector<vertex_descriptor>(1, v),
        aDrift, bDrift)
{
}



// Constructor from a single AssemblyGraph edge.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph3,
    edge_descriptor e,
    double aDrift,
    double bDrift) :
    Tangle(assemblyGraph3,
        vector<vertex_descriptor>({source(e, assemblyGraph3), target(e, assemblyGraph3)}),
        aDrift, bDrift)
{
}



// Constructor from a set of AssemblyGraph vertices.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph3,
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
        BGL_FORALL_INEDGES(v0, e, assemblyGraph3, AssemblyGraph) {
            const vertex_descriptor v1 = source(e, assemblyGraph3);
            if(not isTangleVertex(v1)) {
                entranceEdges.push_back(e);
            }
        }
    }

    // Find the exit edges.
    vector<edge_descriptor> exitEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph3, AssemblyGraph) {
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



void Tangle::connect(uint64_t iEntrance, uint64_t iExit) {
    SHASTA_ASSERT(iEntrance < tangleMatrix->entrances.size());
    SHASTA_ASSERT(iExit < tangleMatrix->exits.size());
    connectList.push_back({iEntrance, iExit});
}



void Tangle::detangle()
{
    // Make a copy of each entrance edge, with the target vertex replaced by a new vertex
    // with the same AnchorId.
    vector<vertex_descriptor> newEntranceVertices;
    for(const auto& entrance: tangleMatrix->entrances) {
        const edge_descriptor eOld = entrance.e;
        const vertex_descriptor v0Old = source(eOld, assemblyGraph3);
        AssemblyGraphEdge& edgeOld = assemblyGraph3[eOld];
        const AnchorId lastAnchorId = edgeOld.back().anchorPair.anchorIdB;

        // Create the new vertex.
        const vertex_descriptor v1 = add_vertex(
            AssemblyGraphVertex(lastAnchorId, assemblyGraph3.nextVertexId++), assemblyGraph3);
        newEntranceVertices.push_back(v1);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0Old, v1, AssemblyGraphEdge(edgeOld.id), assemblyGraph3);
        AssemblyGraphEdge& edgeNew = assemblyGraph3[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }



    // Make a copy of each exit edge, with the source vertex replaced by a new vertex
    // with the same AnchorId.
    vector<vertex_descriptor> newExitVertices;
    for(const auto& exit: tangleMatrix->exits) {
        const edge_descriptor eOld = exit.e;
        const vertex_descriptor v1Old = target(eOld, assemblyGraph3);
        AssemblyGraphEdge& edgeOld = assemblyGraph3[eOld];
        const AnchorId firstAnchorId = edgeOld.front().anchorPair.anchorIdA;

        // Create the new vertex.
        const vertex_descriptor v0 = add_vertex(
            AssemblyGraphVertex(firstAnchorId, assemblyGraph3.nextVertexId++), assemblyGraph3);
        newExitVertices.push_back(v0);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1Old, AssemblyGraphEdge(edgeOld.id), assemblyGraph3);
        AssemblyGraphEdge& edgeNew = assemblyGraph3[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }



    // Now connect these vertices as prescribed by the connectList.
    // The new edges get a single AnchorPair obtained from the
    // corresponding entry of the TangleMatrix.
    for(const auto& p: connectList) {
        const uint64_t iEntrance = p.first;
        const uint64_t iExit = p.second;

        const vertex_descriptor v0 = newEntranceVertices[iEntrance];
        const vertex_descriptor v1 = newExitVertices[iExit];

        const AnchorPair& anchorPair = (*tangleMatrix).tangleMatrix[iEntrance][iExit];
        const uint64_t offset = anchorPair.getAverageOffset(assemblyGraph3.anchors);

        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyGraphEdge(assemblyGraph3.nextEdgeId++), assemblyGraph3);
        AssemblyGraphEdge& edge = assemblyGraph3[e];

        edge.push_back(AssemblyGraphEdgeStep(anchorPair, offset));
    }

    // Now we can remove all the tangle vertices and their edges.
    // This also removes the old entrances and exits.
    removedVertices = tangleVertices;
    for(const vertex_descriptor v: tangleVertices) {
        clear_vertex(v, assemblyGraph3);
        remove_vertex(v, assemblyGraph3);
    }
}
