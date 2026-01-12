// Shasta.
#include "Tangle.hpp"
#include "Anchor.hpp"
#include "TangleMatrix.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "tuple.hpp"



// Constructor from a single AssemblyGraph vertex.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    vertex_descriptor v,
    uint64_t maxTrim,
    double aDrift,
    double bDrift) :
    Tangle(assemblyGraph,
        vector<vertex_descriptor>(1, v),
        maxTrim, aDrift, bDrift)
{
}



// Constructor from a single AssemblyGraph edge.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    edge_descriptor e,
    uint64_t maxTrim,
    double aDrift,
    double bDrift) :
    Tangle(assemblyGraph,
        vector<vertex_descriptor>({source(e, assemblyGraph), target(e, assemblyGraph)}),
        maxTrim, aDrift, bDrift)
{
}



// Constructor from a set of AssemblyGraph vertices.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVerticesArgument,
    uint64_t maxTrim,
    double aDrift,
    double bDrift) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangleVertices so we can do binary searches in it using isTangleVertex.
    sort(tangleVertices.begin(), tangleVertices.end(), assemblyGraph.orderById);

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
    sort(entranceEdges.begin(), entranceEdges.end(), assemblyGraph.orderById);

    // Find the exit edges.
    vector<edge_descriptor> exitEdges;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                exitEdges.push_back(e);
            }
        }
    }
    sort(exitEdges.begin(), exitEdges.end(), assemblyGraph.orderById);

    // Now we can create the TangleMatrix.
    tangleMatrix = make_shared<TangleMatrix>(
        assemblyGraph, entranceEdges, exitEdges, maxTrim, aDrift, bDrift);
}



void Tangle::connect(uint64_t iEntrance, uint64_t iExit) {
    SHASTA2_ASSERT(iEntrance < tangleMatrix->entrances.size());
    SHASTA2_ASSERT(iExit < tangleMatrix->exits.size());
    connectList.push_back({iEntrance, iExit});
}



void Tangle::detangle()
{
    const bool debug = false;

    // Make a copy of each entrance edge, with the target vertex replaced by a new vertex
    // with the same AnchorId. Remove trim steps at its end.
    vector<vertex_descriptor> newEntranceVertices;
    for(const auto& entrance: tangleMatrix->entrances) {
        const edge_descriptor eOld = entrance.e;
        const vertex_descriptor v0Old = source(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId lastAnchorId = edgeOld[edgeOld.size() - 1 - entrance.trim].anchorPair.anchorIdB;

        // Create the new vertex.
        const vertex_descriptor v1 = add_vertex(
            AssemblyGraphVertex(lastAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newEntranceVertices.push_back(v1);

        // Create the new edge, with the same steps and id as the old one except for the last trim removed.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0Old, v1, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
        edgeNew.resize(edgeNew.size() - entrance.trim);
    }



    // Make a copy of each exit edge, with the source vertex replaced by a new vertex
    // with the same AnchorId. Remove trim steps at its beginning.
    vector<vertex_descriptor> newExitVertices;
    for(const auto& exit: tangleMatrix->exits) {
        const edge_descriptor eOld = exit.e;
        const vertex_descriptor v1Old = target(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId firstAnchorId = edgeOld[exit.trim].anchorPair.anchorIdA;

        // Create the new vertex.
        const vertex_descriptor v0 = add_vertex(
            AssemblyGraphVertex(firstAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newExitVertices.push_back(v0);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1Old, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.wasAssembled = edgeOld.wasAssembled;
        copy(edgeOld.begin() + exit.trim, edgeOld.end(), back_inserter(edgeNew));
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
        SHASTA2_ASSERT(anchorPair.orientedReadIds.size() > 0);
        const uint64_t offset = anchorPair.getAverageOffset(assemblyGraph.anchors);

        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];

        edge.push_back(AssemblyGraphEdgeStep(anchorPair, offset));

        if(debug) {
            cout << "Created edge " << edge.id << " to connect entrances " << iEntrance << " " << iExit << endl;
        }
    }

    // Now we can remove all the tangle vertices and their edges.
    // This also removes the old entrances and exits.
    removedVertices = tangleVertices;
    for(const vertex_descriptor v: tangleVertices) {
        if(debug) {
            cout << "Removing for detangling:";
            BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
        }
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }
}



void Tangle::computeExtendedTangleMatrix(
    vector< vector<double> >& extendedTangleMatrix) const
{
    vector<edge_descriptor> entranceEdges;
    for(const auto& entrance: tangleMatrix->entrances) {
        entranceEdges.push_back(entrance.e);
    }
    vector<edge_descriptor> exitEdges;
    for(const auto& exit: tangleMatrix->exits) {
        exitEdges.push_back(exit.e);
    }

    ostream html(0);
    assemblyGraph.computeExtendedTangleMatrix(entranceEdges, exitEdges, extendedTangleMatrix, html);
}
