// Shasta.
#include "Tangle.hpp"
using namespace shasta;

// Boost libraries,
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <tuple.hpp>


// Constructor from a set of AssemblyGraph vertices.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVertices) :
    assemblyGraph(assemblyGraph),
    sorter(assemblyGraph),
    tangleVertices(tangleVertices)
{
    construct();
}



// Constructor from a single AssemblyGraph vertex.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    vertex_descriptor v) :
    assemblyGraph(assemblyGraph),
    sorter(assemblyGraph)
{
    tangleVertices.push_back(v);
    construct();
}



// Constructor from a single AssemblyGraph edge.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    edge_descriptor e) :
    assemblyGraph(assemblyGraph),
    sorter(assemblyGraph)
{
    const vertex_descriptor v0 = source(e, assemblyGraph);
    const vertex_descriptor v1 = target(e, assemblyGraph);

    tangleVertices.push_back(v0);
    if(v1 != v0) {
        tangleVertices.push_back(v1);
    }

    construct();
}



void Tangle::construct()
{
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



void Tangle::detangle()
{

    // For each Entrance we create a new edge with the last AnchorId removed
    // and store the last vertex of that new edge. If the Entrance already has
    // a just two AnchorIds, we just store the source vertex of that Entrance.
    // These are the vertices that will be connected to similar vertices created below
    // for each exit.
    vector<vertex_descriptor> entranceNewVertices(tangleMatrix.entrances.size());
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        const edge_descriptor eOld = tangleMatrix.entrances[iEntrance].e;
        const AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
        const vertex_descriptor v0 = source(eOld, assemblyGraph);
        SHASTA_ASSERT(oldEdge.size() > 1);
        if(oldEdge.size() == 2) {
            entranceNewVertices[iEntrance] = v0;
        } else {
            const vertex_descriptor v1 = add_vertex(assemblyGraph);
            entranceNewVertices[iEntrance] = v1;
            assemblyGraph[v1].anchorId = oldEdge.secondToLast();
            const edge_descriptor eNew = assemblyGraph.addEdge(v0, v1);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            copy(oldEdge.begin(), oldEdge.end()-1, back_inserter(newEdge));
        }
    }



    // Same as above, but for the Exits. Here we remove the first AnchorId.
    vector<vertex_descriptor> exitNewVertices(tangleMatrix.exits.size());
    for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
        const edge_descriptor eOld = tangleMatrix.exits[iExit].e;
        const AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
        const vertex_descriptor v0 = target(eOld, assemblyGraph);
        if(oldEdge.size() == 2) {
            exitNewVertices[iExit] = v0;
        } else {
            const vertex_descriptor v1 = add_vertex(assemblyGraph);
            exitNewVertices[iExit] = v1;
            assemblyGraph[v1].anchorId = oldEdge.second();
            const edge_descriptor eNew = assemblyGraph.addEdge(v1, v0);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            copy(oldEdge.begin()+1, oldEdge.end(), back_inserter(newEdge));
        }
    }



    // Now for each pair to be connected we connect the above vertices.
    // The connecting edge has just two AnchorIds.
    for(const auto& p: connectList) {
        const uint64_t iEntrance = p.first;
        const uint64_t iExit = p.second;

        // Create the connecting edge.
        const vertex_descriptor v0 = entranceNewVertices[iEntrance];
        const vertex_descriptor v1 = exitNewVertices[iExit];
        const edge_descriptor eNew = assemblyGraph.addEdge(v0, v1);
        AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
        newEdge.push_back(assemblyGraph[v0].anchorId);
        newEdge.push_back(assemblyGraph[v1].anchorId);
    }

    // Store the edges that will be removed.
    removedEdges.clear();
    // Edges internal to the Tangle and exit edges.
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            removedEdges.push_back(e);
        }
    }
    // Entrance edges.
    for(const TangleMatrix::Entrance& entrance: tangleMatrix.entrances) {
        removedEdges.push_back(entrance.e);
    }


    // Now we can remove the tangle vertices. This also removes all edges internal
    // to the tangles plus the entrances and exits.
    for(const vertex_descriptor v: tangleVertices) {
        boost::clear_vertex(v, assemblyGraph);
        boost::remove_vertex(v, assemblyGraph);
    }
}

