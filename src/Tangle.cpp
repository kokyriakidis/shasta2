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
    tangleVertices.push_back(source(e, assemblyGraph));
    tangleVertices.push_back(target(e, assemblyGraph));
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

    // For each Entrance we create a new edge with the last AssembyGraphStep removed
    // and store the last vertex of that new edge. If the Entrance already has
    // a single AssemblyGraphStep, we just store the source vertex of that Entrance.
    // These are the vertices that will be connected to similar vertices created below
    // for each exit.
    vector<vertex_descriptor> entranceNewVertices(tangleMatrix.entrances.size());
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        const edge_descriptor eOld = tangleMatrix.entrances[iEntrance].e;
        const AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
        const vertex_descriptor v0 = source(eOld, assemblyGraph);
        if(oldEdge.size() == 1) {
            entranceNewVertices[iEntrance] = v0;
        } else {
            const vertex_descriptor v1 = add_vertex(assemblyGraph);
            entranceNewVertices[iEntrance] = v1;
            assemblyGraph[v1].anchorId = oldEdge.back().anchorPair.anchorIdA;
            const edge_descriptor eNew = assemblyGraph.addEdge(v0, v1);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            copy(oldEdge.begin(), oldEdge.end()-1, back_inserter(newEdge));
            const AssemblyGraphStep& lastOldStep = oldEdge.back();
            newEdge.averageOffset = oldEdge.averageOffset - lastOldStep.averageOffset;
            newEdge.minOffset = oldEdge.minOffset - lastOldStep.minOffset;
            newEdge.maxOffset = oldEdge.maxOffset - lastOldStep.maxOffset;
        }
    }



    // Same as above, but for the Exits. Here we remove the first AssemblyGraphStep.
    vector<vertex_descriptor> exitNewVertices(tangleMatrix.exits.size());
    for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
        const edge_descriptor eOld = tangleMatrix.exits[iExit].e;
        const AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
        const vertex_descriptor v0 = target(eOld, assemblyGraph);
        if(oldEdge.size() == 1) {
            exitNewVertices[iExit] = v0;
        } else {
            const vertex_descriptor v1 = add_vertex(assemblyGraph);
            exitNewVertices[iExit] = v1;
            assemblyGraph[v1].anchorId = oldEdge.front().anchorPair.anchorIdB;
            const edge_descriptor eNew = assemblyGraph.addEdge(v1, v0);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            copy(oldEdge.begin()+1, oldEdge.end(), back_inserter(newEdge));
            const AssemblyGraphStep& firstOldStep = oldEdge.front();
            newEdge.averageOffset = oldEdge.averageOffset - firstOldStep.averageOffset;
            newEdge.minOffset = oldEdge.minOffset - firstOldStep.minOffset;
            newEdge.maxOffset = oldEdge.maxOffset - firstOldStep.maxOffset;
        }
    }



    // Now for each pair to be connected we connect the above vertices.
    // The connecting edge has a single AssemblyGraphStep obtaining
    // by "mixing" the last AssemblyGraphStep of the Entrance with the
    // first AssemblyGraphStep of the Exit being connected.
    for(const auto& p: connectList) {
        const uint64_t iEntrance = p.first;
        const uint64_t iExit = p.second;

        // Create the connecting edge.
        const vertex_descriptor v0 = entranceNewVertices[iEntrance];
        const vertex_descriptor v1 = exitNewVertices[iExit];
        const edge_descriptor eNew = assemblyGraph.addEdge(v0, v1);
        AssemblyGraphEdge& newEdge = assemblyGraph[eNew];

        // Fill in the single AssemblyGraphStep of the new edge.
        newEdge.resize(1);
        AssemblyGraphStep& newStep = newEdge.front();

        const edge_descriptor eOldEntrance = tangleMatrix.entrances[iEntrance].e;
        const edge_descriptor eOldExit = tangleMatrix.exits[iExit].e;
        const AssemblyGraphStep& oldStepEntrance = assemblyGraph[eOldEntrance].back();
        const AssemblyGraphStep& oldStepExit = assemblyGraph[eOldExit].front();

        newStep.anchorPair.anchorIdA = oldStepEntrance.anchorPair.anchorIdA;
        newStep.anchorPair.anchorIdB = oldStepExit.anchorPair.anchorIdB;
        set_intersection(
            oldStepEntrance.anchorPair.orientedReadIds.begin(), oldStepEntrance.anchorPair.orientedReadIds.end(),
            oldStepExit.anchorPair.orientedReadIds.begin(), oldStepExit.anchorPair.orientedReadIds.end(),
            back_inserter(newStep.anchorPair.orientedReadIds));
        SHASTA_ASSERT(not newStep.anchorPair.orientedReadIds.empty());

        // Compute the offsets for the new AssemblyGraphStep and the new AssemblyGraphEdge.
        // They are the same because the edge consists of a single step.
        newStep.anchorPair.getOffsetStatistics(
            assemblyGraph.anchors,
            newStep.averageOffset,
            newStep.minOffset,
            newStep.maxOffset);
        newEdge.averageOffset = newStep.averageOffset;
        newEdge.minOffset = newStep.minOffset;
        newEdge.maxOffset = newStep.maxOffset;
    }

    // Now we can remove the tangle vertices. This also removes all edges internal
    // to the tangles plus the entrances and exits.
    for(const vertex_descriptor v: tangleVertices) {
        boost::clear_vertex(v, assemblyGraph);
        boost::remove_vertex(v, assemblyGraph);
    }
}
