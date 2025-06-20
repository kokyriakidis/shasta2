// Shasta.
#include "Superbubble.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>



Superbubble::Superbubble(
    const AssemblyGraph& assemblyGraph,
    vertex_descriptor sourceVertex,
    vertex_descriptor targetVertex) :
    assemblyGraph(assemblyGraph),
    sourceVertex(sourceVertex),
    targetVertex(targetVertex)
{
    gatherInternalVertices();
    gatherEdges();
}



void Superbubble::gatherInternalVertices()
{
    // Do a BFS starting at the sourceVertex and stopping at the targetVertex.
    std::queue<vertex_descriptor> q;
    std::set<vertex_descriptor> internalVerticesSet;
    q.push(sourceVertex);
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        // cout << "Dequeued " << assemblyGraph[v0].id << endl;
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(v1 != targetVertex) {
                if(not internalVerticesSet.contains(v1)) {
                    internalVerticesSet.insert(v1);
                    q.push(v1);
                    // cout << "Enqueued " << assemblyGraph[v1].id << endl;
                }
            }
        }
    }

    copy(internalVerticesSet.begin(), internalVerticesSet.end(), back_inserter(internalVertices));
}



// Gather the source, target, and internal edges.
// Internal edges include all source and target edges.
void Superbubble::gatherEdges()
{
    // The source edges are the out-edges of the source vertex.
    BGL_FORALL_OUTEDGES(sourceVertex, e, assemblyGraph, AssemblyGraph) {
        sourceEdges.push_back(e);
    }
    assemblyGraph.sortEdgeDescriptors(sourceEdges);

    // The target edges are the in-edges of the target vertex.
    BGL_FORALL_INEDGES(targetVertex, e, assemblyGraph, AssemblyGraph) {
        targetEdges.push_back(e);
    }
    assemblyGraph.sortEdgeDescriptors(targetEdges);

    // The internal edges are the out-edges of the source plus the
    // out-edges of all internal vertices.
    BGL_FORALL_OUTEDGES(sourceVertex, e, assemblyGraph, AssemblyGraph) {
        internalEdges.push_back(e);
    }
    for(const vertex_descriptor v: internalVertices) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            internalEdges.push_back(e);
        }
    }
    sort(internalEdges.begin(), internalEdges.end());
    assemblyGraph.sortEdgeDescriptors(internalEdges);

}

