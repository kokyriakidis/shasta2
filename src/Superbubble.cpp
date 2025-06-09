// Shasta.
#include "Superbubble.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>



Superbubble::Superbubble(
    const AssemblyGraph& assemblyGraph,
    vertex_descriptor source,
    vertex_descriptor target) :
    assemblyGraph(assemblyGraph),
    source(source),
    target(target)
{
    gatherInternalVertices();
    gatherInternalEdges();
}



void Superbubble::gatherInternalVertices()
{
    // Do a BFS starting at the source and stopping at the target.
    std::queue<vertex_descriptor> q;
    std::set<vertex_descriptor> internalVerticesSet;
    q.push(source);
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        // cout << "Dequeued " << assemblyGraph[v0].id << endl;
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = boost::target(e, assemblyGraph);
            if(v1 != target) {
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



void Superbubble::gatherInternalEdges()
{

    // Add the out-edges of the source.
    BGL_FORALL_OUTEDGES(source, e, assemblyGraph, AssemblyGraph) {
        internalEdges.push_back(e);
    }

    // Add the out-edges of the internal vertices.
    for(const vertex_descriptor v: internalVertices) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            internalEdges.push_back(e);
        }
    }

}

