// Shasta.
#include "Tangle1.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// Constructor from a set of AssemblyGraph vertices.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVerticesArgument) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangle vertices by id.
    std::ranges::sort(tangleVertices, assemblyGraph.orderById);

    // Find the edges.
    findEntrances();
    findExits();
    findTangleEdges();

    // Compute the tangle matrix.
    ostream html(0);
    tangleMatrix = make_shared<TangleMatrix1>(assemblyGraph, entrances, exits, html);
}



// Constructor from a single vertex.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    vertex_descriptor v) :
    Tangle1(
        assemblyGraph,
        vector<vertex_descriptor>(1, v))
{}



// Constructor from an edge.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    edge_descriptor e) :
    Tangle1(
        assemblyGraph,
        vector<vertex_descriptor>({source(e, assemblyGraph), target(e, assemblyGraph)}))
{}



void Tangle1::findEntrances()
{
    entrances.clear();
    for(const vertex_descriptor v1: tangleVertices) {
        BGL_FORALL_INEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            if(not isTangleVertex(v0)) {
                entrances.push_back(e);
            }
        }
    }
}



void Tangle1::findExits()
{
    entrances.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                exits.push_back(e);
            }
        }
    }
}



void Tangle1::findTangleEdges()
{
    tangleEdges.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(isTangleVertex(v1)) {
                tangleEdges.push_back(e);
            }
        }
    }
}



bool Tangle1::isTangleVertex(vertex_descriptor v) const
{
    return std::ranges::binary_search(tangleVertices, v, assemblyGraph.orderById);
}



void Tangle1::connect(uint64_t /* iEntrance */, uint64_t /* iExit */)
{
    SHASTA_ASSERT(0);
}



void Tangle1::detangle()
{
    SHASTA_ASSERT(0);
}

