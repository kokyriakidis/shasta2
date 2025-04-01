// Shasta.
#include "AssemblyGraphPostprocessor.hpp"
#include "Detangler.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const Anchors& anchors,
    const string& assemblyStage) :
    AssemblyGraph(anchors, assemblyStage)
{
    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AnchorId anchorId = assemblyGraph[v].anchorId;
        SHASTA_ASSERT(not vertexMap.contains(anchorId));
        vertexMap.insert(make_pair(anchorId, v));
    }

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const uint64_t id = assemblyGraph[e].id;
        SHASTA_ASSERT(not edgeMap.contains(id));
        edgeMap.insert(make_pair(id, e));
    }
}



// Get the vertex_descriptor corresponding to an AnchorId.
AssemblyGraph::vertex_descriptor AssemblyGraphPostprocessor::getVertexDescriptor(AnchorId anchorId) const
{
    auto it = vertexMap.find(anchorId);
    SHASTA_ASSERT(it != vertexMap.end());
    return it->second;
}



// Get the vertex corresponding to a vertex_descriptor.
AssemblyGraphVertex& AssemblyGraphPostprocessor::getVertex(vertex_descriptor v)
{
    return (*this)[v];
}



// Get the edge_descriptor corresponding to a segmentId.
AssemblyGraph::edge_descriptor AssemblyGraphPostprocessor::getEdgeDescriptor(uint64_t segmentId) const
{
    auto it = edgeMap.find(segmentId);
    SHASTA_ASSERT(it != edgeMap.end());
    return it->second;
}



// Get the edge corresponding to an edge_descriptor.
AssemblyGraphEdge& AssemblyGraphPostprocessor::getEdge(edge_descriptor e)
{
    return (*this)[e];
}



// These are needed to simplify the Python API.
void AssemblyGraphPostprocessor::detangleVertices(TrivialDetangler& detangler) {
    AssemblyGraph::detangleVertices(detangler);
}
void AssemblyGraphPostprocessor::detangleEdges(TrivialDetangler& detangler) {
    AssemblyGraph::detangleEdges(detangler);
}

