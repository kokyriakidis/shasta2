#pragma once

#include "AssemblyGraph.hpp"


namespace shasta {
    class AssemblyGraphPostprocessor;

    class PermutationDetangler;
    class TrivialDetangler;
}



// AssemblyGraph functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class shasta::AssemblyGraphPostprocessor : public AssemblyGraph {
public:
    AssemblyGraphPostprocessor(
        const AssemblerOptions&,
        const Anchors&,
        const string& assemblyStage);



    // Some accessors used in the Python API.

    // Get the vertex_descriptor corresponding to an AnchorId.
    vertex_descriptor getVertexDescriptor(AnchorId) const;

    // Get the vertex corresponding to a vertex_descriptor.
    AssemblyGraphVertex& getVertex(vertex_descriptor);

    // Get the edge_descriptor corresponding to a segmentId.
    edge_descriptor getEdgeDescriptor(uint64_t segmentId) const;

    // Get the edge corresponding to an edge_descriptor.
    AssemblyGraphEdge& getEdge(edge_descriptor);



    // These are needed to simplify the Python API.
    void detangleVertices(TrivialDetangler&);
    void detangleEdges(TrivialDetangler&);
    void detangleVertices(PermutationDetangler&);
    void detangleEdges(PermutationDetangler&);

    // Map from AnchorId to vertex_descriptor.
    std::map<AnchorId, vertex_descriptor> vertexMap;

    // Map from edge id toedge_descriptor.
    std::map<uint64_t, edge_descriptor> edgeMap;

};

