#pragma once

// Shasta.
#include "AssemblyGraph.hpp"
#include "TrivialDetangler.hpp"

// Standard library.
#include <span.hpp>


namespace shasta {
    class AssemblyGraphPostprocessor;
}



// AssemblyGraph functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class shasta::AssemblyGraphPostprocessor : public AssemblyGraph {
public:
    AssemblyGraphPostprocessor(
        const Anchors&,
        const Journeys&,
        const AssemblerOptions&,
        const string& assemblyStage);

    // Map from edge id to edge_descriptor.
    std::map<uint64_t, edge_descriptor> edgeMap;



    // Annotations of where each AnchorId is used in the current state
    // of AssemblyGraph. An AnchorId can be used in three ways:
    // - In an AssemblyGraphVertex.
    // - In an AssemblyGraphEdge::anchorPair, as anchorIdA.
    // - In an AssemblyGraphEdge::anchorPair, as anchorIdB.
    // An AnchorId can be used any number of times (0, 1, or multiple times).
    // For the case of an AssemblyGraphVertex, the vertex_descriptor v is not
    // null vertex and contains the vertex descriptor that uses the AnchorId.
    // The other fields are invalid.
    // For the other two cases, the vertex_descriptor is null_vertex.
    class Annotation {
    public:
        AnchorId anchorId;
        vertex_descriptor v = null_vertex();
        edge_descriptor e;
        uint64_t step = invalid<uint64_t>;
        bool isAnchorIdA;

        bool operator<(const Annotation& that) const
        {
            return anchorId < that.anchorId;
        }

        Annotation(AnchorId anchorId) :
            anchorId(anchorId)
        {}
        Annotation(AnchorId anchorId,
            vertex_descriptor v) :
            anchorId(anchorId),
            v(v)
        {}
        Annotation(AnchorId anchorId,
            edge_descriptor e,
            uint64_t step,
            bool isAnchorIdA) :
            anchorId(anchorId),
            e(e),
            step(step),
            isAnchorIdA(isAnchorIdA)
        {}
    };

    // The vector of annotations is kept sorted so we can do searches.
    vector<Annotation> annotations;
    void computeAnnotations();
    span<const Annotation> getAnnotations(AnchorId) const;

    // Return true if an AnchorId has one or more vertex annotations.
    bool hasVertexAnnotation(AnchorId) const;

    // Find the edges that an AnchorId has annotations for.
    void findAnnotationEdges(AnchorId, vector<edge_descriptor>&) const;

};

