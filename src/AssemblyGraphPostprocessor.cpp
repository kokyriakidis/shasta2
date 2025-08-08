// Shasta
#include "AssemblyGraphPostprocessor.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const Anchors& anchors,
    const Journeys& journeys,
    const Options& options,
    const string& assemblyStage) :
    AssemblyGraph(anchors, journeys, options, assemblyStage)
{
    const AssemblyGraph& assemblyGraph = *this;

    // Fill in the vertex map.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert(make_pair(assemblyGraph[v].id, v));
    }

    // Fill in the edge map.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgeMap.insert(make_pair(assemblyGraph[e].id, e));
    }

    computeAnnotations();
}


void AssemblyGraphPostprocessor::computeAnnotations()
{
    const AssemblyGraph& assemblyGraph = *this;
    annotations.clear();

    // Vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        annotations.emplace_back(assemblyGraph[v].anchorId, v);
    }

    // Edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t step=0; step<edge.size(); step++) {
            const AnchorPair& anchorPair = edge[step].anchorPair;
            annotations.emplace_back(anchorPair.anchorIdA, e, step, true);
            annotations.emplace_back(anchorPair.anchorIdB, e, step, false);
        }
    }

    sort(annotations.begin(), annotations.end());
}



span<const AssemblyGraphPostprocessor::Annotation>
    AssemblyGraphPostprocessor::getAnnotations(AnchorId anchorId) const
{
    const auto p = std::equal_range(annotations.begin(), annotations.end(), Annotation(anchorId));
    return span<const Annotation>(p.first, p.second);
}



// Return true if an AnchorId has one or more vertex annotations.
bool AssemblyGraphPostprocessor::hasVertexAnnotation(AnchorId anchorId) const
{
    const auto annotations = getAnnotations(anchorId);
    for(const Annotation& annotation: annotations) {
        if(annotation.v != null_vertex()) {
            return true;
        }
    }
    return false;
}



// Find the edges that an AnchorId has annotations for.
void AssemblyGraphPostprocessor::findAnnotationEdges(
    AnchorId anchorId,
    vector<edge_descriptor>& edges) const
{
    edges.clear();

    const auto annotations = getAnnotations(anchorId);
    for(const Annotation& annotation: annotations) {
        if(annotation.v == null_vertex()) {
            edges.push_back(annotation.e);
        }
    }
    deduplicate(edges);
}
