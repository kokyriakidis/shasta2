// Shasta
#include "AssemblyGraphPostprocessor.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const Anchors& anchors,
    const AssemblerOptions& assemblerOptions,
    const string& assemblyStage) :
    AssemblyGraph(anchors, assemblerOptions, assemblyStage)
{
    const AssemblyGraph& assemblyGraph3 = *this;

    // Fill in the edge map.
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph) {
        edgeMap.insert(make_pair(assemblyGraph3[e].id, e));
    }

    computeAnnotations();
}


void AssemblyGraphPostprocessor::computeAnnotations()
{
    const AssemblyGraph& assemblyGraph3 = *this;
    annotations.clear();

    // Vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph3, AssemblyGraph) {
        annotations.emplace_back(assemblyGraph3[v].anchorId, v);
    }

    // Edges.
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph3[e];
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
