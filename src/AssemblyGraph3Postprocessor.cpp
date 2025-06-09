// Shasta
#include "AssemblyGraph3Postprocessor.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraph3Postprocessor::AssemblyGraph3Postprocessor(
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


void AssemblyGraph3Postprocessor::computeAnnotations()
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



span<const AssemblyGraph3Postprocessor::Annotation>
    AssemblyGraph3Postprocessor::getAnnotations(AnchorId anchorId) const
{
    const auto p = std::equal_range(annotations.begin(), annotations.end(), Annotation(anchorId));
    return span<const Annotation>(p.first, p.second);
}



// Return true if an AnchorId has one or more vertex annotations.
bool AssemblyGraph3Postprocessor::hasVertexAnnotation(AnchorId anchorId) const
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
void AssemblyGraph3Postprocessor::findAnnotationEdges(
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
