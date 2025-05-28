// Shasta
#include "AssemblyGraph3Postprocessor.hpp"
#include "Base.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraph3Postprocessor::AssemblyGraph3Postprocessor(
    const Anchors& anchors,
    const AssemblerOptions& assemblerOptions,
    const string& assemblyStage) :
    AssemblyGraph3(anchors, assemblerOptions, assemblyStage)
{
    const AssemblyGraph3& assemblyGraph3 = *this;

    // Fill in the edge map.
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        edgeMap.insert(make_pair(assemblyGraph3[e].id, e));
    }

    computeAnnotations();
}


void AssemblyGraph3Postprocessor::computeAnnotations()
{
    const AssemblyGraph3& assemblyGraph3 = *this;
    annotations.clear();

    // Vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph3, AssemblyGraph3) {
        annotations.emplace_back(assemblyGraph3[v].anchorId, v);
    }

    // Edges.
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        const AssemblyGraph3Edge& edge = assemblyGraph3[e];
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
