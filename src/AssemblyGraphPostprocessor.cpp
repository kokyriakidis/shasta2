// Shasta.
#include "AssemblyGraphPostprocessor.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const string& assemblyStage,
    const Anchors& anchors
    ) :
    AssemblyGraph(assemblyStage, anchors)
{
    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        segmentMap.insert(make_pair(assemblyGraph[e].id, e));
    }
}
