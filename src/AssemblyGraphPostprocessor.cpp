// Shasta.
#include "AssemblyGraphPostprocessor.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraphPostprocessor::AssemblyGraphPostprocessor(
    const Anchors& anchors,
    const string& assemblyStage) :
    AssemblyGraph(anchors, assemblyStage)
{
    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        segmentMap.insert(make_pair(assemblyGraph[e].id, e));
    }
}
