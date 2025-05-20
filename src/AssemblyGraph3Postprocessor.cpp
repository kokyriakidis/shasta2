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

}
