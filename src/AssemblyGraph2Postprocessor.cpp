// Shasta.
#include "AssemblyGraph2Postprocessor.hpp"
#include "Base.hpp"
#include "TrivialDetangler.hpp"

using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



AssemblyGraph2Postprocessor::AssemblyGraph2Postprocessor(
    const Anchors& anchors,
    const AssemblerOptions& assemblerOptions,
    const string& assemblyStage) :
    AssemblyGraph2(anchors, assemblerOptions, assemblyStage)
{
    const AssemblyGraph2& assemblyGraph2 = *this;

    // Fill in the vertex map.
    BGL_FORALL_VERTICES(v, assemblyGraph2, AssemblyGraph2) {
        vertexMap.insert(make_pair(assemblyGraph2[v].id, v));
    }

}


// These are needed to simplify the Python API.
void AssemblyGraph2Postprocessor::detangleVertices(TrivialDetangler& detangler) {
    AssemblyGraph2::detangleVertices(detangler);
}
