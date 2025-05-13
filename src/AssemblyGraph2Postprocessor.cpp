// Shasta.
#include "AssemblyGraph2Postprocessor.hpp"
#include "Base.hpp"
using namespace shasta;



AssemblyGraph2Postprocessor::AssemblyGraph2Postprocessor(
    const Anchors& anchors,
    const AssemblerOptions& assemblerOptions,
    const string& assemblyStage) :
    AssemblyGraph2(anchors, assemblerOptions, assemblyStage)
{
}
