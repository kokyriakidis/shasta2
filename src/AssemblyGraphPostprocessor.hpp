#pragma once

#include "AssemblyGraph.hpp"


namespace shasta {
    class AssemblyGraphPostprocessor;
}



// AssemblyGraph functionality needed only during postprocessing.
class shasta::AssemblyGraphPostprocessor : public AssemblyGraph {
public:
    AssemblyGraphPostprocessor(
        const string& assemblyStage,
        const Anchors& anchors);
};

