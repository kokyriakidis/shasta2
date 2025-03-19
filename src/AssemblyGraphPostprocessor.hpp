#pragma once

#include "AssemblyGraph.hpp"


namespace shasta {
    class AssemblyGraphPostprocessor;
}



// AssemblyGraph functionality needed only during postprocessing.
class shasta::AssemblyGraphPostprocessor : public AssemblyGraph {
public:
    AssemblyGraphPostprocessor(
        const Anchors&,
        const string& assemblyStage);

    // Map a segment id to an edge decriptor.
    std::map<uint64_t, edge_descriptor> segmentMap;
};

