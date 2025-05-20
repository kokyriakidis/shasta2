#pragma once

#include "AssemblyGraph3.hpp"


namespace shasta {
    class AssemblyGraph3Postprocessor;
}



// AssemblyGraph3 functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class shasta::AssemblyGraph3Postprocessor : public AssemblyGraph3 {
public:
    AssemblyGraph3Postprocessor(
        const Anchors&,
        const AssemblerOptions&,
        const string& assemblyStage);

    // Map from edge id to edge_descriptor.
    std::map<uint64_t, edge_descriptor> edgeMap;

};

