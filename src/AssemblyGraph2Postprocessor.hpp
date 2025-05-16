#pragma once

#include "AssemblyGraph2.hpp"


namespace shasta {
    class AssemblyGraph2Postprocessor;
    class TrivialDetangler;
}



// AssemblyGraph2 functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class shasta::AssemblyGraph2Postprocessor : public AssemblyGraph2 {
public:
    AssemblyGraph2Postprocessor(
        const Anchors&,
        const AssemblerOptions&,
        const string& assemblyStage);

    // Map from vertex id to vertex_descriptor.

    std::map<uint64_t, vertex_descriptor> vertexMap;


    // These are needed to simplify the Python API.
    void detangleVertices(TrivialDetangler&);
};

