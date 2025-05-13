#pragma once

#include "AssemblyGraph2.hpp"


namespace shasta {
    class AssemblyGraph2Postprocessor;
}



// AssemblyGraph2 functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class shasta::AssemblyGraph2Postprocessor : public AssemblyGraph2 {
public:
    AssemblyGraph2Postprocessor(
        const Anchors&,
        const AssemblerOptions&,
        const string& assemblyStage);

};

