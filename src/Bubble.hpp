#pragma once

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"

// Standard library.
#include "vector.hpp"



namespace shasta2 {

    class Bubble;
}



// A Bubble is a set of parallel edges in the AssemblyGraph.
class shasta2::Bubble {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AssemblyGraphBaseClass::edge_descriptor;

    vertex_descriptor v0;
    vertex_descriptor v1;
    vector<edge_descriptor> edges;
};
