#pragma once

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"

// Standard library.
#include "vector.hpp"



namespace shasta2 {

    class Bubble;
    class BubbleChain;
}



// A Bubble is a set of parallel edges in the AssemblyGraph.
class shasta2::Bubble {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AssemblyGraphBaseClass::edge_descriptor;

    vertex_descriptor v0;
    vertex_descriptor v1;
    vector<edge_descriptor> edges;

    Bubble() {}

    // Generate a Bubble given its source and target vertices.
    Bubble(const AssemblyGraph&, vertex_descriptor v0, vertex_descriptor v1);

    // This returns the maximum of the edge lengths.
    uint64_t maxLength(const AssemblyGraph&) const;
};



class shasta2::BubbleChain : public vector<Bubble> {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;

    // Generate a BubbleChain given a vector of AssemblyGraph vertices.
    BubbleChain(const AssemblyGraph&, const vector<vertex_descriptor>&);

    // This returns the sum of the maxLengths of all the bubbles.
    uint64_t maxLength(const AssemblyGraph&) const;
};
