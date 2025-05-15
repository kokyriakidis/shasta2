#pragma once

// Shasta.
#include "AssemblyGraph2.hpp"

// Standard library.
#include <memory.hpp>
#include <vector.hpp>

namespace shasta {
    class Tangle2;
    class TangleMatrix2;
}



class shasta::Tangle2 {
public:
    using vertex_descriptor = AssemblyGraph2::vertex_descriptor;

    // Constructor from a set of AssemblyGraph2 vertices.
    Tangle2(
        AssemblyGraph2&,
        const vector<vertex_descriptor>&,
        double aDrift,
        double bDrift);

    // Constructor from a single AssemblyGraph2 vertex.
    Tangle2(
        AssemblyGraph2&,
        vertex_descriptor,
        double aDrift,
        double bDrift);

    shared_ptr<TangleMatrix2> tangleMatrix;

    // Detangling instructions.
    // Each pair is (entranceIndex, exitIndex) that are to be connected
    // when detangling.
    // Detangling decisions are not made in Tangle2. They are made by the
    // Detangler2 object.
    vector< pair<uint64_t, uint64_t> > connectList;
    void connect(uint64_t iEntrance, uint64_t iExit);
    void detangle();

    // If the Tangle is detangled successfully, we store the vertices that were removed.
    vector<vertex_descriptor> removedVertices;

    bool debug = false;

private:
    AssemblyGraph2& assemblyGraph2;

    // The Tangle2 vertices.
    vector<vertex_descriptor> tangleVertices;
    bool isTangleVertex(vertex_descriptor v) const
    {
        return std::binary_search(tangleVertices.begin(), tangleVertices.end(), v);
    }

};
