#pragma once

// Shasta.
#include "AssemblyGraph3.hpp"

// Standard library.
#include <memory.hpp>
#include <vector.hpp>

namespace shasta {
    class Tangle3;
    class TangleMatrix3;
}



class shasta::Tangle3 {
public:
    using vertex_descriptor = AssemblyGraph3::vertex_descriptor;
    using edge_descriptor = AssemblyGraph3::edge_descriptor;

    // Constructor from a set of AssemblyGraph3 vertices.
    Tangle3(
        AssemblyGraph3&,
        const vector<vertex_descriptor>&,
        double aDrift,
        double bDrift);

    // Constructor from a single AssemblyGraph3 edge.
    Tangle3(
        AssemblyGraph3&,
        edge_descriptor,
        double aDrift,
        double bDrift);

    shared_ptr<TangleMatrix3> tangleMatrix;

    // Detangling instructions.
    // Each pair is (entranceIndex, exitIndex) that are to be connected
    // when detangling.
    // Detangling decisions are not made in Tangle2. They are made by the
    // Detangler2 object.
    vector< pair<uint64_t, uint64_t> > connectList;
    void connect(uint64_t iEntrance, uint64_t iExit);
    void detangle();

    // If the Tangle3 is detangled successfully, we store the vertices that were removed.
    vector<vertex_descriptor> removedVertices;

    bool debug = false;

private:
    AssemblyGraph3& assemblyGraph3;

    // The Tangle3 vertices.
    vector<vertex_descriptor> tangleVertices;
    bool isTangleVertex(vertex_descriptor v) const
    {
        return std::binary_search(tangleVertices.begin(), tangleVertices.end(), v);
    }

};
