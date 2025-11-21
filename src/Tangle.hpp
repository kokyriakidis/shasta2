#pragma once

// Shasta.
#include "AssemblyGraph.hpp"

// Standard library.
#include "memory.hpp"
#include "vector.hpp"

namespace shasta {
    class Tangle;
    class TangleMatrix;
}



class shasta::Tangle {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    // Constructor from a set of AssemblyGraph vertices.
    Tangle(
        AssemblyGraph&,
        const vector<vertex_descriptor>&,
        uint64_t maxTrim,
        double aDrift,
        double bDrift);

    // Constructor from a single AssemblyGraph vertex.
    Tangle(
        AssemblyGraph&,
        vertex_descriptor,
        uint64_t maxTrim,
        double aDrift,
        double bDrift);

    // Constructor from a single AssemblyGraph edge.
    Tangle(
        AssemblyGraph&,
        edge_descriptor,
        uint64_t maxTrim,
        double aDrift,
        double bDrift);

    shared_ptr<TangleMatrix> tangleMatrix;

    // Detangling instructions.
    // Each pair is (entranceIndex, exitIndex) that are to be connected
    // when detangling.
    // Detangling decisions are not made in Tangle. They are made by the
    // Detangler object.
    vector< pair<uint64_t, uint64_t> > connectList;
    void connect(uint64_t iEntrance, uint64_t iExit);
    void detangle();

    // If the Tangle is detangled successfully, we store the vertices that were removed.
    vector<vertex_descriptor> removedVertices;

    bool debug = false;

    AssemblyGraph& assemblyGraph;

    // The Tangle vertices, sorted by assemblyGraph.orderById.
    vector<vertex_descriptor> tangleVertices;
    bool isTangleVertex(vertex_descriptor v) const
    {
        return std::binary_search(tangleVertices.begin(), tangleVertices.end(), v, assemblyGraph.orderById);
    }

    void computeExtendedTangleMatrix(vector< vector<double> >&) const;

};
