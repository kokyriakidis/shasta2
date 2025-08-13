#pragma once

// Shasta.
#include "AssemblyGraph.hpp"

// Standard library.
#include <memory.hpp>
#include <vector.hpp>

namespace shasta {
    class Tangle1;
    class TangleMatrix1;
}



class shasta::Tangle1 {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    // Constructor from a set of AssemblyGraph vertices.
    Tangle1(AssemblyGraph&, const vector<vertex_descriptor>&);

    // Constructor from a single AssemblyGraph vertex.
    Tangle1(AssemblyGraph&, vertex_descriptor);

    // Constructor from a single AssemblyGraph edge.
    Tangle1(AssemblyGraph&, edge_descriptor);

    AssemblyGraph& assemblyGraph;

    // All vectors of vertices and edges are stored sorted by id.

    // The tangle vertices.
    vector<vertex_descriptor> tangleVertices;
    bool isTangleVertex(vertex_descriptor) const;

    // Entrances are edges with source outside the tangle
    // and target inside the tangle.
    vector<edge_descriptor> entrances;
    void findEntrances();

    // Exits are edges with source inside the tangle
    // and target outside the tangle.
    vector<edge_descriptor> exits;
    void findExits();

    // Internal edges are edges with source and target inside the tangle.
    vector<edge_descriptor> tangleEdges;
    void findTangleEdges();

    // The tangle matrix constructed using the entrances and exits.
    shared_ptr<TangleMatrix1> tangleMatrixPointer;
    const TangleMatrix1& tangleMatrix() const
    {
        return *tangleMatrixPointer;
    }

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

    void rerouteEntrances(vector<vertex_descriptor>& newEntranceVertices) const;
    void rerouteExits(vector<vertex_descriptor>& newExitVertices) const;
    void reconnect(
        uint64_t iEntrance,
        uint64_t iExit,
        vertex_descriptor v0,
        vertex_descriptor v1
        ) const;
};
