#pragma once

// A Superbubble is similar a set of AssemblyGraph vertices such that:
// - All the in-edges of Superbubble vertices are into a single target vertex,
//   which is the source of the Superbubble.
// - All the out-edges of Superbubble vertices are out of a source vertex,
//   which is the target of the Superbubble.
// Usually a Superbubble has only one entrance edge (into the source)
// and one exit edge (from the target), but this is not necessarily the case.
// If the Superbubble has no internal vertices, it is a bubble.
// In that case the out-degree of the source must be equal to the
// in-degree of the target, and defines the ploidy of the bubble.

#include "AssemblyGraph.hpp"

namespace shasta {
    class Superbubble;
}


class shasta::Superbubble {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    const AssemblyGraph& assemblyGraph;

    vertex_descriptor sourceVertex;
    vertex_descriptor targetVertex;

    // The internal vertices are stored sorted by id so we can do binary searches on it.
    // They do not include the source and target vertices.
    vector<vertex_descriptor> internalVertices;
    void gatherInternalVertices();
    bool contains(vertex_descriptor v) const {
        return std::binary_search(internalVertices.begin(), internalVertices.end(), v);
    }

    // The source edges are the out-edges of the source vertex, stored sorted by id.
    // The target edges are the in-edges of the target vertex, stored sorted by id.
    // All source and target edges are also stored in the internalEdges vector.
    // An edge can be at the same time a source and edge and a target edge.
    vector<edge_descriptor> sourceEdges;
    vector<edge_descriptor> targetEdges;

    // The internal edges are stored sorted by id.
    // They include the source and target edges.
    vector<edge_descriptor> internalEdges;

    // Gather the source, target, and internal edges.
    void gatherEdges();

    // If there are no internal vertices, this Superbubble is a bubble.
    bool isBubble() const {
        return internalVertices.empty();
    }

    const uint64_t sourcePloidy() const
    {
        return sourceEdges.size();
    }

    const uint64_t targetPloidy() const
    {
        return targetEdges.size();
    }

    const uint64_t ploidy() const {
        SHASTA_ASSERT(isBubble());
        const uint64_t ploidyAtSource = sourcePloidy();
        const uint64_t ploidyAtTarget = targetPloidy();
        SHASTA_ASSERT(ploidyAtSource == ploidyAtTarget);
        return ploidyAtSource;
    }

    bool isTrivial() const
    {
        return internalVertices.empty() and internalEdges.size() == 1;
    }

    Superbubble(
        const AssemblyGraph&,
        vertex_descriptor sourceVertex,
        vertex_descriptor targetVertex);

};
