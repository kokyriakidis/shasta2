#pragma once


/*******************************************************************************

A Tangle is an induced subgraph of the AssemblyGraph.
That is, a subset of its vertices and the edges whose source and target
vertex are both in that subset.

Edges that have only the source or target vertex in the Tangle are not part
of the Tangle.
 - Edges whose source vertex is outside the Tangle and target vertex is in the
   Tangle are the Entrances of the Tangle.
 - Edges whose source vertex is in the Tangle and target vertex is outside the
   Tangle are the Exits of the Tangle.

The set of vertices that define the Tangle is arbitrary. However, in useful
situations the following is true:
 - The number of vertices is in the Tangle is small (as small as a single vertex).
 - The Tangle is connected.
 - There are at least two Entrances and two Exits. The most useful case
   is when there are exactly two Entrances and two Exits.

The Entrances and Exits are used to compute a TangleMatrix.

Detangling decisions are not made by the Tangle, but by a higher level
object that calls Tangle::connect to specify Entrance/Exit pairs that should
be connected. After all such pairs have been specified, the higher level
object calls Tangle::detangle to detangle as described by these pairs.
For each pair, the last AssemblyGrapStep of the Entrance and the
first AssemblyGraphStep are merged.

Two important special cases:
 - A Tangle consisting of a single vertex.
 - A Tangle consisting of the source and target vertices of a given edge.
   The Tangle then by definition also includes that edge.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"
#include "TangleMatrix.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class Tangle;
}



class shasta::Tangle {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    // Constructor from a set of AssemblyGraph vertices.
    Tangle(
        const AssemblyGraph&,
        const vector<vertex_descriptor>&);

    // Constructor from a single AssemblyGraph vertex.
    Tangle(
        const AssemblyGraph&,
        vertex_descriptor);

    // Constructor from a single AssemblyGraph edge.
    Tangle(
        const AssemblyGraph&,
        edge_descriptor);

    TangleMatrix tangleMatrix;

private:
    void construct(
        const AssemblyGraph&,
        const vector<vertex_descriptor>&);

    // The Tangle vertices, sorted using AssemblyGraphVertexOrderByAnchorId.
    vector<vertex_descriptor> tangleVertices;
    AssemblyGraphVertexOrderByAnchorId sorter;
    bool isTangleVertex(vertex_descriptor v) const
    {
        return std::binary_search(tangleVertices.begin(), tangleVertices.end(), v, sorter);
    }

};
