#pragma once

/*******************************************************************************

A TangleMatrix is defined by two sets of AssemblyGraph edges (segments):
- The entrances.
- The exits.
Any two sets of edges can be used to define a TangleMatrix,
regardless of their location and connectivity in the AssemblyGraph.
The two sets are not required to be disjoint.

We use the last AssemblyGraphStep of each entrances and the first
AssemblyGraphStep of each exit to count common reads.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class TangleMatrix;

}


class shasta::TangleMatrix {
public:
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    TangleMatrix(
        const AssemblyGraph&,
        vector<edge_descriptor> entranceEdges,
        vector<edge_descriptor> exitEdges);

    TangleMatrix() {}

    void construct(
        const AssemblyGraph&,
        vector<edge_descriptor> entranceEdges,
        vector<edge_descriptor> exitEdges);

    void writeHtml(
        const AssemblyGraph&,
        ostream&) const;

    class Entrance {
        public:
        edge_descriptor e;

        // The last AssemblyGraphStep of this AssemblyGraphEdge,
        // after removing OrientedReadIds that also appear in other Entrances.
        AssemblyGraphStep step;

        Entrance(
            edge_descriptor e,
            const AssemblyGraphStep& step,
            const vector<OrientedReadId>& duplicateOrientedReadIdsOnEntrances) :
            e(e), step(step, duplicateOrientedReadIdsOnEntrances) {}

    };
    vector<Entrance> entrances;



    class Exit {
        public:
        edge_descriptor e;

        // The first AssemblyGraphStep of this AssemblyGraphEdge,
        // after removing OrientedReadIds that also appear in other Exits.
        AssemblyGraphStep step;

        Exit(
            edge_descriptor e,
            const AssemblyGraphStep& step,
            const vector<OrientedReadId>& duplicateOrientedReadIdsOnExits) :
            e(e), step(step, duplicateOrientedReadIdsOnExits) {}

        // The OrientedReadIds in this exit that don't appear in other exits.
        vector<OrientedReadId> orientedReadIds;
    };
    vector<Exit> exits;

    // The OrientedReadIds that appear in more than one Entrance or Exit.
    // These are not counted when constructing the TangleMatrix.
    vector<OrientedReadId> duplicateOrientedReadIdsOnEntrances;
    vector<OrientedReadId> duplicateOrientedReadIdsOnExits;

    // The tangle matrix is indexed by [iEntrance][iExit].
    vector < vector<uint64_t> > tangleMatrix;
};
