#pragma once

/*******************************************************************************

A TangleMatrix is defined by two sets of AssemblyGraph2 vertices (segments):
- The entrances.
- The exits.

The entrances and exits cannot be adjacent.

The TangleMatrix is constructed using the last step of each entrance
and the first step of each exit.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph2.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class TangleMatrix2;

}


class shasta::TangleMatrix2 {
public:
    using vertex_descriptor = AssemblyGraph2::vertex_descriptor;

    TangleMatrix2(
        const AssemblyGraph2&,
        vector<vertex_descriptor> entranceVertices,
        vector<vertex_descriptor> exitVertices,
        double aDrift,
        double bDrift);

    void writeHtml(
        const AssemblyGraph2&,
        ostream&) const;

    class EntranceOrExit {
        public:
        vertex_descriptor v;

        // The last step of this AssemblyGraph2Vertex (for an Entrance).
        // The first step of thisAssemblyGraph2Vertex (for an Exit).
        const AssemblyGraph2VertexStep& step;

        // Common coverage for this entrance or exit.
        // This is the sum of tangle matrix entries for this entrance or exit.
        uint64_t commonCoverage = 0;

        EntranceOrExit(
            vertex_descriptor v,
            const AssemblyGraph2VertexStep& step) :
            v(v), step(step) {}

        uint64_t coverage() const
        {
            return step.anchorPair.orientedReadIds.size();
        }
    };
    vector<EntranceOrExit> entrances;
    vector<EntranceOrExit> exits;

    // The tangle matrix is indexed by [iEntrance][iExit].
    // It contains the "bridge" AnchorPair between each entrance and exit.
    vector < vector<AnchorPair> > tangleMatrix;
};
