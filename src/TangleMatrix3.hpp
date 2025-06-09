#pragma once

/*******************************************************************************

A TangleMatrix is defined by two sets of AssemblyGraphe edges (segments):
- The entrances.
- The exits.

The TangleMatrix is constructed using the last step of each entrance
and the first step of each exit.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph3.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class TangleMatrix3;

}


class shasta::TangleMatrix3 {
public:
    using edge_descriptor = AssemblyGraph3::edge_descriptor;

    TangleMatrix3(
        const AssemblyGraph3&,
        vector<edge_descriptor> entranceVertices,
        vector<edge_descriptor> exitVertices,
        double aDrift,
        double bDrift);

    void writeHtml(
        const AssemblyGraph3&,
        ostream&) const;

    class EntranceOrExit {
        public:
        edge_descriptor e;

        // The last step of this AssemblyGraphEdge (for an Entrance).
        // The first step of thisAssemblyGraph3Edge (for an Exit).
        const AssemblyGraphEdgeStep& step;

        // Common coverage for this entrance or exit.
        // This is the sum of tangle matrix entries for this entrance or exit.
        uint64_t commonCoverage = 0;

        EntranceOrExit(
            edge_descriptor e,
            const AssemblyGraphEdgeStep& step) :
            e(e), step(step) {}

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
