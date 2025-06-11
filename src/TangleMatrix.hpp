#pragma once

/*******************************************************************************

A TangleMatrix is defined by two sets of AssemblyGraphe edges (segments):
- The entrances.
- The exits.

The TangleMatrix is constructed using the last step of each entrance
and the first step of each exit.

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
        vector<edge_descriptor> entranceVertices,
        vector<edge_descriptor> exitVertices,
        double aDrift,
        double bDrift);

    void writeHtml(
        const AssemblyGraph&,
        ostream&) const;

    class EntranceOrExit {
        public:
        edge_descriptor e;

        // The last step of this AssemblyGraphEdge (for an Entrance).
        // The first step of this AssemblyGraphEdge (for an Exit).
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



    // Likelihood ratio test of the tangle matrix (G test).
    // https://en.wikipedia.org/wiki/G-test
    class Hypothesis {
    public:
        vector< vector<bool> > connectivityMatrix;
        double G;

        Hypothesis(
            const vector< vector<bool> >& connectivityMatrix,
            double G) :
            connectivityMatrix(connectivityMatrix),
            G(G)
            {}

        // Sort by G2.
        bool operator<(const Hypothesis& that) const {
            return G < that.G;
        }
    };
    vector<Hypothesis> hypotheses;
    bool gTest(double epsilon);
};
