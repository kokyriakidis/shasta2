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
        vector<edge_descriptor> entranceEdges,
        vector<edge_descriptor> exitEdges,
        uint64_t maxTrim,
        double aDrift,
        double bDrift);

    void writeHtml(
        const AssemblyGraph&,
        ostream&) const;

    class EntranceOrExit {
        public:
        edge_descriptor e;

        // The number of steps that were trimmed at the end of an entrance
        // or at the beginning of an exit.
        uint64_t trim;

        // The last step of this AssemblyGraphEdge (for an Entrance).
        // The first step of this AssemblyGraphEdge (for an Exit).
        const AssemblyGraphEdgeStep& step;

        // Common coverage for this entrance or exit.
        // This is the sum of tangle matrix entries for this entrance or exit.
        uint64_t commonCoverage = 0;

        EntranceOrExit(
            edge_descriptor e,
            uint64_t trim,
            const AssemblyGraphEdgeStep& step) :
            e(e), trim(trim), step(step) {}

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

    void getTangleMatrixCoverage(vector< vector<uint64_t> >& tangleMatrixCoverage) const
    {
        tangleMatrixCoverage.resize(entrances.size(), vector<uint64_t>(exits.size()));
        for(uint64_t i=0; i<entrances.size(); i++) {
            for(uint64_t j=0; j<exits.size(); j++) {
                tangleMatrixCoverage[i][j] = tangleMatrix[i][j].size();
            }
        }
    }



    // Read following on the entrances/exits.
    class StepIdentifier {
    public:
        uint64_t iEntrance = invalid<uint64_t>;
        uint64_t iExit = invalid<uint64_t>;
        uint64_t stepId;
    };
    std::map<OrientedReadId, vector<StepIdentifier> > readFollowingMap;
    void readFollowing(const AssemblyGraph&);
};
