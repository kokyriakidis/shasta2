#pragma once

/*******************************************************************************

A TangleMatrix1 is defined by two sets of AssemblyGraphe edges (segments):
- The entrances.
- The exits.

The TangleMatrix1 is constructed using AssemblyGraph::orientedReadSegments
and AssemblyGraphEdge::transitioningOrientedReadIds.
These are filled in by  AssemblyGraph::countOrientedReadStepsBySegment
and should not be out of date.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class TangleMatrix1;
}


class shasta::TangleMatrix1 {
public:
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    TangleMatrix1(
        const AssemblyGraph&,
        vector<edge_descriptor> entrances,
        vector<edge_descriptor> exits,
        ostream& html);

    const AssemblyGraph& assemblyGraph;
    const vector<edge_descriptor>& entrances;
    const vector<edge_descriptor>& exits;


    // Some information for each of  OrientedReadId that contributes to
    // this tangle matrix. Sorted by OrientedReadId.
    // These oriented reads must appear in at least one entrance and
    // at least one exit.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        vector<uint64_t> entranceStepCount;
        vector<uint64_t> exitStepCount;
        OrientedReadInfo(
            OrientedReadId orientedReadId,
            uint64_t entranceCount = 0,
            uint64_t exitCount = 0);
        bool operator<(const OrientedReadInfo& that) const;

        // The contribution of this oriented read to the total tangle matrix.
        vector< vector<double> > tangleMatrix;
        void computeTangleMatrix();
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();
    void writeOrientedReads(ostream& html) const;

    // Return the index of a given OrientedReadId in the orientedReadInfos vector,
    // or invalid<uint64_t> if not present.
    uint64_t getOrientedReadIdIndex(OrientedReadId orientedReadId) const;

    // The total tangle matrix.
    vector< vector<double> > tangleMatrix;
    void computeTotalTangleMatrix();
    void writeTotalTangleMatrix(ostream& html) const;
};
