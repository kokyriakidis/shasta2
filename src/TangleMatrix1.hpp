#pragma once

/*******************************************************************************

A TangleMatrix1 is defined by two sets of AssemblyGraphe edges (segments):
- The entrances.
- The exits.

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



    // For each entrance, we define a representative region
    // consisting of the last representativeRegionLength
    // step in the entrance edge.
    // For each exit, we define a representative region
    // consisting of the first representativeRegionLength
    // step in the exit edge.
    // We gather information for each oriented read that appear
    // in these representative regions.



    // The common OrientedReadIds are the ones that appear
    // in the representative regions of at least one entrance and
    // at leas one exit. They are the only ones that contribute
    // to the tangle matrix.
    class CommonOrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The number of times this OrientedReadId appears
        // in the representative region of each entrance.
        vector<uint64_t> entranceStepCount;

        // The number of times this OrientedReadId appears
        // in the representative region of each exit.
        vector<uint64_t> exitStepCount;

        CommonOrientedReadInfo(
            OrientedReadId orientedReadId,
            uint64_t entranceCount = 0,
            uint64_t exitCount = 0);

        // Order them by OrientedReadId.
        bool operator<(const CommonOrientedReadInfo& that) const;

        // The contribution of this oriented read to the total tangle matrix.
        vector< vector<double> > tangleMatrix;
        void computeTangleMatrix();
    };
    vector<CommonOrientedReadInfo> commonOrientedReadInfos;
    void gatherCommonOrientedReads();
    void writeCommonOrientedReads(ostream& html) const;



    // Return the index of a given OrientedReadId in the commonOrientedReadInfos vector,
    // or invalid<uint64_t> if not present.
    uint64_t getCommonOrientedReadIdIndex(OrientedReadId orientedReadId) const;

    // The total tangle matrix.
    // It is computed by adding the contributions of each of the common
    // oriented reads.
    vector< vector<double> > tangleMatrix;
    void computeTotalTangleMatrix();
    void writeTotalTangleMatrix(ostream& html) const;
};
