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



    // We define a representative region for each entrance and exit.
    // For an entrance, the representative region
    // consists of the last representativeRegionLength
    // step in the entrance edge.
    // For an exit, the representative region
    // consists of the first representativeRegionLength
    // step in the exit edge.
    // We gather information for each oriented read that appear
    // in these representative regions.
    // These are stored sorted by OrientedReadId.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The number of steps where this OrientedReadId is present
        // in the representative region of the entrance or exit
        // this refers to.
        uint64_t stepCount = 0;

        // Position in journey.
        // For an entrance, this is the position in journey of this orientedReadId
        // at anchorIdB of the last step in which this oriented read appears.
        // For an exit, this is the position in journey of this orientedReadId
        // at anchorIdA of the first step in which this oriented read appears.
        // This information is used in detangling.
        uint32_t positionInJourney;

        bool operator<(const OrientedReadInfo&) const;
    };
    vector< vector<OrientedReadInfo> > entranceOrientedReadInfos;
    vector< vector<OrientedReadInfo> > exitOrientedReadInfos;
    void gatherOrientedReads(uint64_t representativeRegionLength);
    void gatherEntranceOrientedReads(uint64_t iEntrance, uint64_t representativeRegionLength);
    void gatherExitOrientedReads(uint64_t iExit, uint64_t representativeRegionLength);
    void writeOrientedReads(ostream& html) const;



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
        bool operator<(const CommonOrientedReadInfo&) const;

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
