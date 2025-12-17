#pragma once

/******************************************************************

In shasta2, the ReadGraph is only used for strand separation,
during initial generation of the anchors.
Each vertex represents an OrientedReadId.
An undirected edge v0--v1 is created if the
corresponding OrientedReadIds appear together in a sufficient
number of anchors.

The ReadGraph is generated using as input the initial Anchors,
which include both strands. As a result, the ReadGraph is exactly
strand symmetric. This allows us to store only half of the edges.
We store the edges in which the lowest numbered OrientedReadId
is on strand 0. Each stored edge correspond to a pair
of reverse complemented edges.

******************************************************************/

// Shasta.
#include <invalid.hpp>
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include "memory.hpp"
#include "tuple.hpp"

namespace shasta2 {
    class ReadGraph;
    class Anchors;
}



class shasta2::ReadGraph :
    public MappedMemoryOwner,
    public MultithreadedObject<ReadGraph>
{
public:

    // Initial construction of the ReadGraph from Anchors.
    ReadGraph(
        const Anchors&,
        uint64_t threadCount);

    // Access the ReadGraph from binary data.
    ReadGraph(const Anchors&);

    const Anchors& anchors;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 1;
    const uint64_t maxKargerIterationCount = 10;

    // During creation, we store all the OrientedReadIds,
    // including duplicates, encountered for each readId0.
    MemoryMapped::VectorOfVectors<OrientedReadId, ReadId> orientedReadIds;


    // Edge Pairs are stored with readId0 < readId1.
    class EdgePair {
    public:
        ReadId readId0 = invalid<ReadId>;
        ReadId readId1 = invalid<ReadId>;
        bool isSameStrand = false;
        bool isCrossStrand = false;
        uint16_t coverage = invalid<uint16_t>;
        EdgePair() {}
        EdgePair(
            ReadId readId0,
            ReadId readId1,
            bool isSameStrand,
            uint16_t coverage) :
            readId0(readId0),
            readId1(readId1),
            isSameStrand(isSameStrand),
            coverage(coverage)
        {}

        ReadId getOther(ReadId readId) const
        {
            if(readId == readId0) {
                return readId1;
            } else if(readId == readId1) {
                return readId0;
            } else {
                SHASTA2_ASSERT(0);
            }
        }
    };
    MemoryMapped::Vector<EdgePair> edgePairs;
    vector< shared_ptr<MemoryMapped::Vector<EdgePair> > > threadEdgePairs;

    void threadFunctionPass1(uint64_t threadId);
    void threadFunctionPass2(uint64_t threadId);
    void threadFunctionPass12(uint64_t pass);
    void threadFunctionPass3(uint64_t threadId);

    // Store the indexes of the edge pairs that each ReadId is involved in.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> connectivityTable;
    void threadFunctionPass4(uint64_t threadId);
    void threadFunctionPass5(uint64_t threadId);
    void threadFunctionPass45(uint64_t pass);
    void threadFunctionPass6(uint64_t threadId);
    void writeConnectivityTable() const;



    // Find strand-symmetric quadrilaterals in the read graph.
    // A strand-symmetric quadrilateral is defined by two EdgePairs
    // with the same readId0  and readId1 (with readId0 != readId1).
    // One of the two EdgePairs has isSameStrand set to true
    // and the other one had isSameStrand set to false.
    // This implicitly defines 4 vertices and 4 edges of the ReadGraph,
    // The vertices correspond to the OrientedReadIds:
    // readId0-0 readId0-1 readId1-0 readId1-1.
    // The edges are:
    // readId0-0 readId1-0
    // readId0-1 readId1-1
    // readId0-0 readId1-1
    // readId0-1 readId1-0
    // The "diagonal" edges cannot exist because we don't allow
    // a ReadId to appear twice in the same Anchor:
    // readId0-0 readId0-1
    // readId1-0 readId1-1
    class StrandSymmetricQuadrilateral {
    public:

        // The lowest number ReadId.
        ReadId readId0;

        // The highest numbered ReadId.
        ReadId readId1;

        // The index and coverage of the EdgePair with these two ReadIds
        // and isSameStrand set to true.
        uint64_t edgePairIndexSameStrand;
        uint64_t coverageSameStrand;

        // The index and coverage of the EdgePair with these two ReadIds
        // and isSameStrand set to false.
        uint64_t edgePairIndexOppositeStrands;
        uint64_t coverageOppositeStrands;

        bool operator<(const StrandSymmetricQuadrilateral& that) const
        {
            return tie(readId0, readId1) < tie(that.readId0, that.readId1);
        }
    };
    vector<StrandSymmetricQuadrilateral> strandSymmetricQuadrilaterals;
    void findStrandSymmetriQuadrilaterals();


    // This creates the keep vector that specifies which AnchorMarkerInfos
    // should be kept for Anchors cleanup.
    // The keep vector is of size anchors.anchorMarkerInfos.totalSize() and
    // is indexed by the global position of the AnchorMarkerInfo
    // in anchorMarkerInfos, that is, &anchorMarkerInfo-anchors.anchorMarkerInfos.begin().
    // This is done by finding ReadGraph edges that are likely to cross strands,
    // then flagging to be removed the corresponding AnchorMarkerInfos.
    void anchorCleanup(vector<bool>& keep);


    // Use self-complementary paths of length 2 to find and flag cross-strand EdgePairs.
    void flagCrossStrandEdgePairs2(vector<uint64_t>& crossStrandEdgePairIndexes);

    // Same, but using self-complementary paths of length m.
    uint64_t flagCrossStrandEdgePairs(uint64_t m);



    void writeGraphviz() const;



    // Class to describe a set of connected components of the ReadGraph.
    // These components don't necessarily include the entire ReadGraph.
    // They can include a subset.
    // Each component is sorted.
    class Component :public vector<OrientedReadId> {
    public:
        bool isDoubleStranded() const;
        bool isSingleStranded() const;
    };
    class Components : public vector<Component> {
    public:

        // The id (index in the base class vector) of the component that
        // each OrientedReadId belongs to.
        // Indexed by OrientedReadId::getValue().
        vector<uint64_t> componentId;

        // This should be called when all the components have been added.
        void fillComponentId(uint64_t orientedReadCount);

    };



    // The connected components of the complete, double-stranded ReadGraph.
    Components components2;
    void computeComponents2();



    // We compute a single-stranded version of the ReadGraph as follows:
    // - For each reverse complemented pair of single-stranded connected
    //   component, we only keep one (the one in which the lowest numbered
    //   OrientedReadId is on strand 0).
    // - For each double-stranded component, we do a strand-aware
    //   approximate min-cut to approximately separate strands.
    //   The cut separates that component in two single-stranded
    //   components and we keep only one of them (the one in which
    //   the lowest numbered OrientedReadId is on strand 0).
    Components components1;
    void computeComponents1();


    // Experiments with local strand-aware min-cuts.
    void localMinCut(ReadId, vector<uint64_t>& cutEdgePairIndexes) const;
};
