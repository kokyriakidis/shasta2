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
    const uint64_t minCoverage = 6;

    // During creation, we store all the OrientedReadIds,
    // including duplicates, encountered for each readId0.
    MemoryMapped::VectorOfVectors<OrientedReadId, ReadId> orientedReadIds;



    class EdgePair {
    public:
        ReadId readId0 = invalid<ReadId>;
        ReadId readId1 = invalid<ReadId>;
        bool isSameStrand = false;
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
    };
    MemoryMapped::Vector<EdgePair> edgePairs;
    vector< shared_ptr<MemoryMapped::Vector<EdgePair> > > threadEdgePairs;

    void threadFunctionPass1(uint64_t threadId);
    void threadFunctionPass2(uint64_t threadId);
    void threadFunctionPass12(uint64_t pass);
    void threadFunctionPass3(uint64_t threadId);

    // Store the indexes of the edge pairs that each ReadId is involved in.
    // These are indexes relative to edgePairs.begin().
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> connectivityTable;
    void threadFunctionPass4(uint64_t threadId);
    void threadFunctionPass5(uint64_t threadId);
    void threadFunctionPass45(uint64_t pass);

    void writeGraphviz() const;

};
