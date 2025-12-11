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
is on strand 0.

******************************************************************/

// Shasta.
#include <invalid.hpp>
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

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

private:
    const Anchors& anchors;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 6;



    // We only store edges in which the lowest number ReadId, readId0,
    // is on strand 0. Indexed by readId0.
    // Each of the stored edges corresponds to a pair of reverse
    // complemented edges in the ReadGraph.
    class Edge {
    public:
        OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);;
        uint32_t coverage = invalid<uint32_t>;
        Edge() {}
        Edge(
            OrientedReadId orientedReadId,
            uint32_t coverage) :
            orientedReadId(orientedReadId),
            coverage(coverage)
        {}
    };
    MemoryMapped::VectorOfVectors<Edge, ReadId> edges;



    // During creation of the ReadGraph, we store all the OrientedReadIds,
    // including duplicates, encountered for each readId0.
    MemoryMapped::VectorOfVectors<OrientedReadId, ReadId> orientedReadIds;

    void threadFunctionPass1(uint64_t threadId);
    void threadFunctionPass2(uint64_t threadId);
    void threadFunctionPass12(uint64_t pass);
    void threadFunctionPass3(uint64_t threadId);
    void threadFunctionPass4(uint64_t threadId);
    void threadFunctionPass34(uint64_t pass);
};
