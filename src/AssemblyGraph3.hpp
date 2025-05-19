#pragma once

// In the AssemblyGraph3, each edge is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of edges in the AnchorGraph.
// Each Assembly3GraphEdge generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>


namespace shasta {

    class AssemblyGraph3;
    class AssemblyGraph3Vertex;
    class AssemblyGraph3Edge;
    class AssemblyGraph3EdgeStep;

    using AssemblyGraph3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraph3Vertex,
        AssemblyGraph3Edge>;

    class AnchorGraph;
    class Anchors;
    class AssemblerOptions;

}



// Tha AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
class shasta::AssemblyGraph3Vertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    AssemblyGraph3Vertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
};



class shasta::AssemblyGraph3EdgeStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;

    AssemblyGraph3EdgeStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraph3EdgeStep()
    {}

};



class shasta::AssemblyGraph3Edge : public vector<AssemblyGraph3EdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraph3Edge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
};



class shasta::AssemblyGraph3 :
    public AssemblyGraph3BaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph3> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph3(
        const Anchors&,
        const AnchorGraph&,
        const AssemblerOptions&);

    // Detangle, phase, assemble sequence, output.
    void run(uint64_t threadCount);

private:
    const Anchors& anchors;
    const AssemblerOptions& assemblerOptions;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;

    void check() const;
};
