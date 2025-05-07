#pragma once

// In the AssemblyGraph2, each vertex is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of vertices in the TransitionGraph.
// Each Assembly2GraphVertex generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <iosfwd.hpp>
#include <string.hpp>
#include <vector.hpp>

namespace shasta {

    class AssemblyGraph2;
    class AssemblyGraph2Vertex;
    class AssemblyGraph2VertexStep;
    class AssemblyGraph2Edge;

    using AssemblyGraph2BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraph2Vertex,
        AssemblyGraph2Edge>;

    class Anchors;
    class TransitionGraph;
}



class shasta::AssemblyGraph2VertexStep {
public:
    AnchorPair anchorPair;
    uint64_t offset;

    AssemblyGraph2VertexStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}
};



class shasta::AssemblyGraph2Vertex : public vector<AssemblyGraph2VertexStep> {
public:
    uint64_t id = invalid<uint64_t>;
    AssemblyGraph2Vertex(uint64_t id) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
};



class shasta::AssemblyGraph2Edge {
public:
};



class shasta::AssemblyGraph2 :
    public AssemblyGraph2BaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph2> {
public:
    AssemblyGraph2(
        const Anchors&,
        const TransitionGraph&);
    uint64_t nextVertexId = 0;

    void check(const Anchors&) const;
    void check(edge_descriptor) const;

    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
};
