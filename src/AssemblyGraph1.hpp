#pragma once

// In the AssemblyGraph1, each edge is a linear chain of AnchorPairs.
// In contrast, in the AssemblyGraph each edge is a linear chain of Anchors.
// The AssemblyGraph1 is created using the TransitionGraph.
// In contrast, the AssemblyGraph is created using the AnchorGraph.

// Shasta.
#include "invalid.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <iosfwd.hpp>
#include <string.hpp>
#include <vector.hpp>

namespace shasta {

    class AssemblyGraph1;
    class AssemblyGraph1Vertex;
    class AssemblyGraph1Edge;

    using AssemblyGraph1BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraph1Vertex,
        AssemblyGraph1Edge>;

    class AnchorGraph;
    class AnchorPair;
    class TransitionGraph;
}


class shasta::AssemblyGraph1Vertex {
public:
    uint64_t vertexId = invalid<uint64_t>;
    AssemblyGraph1Vertex(uint64_t vertexId) : vertexId(vertexId) {}
};



class shasta::AssemblyGraph1Edge : public vector<AnchorPair> {
public:
    uint64_t edgeId = invalid<uint64_t>;
    AssemblyGraph1Edge(uint64_t edgeId) : edgeId(edgeId) {}
};



class shasta::AssemblyGraph1 : public AssemblyGraph1BaseClass {
public:
    AssemblyGraph1(
        const AnchorGraph&,
        const TransitionGraph&);
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;

    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
};
