#pragma once

// In the TangleGraph, each vertex represents a tangle
// and each edge represents a BubbleChain.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "Bubble.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"

namespace shasta2 {

    class TangleGraph;
    class TangleGraphEdge;
    class TangleGraphVertex;

    using TangleGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        TangleGraphVertex,
        TangleGraphEdge>;
}



class shasta2::TangleGraphVertex {
public:
    uint64_t id;
    vector<AssemblyGraphBaseClass::vertex_descriptor> assemblyGraphVertices;

    TangleGraphVertex(
        uint64_t id,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& assemblyGraphVertices) :
        id(id),
        assemblyGraphVertices(assemblyGraphVertices) {}

    TangleGraphVertex(
        uint64_t id,
        AssemblyGraphBaseClass::vertex_descriptor v) :
        id(id),
        assemblyGraphVertices(1, v) {}
};



class shasta2::TangleGraphEdge {
public:
    uint64_t id;    // Also the BubbleChain id.
    BubbleChain bubbleChain;

    TangleGraphEdge(uint64_t id, const BubbleChain& bubbleChain) :
        id(id),
        bubbleChain(bubbleChain)
    {}
};



class shasta2::TangleGraph : public TangleGraphBaseClass {
public:

    TangleGraph(
        const AssemblyGraph&,
        const vector<BubbleChain>&,
        const vector< vector<vertex_descriptor> >& tangles
    );

    const AssemblyGraph& assemblyGraph;

    void writeVertices(const string& fileName) const;
    void writeVertices(ostream&) const;
    void writeEdges(const string& fileName) const;
    void writeEdges(ostream&) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
};
