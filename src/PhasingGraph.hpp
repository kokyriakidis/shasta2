#pragma once


/******************************************************************

In the PhasingGraph, each vertex corresponds to
a non-trivial bubble in the SuperbubbleChain being phased.
A directed edge is created between two vertices if
the corresponding bubbles can be phased relative to each other.

******************************************************************/

// Shasta.
#include "invalid.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <string.hpp>
#include <vector.hpp>



namespace shasta {

    class PhasingGraph;
    class PhasingGraphVertex;
    class PhasingGraphEdge;

    using PhasingGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        PhasingGraphVertex,
        PhasingGraphEdge>;
}



class shasta::PhasingGraphVertex {
public:

    // The position in the SuperbubbleChain of the bubble
    // corresponding to this vertex.
    uint64_t position;

    // The component this vertex belongs to.
    uint64_t componentId = invalid<uint64_t>;

    PhasingGraphVertex(uint64_t position) : position(position) {}
};



class shasta::PhasingGraphEdge {
public:
};



class shasta::PhasingGraph: public PhasingGraphBaseClass {
public:

    void addVertex(uint64_t position);
    void addEdge(uint64_t position0, uint64_t position1);

    // Remove isolated vertices and return the number of such vertices that were removed.
    uint64_t removeIsolatedVertices();

    // Compute connected components consisting of at least two vertices.
    // Each connected component is a vector of positions in the SuperbubbleChain.
    // They are returned sorted by decreasing size.
    vector< vector<uint64_t> > components;
    void computeConnectedComponents();

    void writeGraphviz(const string& fileName) const;

private:
    // The vertexTable contains the vertex descriptor at
    // each position in the SuperbubbleChain.
    // If a position does not correspond to a vertex,
    // it contains null_vertex().
    vector<vertex_descriptor> vertexTable;
};
