#pragma once


/******************************************************************

In the PhasingGraph, each vertex corresponds to
a non-trivial bubble in the SuperbubbleChain being phased.
A directed edge is created between two vertices if
the corresponding bubbles can be phased relative to each other.

******************************************************************/

// Shasta.
#include "invalid.hpp"
#include "GTest.hpp"

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

    // Length of longest path from a source vertex to this vertex.
    // Used when computing longest paths.
    uint64_t pathLength = invalid<uint64_t>;


    PhasingGraphVertex(uint64_t position) : position(position) {}
};



class shasta::PhasingGraphEdge {
public:
    bool isShortestPathEdge = false;
    GTest::Hypothesis bestHypothesis;

    PhasingGraphEdge(const GTest::Hypothesis& bestHypothesis) :
        bestHypothesis(bestHypothesis) {}
};



class shasta::PhasingGraph: public PhasingGraphBaseClass {
public:

    void addVertex(uint64_t position);
    void addEdge(
        uint64_t position0,
        uint64_t position1,
        const GTest::Hypothesis& bestHypothesis);

    // Remove low degree vertices and return the number of such vertices that were removed.
    uint64_t removeLowDegreeVertices(uint64_t minDegree);

    // Compute connected components consisting of at least two vertices.
    // Each connected component is a sorted vector of positions in the SuperbubbleChain.
    // They are returned sorted by decreasing size.
    vector< vector<uint64_t> > components;
    void computeConnectedComponents();

    // Find the longest path in each connected component.
    vector< vector<edge_descriptor> > longestPaths;
    void findLongestPaths();

    void writeGraphviz(const string& fileName) const;

private:
    // The vertexTable contains the vertex descriptor at
    // each position in the SuperbubbleChain.
    // If a position does not correspond to a vertex,
    // it contains null_vertex().
    vector<vertex_descriptor> vertexTable;
};
