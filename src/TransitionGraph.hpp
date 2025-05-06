#pragma once

// In the TransitionGraph, each vertex represents an edge
// of the AncorGraph and its AnchorPair.
// An edge v0->v1 is created if the AnchorGraph edge journeys
// of a sufficient number of OrientedReadIds visit
// v1 immediately after visiting v0.

// Shasta2.
#include "AnchorGraph.hpp"
#include "invalid.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include <string.hpp>



namespace shasta {
    class TransitionGraph;
    class TransitionGraphVertex;
    class TransitionGraphEdge;

    using TransitionGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        TransitionGraphVertex,
        TransitionGraphEdge>;
}



class shasta::TransitionGraphVertex {
public:
    AnchorGraph::edge_descriptor eAnchorGraph;
    uint64_t vertexId;
    TransitionGraphVertex(AnchorGraph::edge_descriptor eAnchorGraph, uint64_t vertexId) :
        eAnchorGraph(eAnchorGraph), vertexId(vertexId) {}
    TransitionGraphVertex() : vertexId(invalid<uint64_t>) {}
};



class shasta::TransitionGraphEdge {
public:
    uint64_t coverage;
    TransitionGraphEdge(uint64_t coverage) : coverage(coverage) {}
};



class shasta::TransitionGraph : public TransitionGraphBaseClass {
public:

    TransitionGraph(
        const Anchors&,
        const AnchorGraph&);

    // A map that gives the TransitionGraph vertex
    // corresponding to a given AnchorGraph edge.
    std::map<AnchorGraph::edge_descriptor, vertex_descriptor> vertexMap;

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

};
