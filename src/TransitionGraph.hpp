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
    AnchorPair anchorPair;
    uint64_t id;
    TransitionGraphVertex(const AnchorPair& anchorPair, uint64_t id) :
        anchorPair(anchorPair), id(id) {}
    TransitionGraphVertex() : id(invalid<uint64_t>) {}
};



// The TransitionGraphEdge contains an AnchorPair to bridge
// between its source and target vertices.
// If its source and target vertices have adjacent AnchorPairs,
// this stores a null AnchorPair (converts to bool false).
class shasta::TransitionGraphEdge {
public:
};



class shasta::TransitionGraph : public TransitionGraphBaseClass {
public:

    TransitionGraph(
        const Anchors&,
        const AnchorGraph&);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // The journey of an oriented read in the TransitionGraph is
    // the sequence of vertices it encounters.
    vector< vector<vertex_descriptor> > journeys;
    void computeJourneys(const Anchors&);

    void check() const;

};
