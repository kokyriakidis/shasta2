#pragma once

// A RestrictedAnchorGraph is a small AnchorGraph constructed
// using only selected portions of the Journeys of a given set of OrientedReadIds.

// Shasta
#include "AnchorPair.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <utility.hpp>
#include <vector.hpp>



namespace shasta {

    class RestrictedAnchorGraph;
    class RestrictedAnchorGraphVertex;
    class RestrictedAnchorGraphEdge;

    using RestrictedAnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        RestrictedAnchorGraphVertex,
        RestrictedAnchorGraphEdge>;
    class RestrictedAnchorGraph;

    class JourneyPortion;

    class TangleMatrix1;
}



class shasta::JourneyPortion {
public:
    OrientedReadId orientedReadId;
    uint32_t begin;
    uint32_t end;
    JourneyPortion(OrientedReadId orientedReadId, uint32_t begin, uint32_t end) :
        orientedReadId(orientedReadId), begin(begin), end(end) {}
};



class shasta::RestrictedAnchorGraphVertex {
public:
    AnchorId anchorId;
    RestrictedAnchorGraphVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}

    // Fields used by approximateTopologicalSort.
    uint64_t color;
    uint64_t rank;
};



class shasta::RestrictedAnchorGraphEdge {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;

    bool isOptimalPathEdge = false;

    // Field used by approximateTopologicalSort.
    bool isDagEdge = false;
};



class shasta::RestrictedAnchorGraph : public RestrictedAnchorGraphBaseClass {
public:
    RestrictedAnchorGraph(
        const Anchors&,
        const Journeys&,
        const TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);

    void create(
        const Anchors&,
        const Journeys&,
        const vector<JourneyPortion>& journeyPortions,
        ostream& html);

    std::map<AnchorId, vertex_descriptor> vertexMap;

    // Return the vertex_descriptor correspponding to an AnchorId,
    // creating the vertex if necessary.
    vertex_descriptor getVertex(AnchorId);

    // Only keep vertices that are forward reachable from the
    // vertex at anchorId0 and backward reachable from the vertex at anchorId1.
    void keepBetween(AnchorId anchorId0, AnchorId anchorId1);

    // Remove cycles by doing an approximate topological ordering,
    // the removing edges that are not DAG edges.
    void removeCycles();

    // Find the longest path.
    // This also sets the isOptimalPathEdge on the edges of the longest path.
    void findLongestPath(vector<edge_descriptor>&);

    // Find the optimal assembly path.
    // This also sets the isOptimalPathEdge on the edges of the optimal path path.
    void findOptimalPath(
        AnchorId anchorId0,
        AnchorId anchorId1,
        vector<edge_descriptor>&);

    void writeGraphviz(const string& fileName, const vector<AnchorId>& highlightVertices) const;
    void writeGraphviz(ostream&, const vector<AnchorId>& highlightVertices) const;
};
