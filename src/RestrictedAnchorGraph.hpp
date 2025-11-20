#pragma once

// A RestrictedAnchorGraph is a small AnchorGraph constructed
// using only selected portions of the Journeys of a given set of OrientedReadIds.

// Shasta
#include "AnchorPair.hpp"
#include "CycleAvoider.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
#include "SimpleMap.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {

    class RestrictedAnchorGraph;
    class RestrictedAnchorGraphVertex;
    class RestrictedAnchorGraphEdge;

    using RestrictedAnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
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



class shasta::RestrictedAnchorGraphVertex : public CycleAvoiderVertex {
public:
    AnchorId anchorId;
    vector<OrientedReadId> orientedReadIds;
    RestrictedAnchorGraphVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}
};



class shasta::RestrictedAnchorGraphEdge {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;

    bool isOptimalPathEdge = false;

    // Field used by approximateTopologicalSort.
    bool isDagEdge = false;

    // Field used by removeLowCoverageEdges.
    bool wasRemoved = false;
};



class shasta::RestrictedAnchorGraph : public RestrictedAnchorGraphBaseClass {
public:

    // Constructor using a TangleMatrix1.
    RestrictedAnchorGraph(
        const Anchors&,
        const Journeys&,
        const TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);
    void constructFromTangleMatrix1(
        const Anchors&,
        const Journeys&,
        const TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);

    // The journey portions that define this RestrictedAnchorGraph.
    vector<JourneyPortion> journeyPortions;

    // Fill the journey portions using a TangleMatrix1.
    void fillJourneyPortions(
        const Journeys&,
        const TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);

    // A SimpleMap containing information about all the N distinct AnchorIds
    // used in this RestrictedAnchorGraph.
    // The anchorIndex of an AnchorId is a number in [0, N) that uniquely
    // identifies each of these distinct AnchorIds.
    // In other words, it is a perfect hash function for these N AnchorIds.
    class AnchorInformation {
    public:
        uint64_t anchorIndex = invalid<uint64_t>;
        vertex_descriptor v = null_vertex();
    };
    SimpleMap<AnchorId, AnchorInformation> anchorMap;
    vector<AnchorId> anchorIds;    // Indexed by anchorIndex.

    void writeAnchorMap() const;
    void writeAnchorIds() const;



    // Create a new vertex and add it to the anchorMap.
    vertex_descriptor addVertex(AnchorId);

    // Return the vertex_descriptor corresponding to an AnchorId.
    // This asserts if there is not such vertex.
    vertex_descriptor getExistingVertex(AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        SHASTA2_ASSERT(p);
        const vertex_descriptor v = p->second.v;
        SHASTA2_ASSERT(v != null_vertex());
        return v;
    }



    // Return the vertex_descriptor corresponding to an AnchorId.
    // This returns null_vertex() there is not such vertex.
    vertex_descriptor getVertex(AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        SHASTA2_ASSERT(p);
        const AnchorInformation& anchorInformation = p->second;
        return anchorInformation.v;
    }



    // Find out if a vertex with the given AnchorId exists.
    bool vertexExists(AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        SHASTA2_ASSERT(p);
        const AnchorInformation& anchorInformation = p->second;
        return anchorInformation.v != null_vertex();
    }



    // A pointer to CycleAvoider, if we are using one.
    // I was not able to get this to compile with a shared_ptr.
    CycleAvoider<RestrictedAnchorGraph>* cycleAvoider = 0;

    // Only keep vertices that are forward reachable from the
    // vertex at anchorId0 and backward reachable from the vertex at anchorId1.
    void keepBetween(AnchorId anchorId0, AnchorId anchorId1);

    // Approximate topological sort.
    void approximateTopologicalSort();

    // Find the optimal assembly path.
    // This also sets the isOptimalPathEdge on the edges of the optimal path path.
    void findOptimalPath(
        AnchorId anchorId0,
        AnchorId anchorId1,
        vector<edge_descriptor>&);

    // Write a table showing which OrientedReadIds are in each vertex.
    // Vertices are written out in rank order.
    void writeOrientedReadsInVertices(ostream& html) const;

    void writeGraphviz(const string& fileName, const vector<AnchorId>& highlightVertices) const;
    void writeGraphviz(ostream&, const vector<AnchorId>& highlightVertices) const;
    void writeHtml(ostream&, const vector<AnchorId>& highlightVertices) const;



    // Gather all the distinct AnchorIds that appear in the JourneyPortions
    // and store them in the anchorMap.
    void gatherAllAnchorIds(const Journeys&);

    uint64_t getAnchorIndex(AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        SHASTA2_ASSERT(p);
        return p->second.anchorIndex;
    }

    // The anchorIndexes for each Anchor of the JourneyPortions.
    vector< vector<uint64_t> > journeyPortionsAnchorIndexes;
    void fillJourneyPortionsAnchorIndexes(const Journeys&);

    // Gather all transitions(anchorIndex0, anchorIndex1) for consecutive
    // anchors in the journey portions. The number of times each
    // transition appears in the journeys is its coverage.
    // Store the transitions in a vector indexed by coverage.
    class Transition {
    public:
        uint64_t anchorIndex0;
        uint64_t anchorIndex1;
    };
    vector< vector<Transition> > transitions;
    void gatherTransitions(ostream& html);

    // Failure modes.
    class NoTransitions : public std::exception {};
    class Unreachable : public std::exception {};

};
