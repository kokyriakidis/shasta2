#pragma once

// A RestrictedAnchorGraph is a small AnchorGraph constructed
// using only selected portions of the Journeys of a given set of OrientedReadIds.

// Shasta
#include "AnchorPair.hpp"
#include "CycleAvoider.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"

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


    // Table that contains pairs (AnchorId, vertex_descriptor) for each vertex.
    // We keep it sorted by AnchorId, so getExistingVertex can use a binary search.
    vector< pair<AnchorId, vertex_descriptor> > vertexTable;
    void sortVertexTable();
    bool vertexTableIsValid = true;

    // Create a new vertex and add it to the vertexTable, without resorting
    // the vertexTable. This invalidates the vertexTable.
    // If this is called with the AnchorId of an existing vertex,
    // the call succeeds but the subsequent call to sortVertexTable will assert.
    vertex_descriptor addVertex(AnchorId);

    // Return the vertex_descriptor corresponding to an AnchorId.
    // This asserts if there is not such vertex.
    vertex_descriptor getExistingVertex(AnchorId anchorId) const
    {
        SHASTA2_ASSERT(vertexTableIsValid);

        const auto it = lower_bound(
            vertexTable.begin(), vertexTable.end(),
            make_pair(anchorId, null_vertex()),
            OrderPairsByFirstOnly<AnchorId, vertex_descriptor>());
        SHASTA2_ASSERT(it != vertexTable.end());
        SHASTA2_ASSERT(it->first == anchorId);
        return it->second;
    }



    // Return the vertex_descriptor corresponding to an AnchorId.
    // This returns null_vertex() there is not such vertex.
    vertex_descriptor getVertex(AnchorId anchorId) const
    {
        SHASTA2_ASSERT(vertexTableIsValid);
        const auto it = lower_bound(
            vertexTable.begin(), vertexTable.end(),
            make_pair(anchorId, null_vertex()),
            OrderPairsByFirstOnly<AnchorId, vertex_descriptor>());

        if((it == vertexTable.end()) or (it->first != anchorId)) {
            return null_vertex();
        } else {
            return it->second;
        }
    }



    // Find out if a vertex with the given AnchorId exists.
    bool vertexExists(AnchorId anchorId) const
    {
        SHASTA2_ASSERT(vertexTableIsValid);

        const auto it = lower_bound(
            vertexTable.begin(), vertexTable.end(),
            make_pair(anchorId, null_vertex()),
            OrderPairsByFirstOnly<AnchorId, vertex_descriptor>());

        return (it != vertexTable.end())  and (it->first == anchorId);
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
    // and store them sorted.
    vector<AnchorId> allAnchorIds;
    void gatherAllAnchorIds(const Journeys&);

    // The index of an AnchorId in the allAnchorIds vector is called "anchorIndex"
    // in constructFromTangleMatrix1 code,
    // and serves as a perfect hash function for these AnchorIds.
    uint64_t getAnchorIndex(AnchorId anchorId) const
    {
        const auto it = find(allAnchorIds.begin(), allAnchorIds.end(), anchorId);
        // SHASTA2_ASSERT(it != allAnchorIds.end());
        // SHASTA2_ASSERT(*it == anchorId);
        return it - allAnchorIds.begin();
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
