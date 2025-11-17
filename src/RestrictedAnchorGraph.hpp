#pragma once

// A RestrictedAnchorGraph is a small AnchorGraph constructed
// using only selected portions of the Journeys of a given set of OrientedReadIds.

// Shasta
#include "AnchorPair.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

// Standard library.
#include <utility.hpp>
#include <vector.hpp>



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



class shasta::RestrictedAnchorGraphVertex {
public:
    AnchorId anchorId;
    vector<OrientedReadId> orientedReadIds;
    RestrictedAnchorGraphVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}

    // Fields used by approximateTopologicalSort.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    // Field used by removeLowCoverageEdges.
    bool wasRemoved = false;

    RestrictedAnchorGraphBaseClass::vertex_descriptor dominator =
        RestrictedAnchorGraphBaseClass::null_vertex();
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

    // Constructor from an AnchorId.
    RestrictedAnchorGraph(
        const Anchors&,
        const Journeys&,
        AnchorId,
        uint32_t distanceInJourney,
        ostream& html);

    // Constructor using a TangleMatrix1.
    RestrictedAnchorGraph(
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
        uint64_t iExit);

    // Create the graph from the journey portions.
    void create(
        const Anchors&,
        const Journeys&,
        ostream& html);


    std::map<AnchorId, vertex_descriptor> vertexMap;

    // Return the vertex_descriptor corresponding to an AnchorId,
    // creating the vertex if necessary.
    vertex_descriptor getVertex(AnchorId);

    // Return the vertex_descriptor corresponding to an AnchorId.
    // This asserts if there is not such vertex.
    vertex_descriptor getExistingVertex(AnchorId) const;

    // Only keep vertices that are forward reachable from the
    // vertex at anchorId0 and backward reachable from the vertex at anchorId1.
    void keepBetween(AnchorId anchorId0, AnchorId anchorId1);

    // Approximate topological sort.
    void approximateTopologicalSort();

    // Remove cycles by doing an approximate topological ordering,
    // the removing edges that are not DAG edges.
    void removeCycles();

    // Find the longest path.
    // This also sets the isOptimalPathEdge on the edges of the longest path.
    void findLongestPath(vector<edge_descriptor>&);

    // Remove low coverage edges without destroying reachability
    // of anchorId1 from anchorId0. In the process, this also
    // removes vertices that become unreachable from anchorId0 (forward)
    // or anchorId1 (backward).
    void removeLowCoverageEdges(AnchorId anchorId0, AnchorId anchorId1);

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


    // A filtered graph that only includes vertices and edges for which the wasRemoved flag is not set.
    class FilteringPredicate {
    public:
        bool operator()(const vertex_descriptor& v) const
        {
            return  not (*restrictedAnchorGraph)[v].wasRemoved;
        }
        bool operator()(const edge_descriptor& e) const
        {
            return  not (*restrictedAnchorGraph)[e].wasRemoved;
        }
        FilteringPredicate(const RestrictedAnchorGraph* restrictedAnchorGraph = 0) :
            restrictedAnchorGraph(restrictedAnchorGraph)
            {}
        const RestrictedAnchorGraph* restrictedAnchorGraph;
    };
    FilteringPredicate filteringPredicate;
    using FilteredGraph = boost::filtered_graph<RestrictedAnchorGraph, FilteringPredicate, FilteringPredicate>;
    FilteredGraph filteredGraph;

    // Compute the dominator tree starting at a given vertex.
    // Store the immediate dominator of each vertex in
    // RestrictedAnchorGraphVertex::immediateDominator.
    void computeDominatorTree(vertex_descriptor);
};
