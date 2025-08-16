// Shasta.
#include "RestrictedAnchorGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "findReachableVertices.hpp"
#include "Journeys.hpp"
#include "longestPath.hpp"
#include "orderPairs.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>



RestrictedAnchorGraph::RestrictedAnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    fillJourneyPortions(journeys, tangleMatrix1, iEntrance, iExit);
    create(anchors, journeys, html);
}



// Fill the journey portions using a TangleMatrix1.
void RestrictedAnchorGraph::fillJourneyPortions(
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit)
{
    using OrientedReadInfo = TangleMatrix1::OrientedReadInfo;

    const vector<OrientedReadInfo>& entranceOrientedReadInfos = tangleMatrix1.entranceOrientedReadInfos[iEntrance];
    const vector<OrientedReadInfo>& exitOrientedReadInfos = tangleMatrix1.exitOrientedReadInfos[iExit];



    // Joint loop over the oriented reads of the entrance and exit.
    auto itEntrance = entranceOrientedReadInfos.begin();
    auto itExit = exitOrientedReadInfos.begin();
    const auto itEntranceEnd = entranceOrientedReadInfos.end();
    const auto itExitEnd = exitOrientedReadInfos.end();
    while(true) {

        // If both iterators are at their end, we are done.
        if((itEntrance==itEntranceEnd) and (itExit==itExitEnd)) {
            break;
        }

        // Case where the OrientedReadId appears only in the entrance.
        // In this case we use a journey portion that begins at the
        // journey position stored in the OrientedReadInfo and
        // ends at the end of the journey.
        if(
            (itEntrance != itEntranceEnd) and
            (itExit==itExitEnd or (itEntrance->orientedReadId < itExit->orientedReadId))) {
            const OrientedReadId orientedReadId = itEntrance->orientedReadId;
            const Journey journey = journeys[orientedReadId];
            const uint32_t begin = itEntrance->positionInJourney;
            const uint32_t end = uint32_t(journey.size());
            journeyPortions.emplace_back(orientedReadId, begin, end);
            // cout << "Entrance only " << orientedReadId << endl;
            ++itEntrance;
        }

        // Case where the OrientedReadId appears only in the exit.
        // In this case we use a journey portion that begins at the
        // beginning of the journey and ends at the
        // journey position stored in the OrientedReadInfo.
        else if(
            (itExit != itExitEnd) and
            (itEntrance==itEntranceEnd or (itExit->orientedReadId < itEntrance->orientedReadId))) {
            const OrientedReadId orientedReadId = itExit->orientedReadId;
            const Journey journey = journeys[orientedReadId];
            const uint32_t begin = 0;
            const uint32_t end = itExit->positionInJourney + 1;
            journeyPortions.emplace_back(orientedReadId, begin, end);
            // cout << "Exit only " << orientedReadId << endl;
            ++itExit;
        }

        // Case where the OrientedReadId appears in both the entrance and the exit.
        // In this case we use a journey portion that begins at the
        // journey position stored in the entrance OrientedReadInfo and
        // ends at the journey position stored in the exit OrientedReadInfo.
        else if(
            (itEntrance != itEntranceEnd) and
            (itExit!=itExitEnd) and
            (itEntrance->orientedReadId == itExit->orientedReadId)) {
            const OrientedReadId orientedReadId = itEntrance->orientedReadId;
            const Journey journey = journeys[orientedReadId];
            const uint32_t begin = itEntrance->positionInJourney;
            const uint32_t end = itExit->positionInJourney + 1;
            journeyPortions.emplace_back(orientedReadId, begin, end);
            // cout << "Both " << orientedReadId << endl;
            ++itEntrance;
            ++itExit;
        }

        // All cases have been covered.
        else {
            SHASTA_ASSERT(0);
        }
    }
}



// Create the graph from the journey portions.
void RestrictedAnchorGraph::create(
    const Anchors& anchors,
    const Journeys& journeys,
    ostream& html)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    if(html) {
        html << "<p>The RestrictedAnchorGraph uses the following journey portions."
            "<table><tr>"
            "<th>Oriented<br>read<br>id"
            "<th>Journey<br>length"
            "<th>Begin<th>End"
            "<th>Journey<br>portion<br>length";
        for(const JourneyPortion& journeyPortion: journeyPortions) {
            html << "<tr>" <<
                "<td class=centered>" << journeyPortion.orientedReadId <<
                "<td class=centered>" << journeys[journeyPortion.orientedReadId].size() <<
                "<td class=centered>" << journeyPortion.begin <<
                "<td class=centered>" << journeyPortion.end <<
                "<td class=centered>" << journeyPortion.end -journeyPortion.begin;
        }
        html << "</table>";
    }



    // Loop over the JourneyPortions to create the vertices.
    for(const JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i=begin; i<end; i++) {
            const AnchorId anchorId = journey[i];
            const vertex_descriptor v = getVertex(anchorId);
            graph[v].orientedReadIds.push_back(orientedReadId);
        }
    }



    // Loop over the JourneyPortions to create the edges.
    for(const JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i1=begin+1; i1<end; i1++) {
            const uint32_t i0 = i1 - 1;
            const AnchorId anchorId0 = journey[i0];
            const AnchorId anchorId1 = journey[i1];
            const vertex_descriptor v0 = getVertex(anchorId0);
            const vertex_descriptor v1 = getVertex(anchorId1);

            // See if we already have this edge.
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, graph);
            if(not edgeExists)  {
                tie(e, ignore) = add_edge(v0, v1, *this);
                auto& edge = graph[e];
                edge.anchorPair.anchorIdA = anchorId0;
                edge.anchorPair.anchorIdB = anchorId1;
            }
            auto& edge = graph[e];
            edge.anchorPair.orientedReadIds.push_back(orientedReadId);
        }
    }

    // Fill in the offsets.
    BGL_FORALL_EDGES(e, graph, Graph) {
        RestrictedAnchorGraphEdge& edge = graph[e];
        edge.offset = edge.anchorPair.getAverageOffset(anchors);
    }
}



RestrictedAnchorGraph::vertex_descriptor RestrictedAnchorGraph::getVertex(AnchorId anchorId)
{
    const auto it = vertexMap.find(anchorId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(RestrictedAnchorGraphVertex(anchorId), *this);
        vertexMap.insert({anchorId, v});
        return v;
    } else {
        return it->second;
    }
}



void RestrictedAnchorGraph::writeGraphviz(
    const string& fileName,
    const vector<AnchorId>& highlightVertices) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, highlightVertices);
}



void RestrictedAnchorGraph::writeGraphviz(
    ostream& dot,
    const vector<AnchorId>& highlightVertices) const
{
    using Graph = RestrictedAnchorGraph;
    const Graph& graph = *this;

    dot << "digraph RestrictedAnchorGraph {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const AnchorId anchorId = graph[v].anchorId;
        dot <<
            "\"" << anchorIdToString(anchorId) << "\""
            " [label=\"" << anchorIdToString(anchorId) << "\\n" << graph[v].orientedReadIds.size() << "\"";
        if(find(highlightVertices.begin(), highlightVertices.end(), anchorId) != highlightVertices.end()) {
            dot << " style=filled fillcolor=pink";
        }
        dot << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AnchorId anchorId0 = graph[v0].anchorId;
        const AnchorId anchorId1 = graph[v1].anchorId;
        const RestrictedAnchorGraphEdge& edge = graph[e];
        const uint64_t coverage = edge.anchorPair.size();
        const uint64_t offset = edge.offset;

        dot << "\"" << anchorIdToString(anchorId0) << "\"->\"" <<
            anchorIdToString(anchorId1) << "\""
            "["
            "penwidth=" << std::setprecision(2) << 0.5 * double(coverage) <<
            " label=\"" << coverage << "\\n" << offset << "\"";

        if(edge.isOptimalPathEdge) {
            dot << " color=green";
        }

        dot << "];\n";
    }

    dot << "}\n";
}



// Only keep vertices that are forward reachable from the
// vertex at anchorId0 and backward reachable from the vertex at anchorId1.
void RestrictedAnchorGraph::keepBetween(AnchorId anchorId0, AnchorId anchorId1)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    const auto it0 = vertexMap.find(anchorId0);
    SHASTA_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;

    const auto it1 = vertexMap.find(anchorId1);
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    std::set<vertex_descriptor> reachableVertices;
    vector<vertex_descriptor> verticesToBeRemoved;

    // Remove vertices that are not forward reachable from v0.
    findReachableVertices(graph, v0, 0, reachableVertices);
    SHASTA_ASSERT(reachableVertices.contains(v1));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }
    reachableVertices.clear();
    verticesToBeRemoved.clear();

    // Remove vertices that are not backward reachable from v1.
    findReachableVertices(graph, v1, 1, reachableVertices);
    SHASTA_ASSERT(reachableVertices.contains(v0));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }
}



// Remove cycles by doing an approximate topological ordering,
// the removing edges that are not DAG edges.
void RestrictedAnchorGraph::removeCycles()
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    // Gather the edges with their coverage.
    vector< pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, Graph)
    {
        const uint64_t coverage = graph[e].anchorPair.size();
        edgeTable.push_back(make_pair(e, coverage));
    }

    // Sort by decreasing coverage.
    std::ranges::sort(edgeTable, OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
    vector<edge_descriptor> edgesSortedByDecreasingCoverage;
    for(const auto& p: edgeTable) {
        edgesSortedByDecreasingCoverage.push_back(p.first);
    }

    // Do the approxinate topological sort.
    approximateTopologicalSort(graph, edgesSortedByDecreasingCoverage);

    // Gather the edges to be removed.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph)
    {
        if(not graph[e].isDagEdge) {
            edgesToBeRemoved.push_back(e);
        }
    }

    // Remove these edges.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }


}



// Find the longest path.
void RestrictedAnchorGraph::findLongestPath(vector<edge_descriptor>& longestPath)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    shasta::longestPath(graph, longestPath);

    for(const edge_descriptor e: longestPath) {
        graph[e].isOptimalPathEdge = true;
    }
}



// Find the optimal assembly path.
// This finds the shortest path from anchorId0 to anchorId1,
// with lengh of an edge defined to estimate average
// number of assembly errors in the edge.
void RestrictedAnchorGraph::findOptimalPath(
    AnchorId anchorId0,
    AnchorId anchorId1,
    vector<edge_descriptor>& optimalPath)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    // Find the vertices corresponding to anchorId0 anc anchorId1.
    const auto it0 = vertexMap.find(anchorId0);
    SHASTA_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const auto it1 = vertexMap.find(anchorId1);
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Define the weight map, which estimates the "cost" of each edge
    // (number of likely errors).
    const double epsilon = 0.2;
    std::map<edge_descriptor, double> weightMap;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const RestrictedAnchorGraphEdge& edge = graph[e];
        const double coverage = double(edge.anchorPair.size());
        const double offset = double(edge.offset);
        const double w = offset * std::pow(epsilon, coverage);
        weightMap.insert(make_pair(e, w));
    }

    // The predecessor map will record the predecessor of each vertex
    // in the shortest paths.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

    // Compute shortest paths starting at v0.
    using boost::dijkstra_shortest_paths;
    using boost::vertex_index_map;
    using boost::make_assoc_property_map;
    using boost::weight_map;
    using boost::predecessor_map;
    dijkstra_shortest_paths(graph, v0,
        vertex_index_map(make_assoc_property_map(vertexIndexMap)).
        weight_map(make_assoc_property_map(weightMap)).
        predecessor_map(make_assoc_property_map(predecessorMap))
    );

    // Now we can use the predecessor map to assemble the optimal path from v0 to v1.
    optimalPath.clear();
    vertex_descriptor v = v1;
    while(v != v0) {
        const auto it = predecessorMap.find(v);
        SHASTA_ASSERT(it != predecessorMap.end());
        const vertex_descriptor vPrevious = it->second;
        edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(vPrevious, v, graph);
        SHASTA_ASSERT(edgeWasFound);
        optimalPath.push_back(e);
        graph[e].isOptimalPathEdge = true;
        v = vPrevious;
    }
    std::ranges::reverse(optimalPath);
}



// Write a table showing which OrientedReadIds are in each vertex.
// Vertices are written out in rank order.
void RestrictedAnchorGraph::writeOrientedReadsInVertices(ostream& html) const
{
    using Graph = RestrictedAnchorGraph;
    const Graph& graph = *this;

    // Gather the vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > verticesSortedByRank;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        verticesSortedByRank.push_back(make_pair(v, graph[v].rank));
    }
    std::ranges::sort(verticesSortedByRank, OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());

    // Write the table header.
    html << "<table><tr><th>Oriented<br>read<br<id>";
    for(const auto& p: verticesSortedByRank) {
        const vertex_descriptor v = p.first;
        const AnchorId anchorId = graph[v].anchorId;
        html << "<th>" << anchorIdToString(anchorId);
    }



    // Write one row for each OrientedReadId.
    for(const JourneyPortion& journeyPortion: journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;

        html << "<tr><th>" << orientedReadId;

        for(const auto& p: verticesSortedByRank) {
            const vertex_descriptor v = p.first;
            const bool vertexContainsOrientedReadId =
                std::ranges::binary_search(graph[v].orientedReadIds, orientedReadId);

            html << "<td";
            if(vertexContainsOrientedReadId) {
                html << " style='background-color:LightGreen'";
            }
            html << ">" << int(vertexContainsOrientedReadId);
        }
    }


    html << "</table>";
}
