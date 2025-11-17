// Shasta.
#include "RestrictedAnchorGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "dominatorTree.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "longestPath.hpp"
#include "orderPairs.hpp"
#include "TangleMatrix1.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>



RestrictedAnchorGraph::RestrictedAnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html) :
    filteringPredicate(this),
    filteredGraph(*this, filteringPredicate, filteringPredicate)
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
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Journey journey = journeys[orientedReadId];
                const uint32_t begin = itEntrance->positionInJourney;
                const uint32_t end = uint32_t(journey.size());
                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Entrance only " << orientedReadId << endl;
            }
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
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Journey journey = journeys[orientedReadId];
                const uint32_t begin = 0;
                const uint32_t end = itExit->positionInJourney + 1;
                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Exit only " << orientedReadId << endl;
            }
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
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Journey journey = journeys[orientedReadId];
                const uint32_t begin = itEntrance->positionInJourney;
                const uint32_t end = itExit->positionInJourney + 1;
                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Both " << orientedReadId << endl;
            }
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
            const Journey journey = journeys[journeyPortion.orientedReadId];
            html << "<tr>" <<
                "<td class=centered>" << journeyPortion.orientedReadId <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << journeyPortion.begin <<
                "<td class=centered>" << journeyPortion.end <<
                "<td class=centered>" << journeyPortion.end -journeyPortion.begin;
#if 0
            html <<    "<td class=centered>";
            for(uint64_t position=journeyPortion.begin; position!=journeyPortion.end; ++position) {
                if(position != journeyPortion.begin) {
                    html << " ";
                }
                html << anchorIdToString(journey[position]);
            }
#endif
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



void RestrictedAnchorGraph::writeHtml(
    ostream& html,
    const vector<AnchorId>& highlightVertices) const
{
    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, highlightVertices);


    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);
}



// Only keep vertices that are forward reachable from the
// vertex at anchorId0 and backward reachable from the vertex at anchorId1.
// To permit future optimizations, we don't really remove
// the vertices - we only disconnect them for the rest of the graph.
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
    // To permit future optimizations, we don't really remove
    // the vertices - we only disconnect them for the rest of the graph.
    findReachableVertices(graph, v0, 0, reachableVertices);
    SHASTA_ASSERT(reachableVertices.contains(v1));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        // vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        // boost::remove_vertex(v, graph);
    }
    reachableVertices.clear();
    verticesToBeRemoved.clear();

    // Remove vertices that are not backward reachable from v1.
    // To permit future optimizations, we don't really remove
    // the vertices - we only disconnect them for the rest of the graph.
    findReachableVertices(graph, v1, 1, reachableVertices);
    SHASTA_ASSERT(reachableVertices.contains(v0));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        // vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        // boost::remove_vertex(v, graph);
    }
}




void RestrictedAnchorGraph::approximateTopologicalSort()
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
    shasta::approximateTopologicalSort(graph, edgesSortedByDecreasingCoverage);
}




// Remove cycles by doing an approximate topological ordering,
// the removing edges that are not DAG edges.
void RestrictedAnchorGraph::removeCycles()
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    approximateTopologicalSort();

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


#if 0
// Find the optimal assembly path.
// This finds the shortest path from anchorId0 to anchorId1,
// with length of an edge defined to estimate average
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
#endif



// Find the optimal assembly path.
void RestrictedAnchorGraph::findOptimalPath(
    AnchorId anchorId0,
    AnchorId anchorId1,
    vector<edge_descriptor>& optimalPath)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    // Find the vertices corresponding to anchorId0 and anchorId1.
    const auto it0 = vertexMap.find(anchorId0);
    SHASTA_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const auto it1 = vertexMap.find(anchorId1);
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    // Compute the longest path.
    shasta::longestPath(graph, optimalPath);
    SHASTA_ASSERT(not optimalPath.empty());
    SHASTA_ASSERT(source(optimalPath.front(), graph) == v0);
    SHASTA_ASSERT(target(optimalPath.back(), graph) == v1);

    // Set the isOptimalPathEdge flags.
    for(const edge_descriptor e: optimalPath) {
        graph[e].isOptimalPathEdge = true;
    }


}



// Remove low coverage edges without destroying reachability
// of anchorId1 from anchorId0. In the process, this also
// removes vertices that become unreachable from anchorId0 (forward)
// or anchorId1 (backward).
// To permit future optimizations, we don't really remove
// the vertices - we only disconnect them for the rest of the graph.
void RestrictedAnchorGraph::removeLowCoverageEdges(
    AnchorId anchorId0,
    AnchorId anchorId1)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    // Find the vertices corresponding to anchorId0 and anchorId1.
    const auto it0 = vertexMap.find(anchorId0);
    SHASTA_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const auto it1 = vertexMap.find(anchorId1);
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    // First, clear all "wasRemoved" flags.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        graph[v].wasRemoved = false;
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        graph[e].wasRemoved = false;
    }



    // Loop over increasing values of minCoverage.
    for(uint64_t minCoverage=1; ; ++minCoverage) {

        // Tentatively flag as removed the edges with coverage less than minCoverage.
        // Keep track of them.
        vector<edge_descriptor> edgesTentativelyRemoved;
        BGL_FORALL_EDGES(e, graph, Graph) {
            RestrictedAnchorGraphEdge& edge = graph[e];
            if(edge.anchorPair.size() < minCoverage) {
                edge.wasRemoved = true;
                edgesTentativelyRemoved.push_back(e);
            }
        }



        // If v1 is still reachable from v0 (on the filteredGraph),
        // remove from the RestrictedAnchorGraph all edges
        // we tentatively marked as removed, plus the vertices
        // that are no longer reachable from anchorId0 (forward)
        // and anchorId1 (backward).
        // Then try increasing minCoverage more.
        if(isReachable(filteredGraph, v0, v1, 0)) {

            for(const edge_descriptor e: edgesTentativelyRemoved) {
                boost::remove_edge(e, graph);
            }

            // Find vertices that are still forward reachable from v0.
            std::set<vertex_descriptor> reachableVertices0;
            findReachableVertices(graph, v0, 0, reachableVertices0);

            // Find vertices that are still backward reachable from v1.
            std::set<vertex_descriptor> reachableVertices1;
            findReachableVertices(graph, v1, 1, reachableVertices1);

            // Remove vertices that are no longer reachable in both directions.
            vector<vertex_descriptor> verticesToBeRemoved;
            BGL_FORALL_VERTICES(v, graph, Graph) {
                if(not (reachableVertices0.contains(v) and reachableVertices1.contains(v))) {
                    SHASTA_ASSERT(v != v0);
                    SHASTA_ASSERT(v != v1);
                    verticesToBeRemoved.push_back(v);
                }
            }
            for(const vertex_descriptor v: verticesToBeRemoved) {
                // vertexMap.erase(graph[v].anchorId);
                boost::clear_vertex(v, graph);
                // boost::remove_vertex(v, graph);
            }
        }

        else {
            // If v1 is no longer reachable from v0 (on the filteredGraph),
            // put back the edgesTentativelyRemoved in the filtered graph
            // and exit the loop over minCoverage.
            for(const edge_descriptor e: edgesTentativelyRemoved) {
                graph[e].wasRemoved = false;
            }
            break;
        }
    }


#if 0

    // Compute the dominator tree starting at v0.
    computeDominatorTree(v0);
#if 0
    cout << "Dominator tree:" << endl;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const vertex_descriptor dominator = graph[v].dominator;
        if(dominator != null_vertex()) {
            cout << anchorIdToString(graph[v].anchorId) << " " <<
                anchorIdToString(graph[dominator].anchorId) << endl;
        }
    }
#endif


    // Walk back the dominator tree starting at v1.
    vector<vertex_descriptor> dominatorSequence({v1});
    while(true) {
        const vertex_descriptor v = dominatorSequence.back();
        const vertex_descriptor dominator = graph[v].dominator;
        if(dominator == null_vertex()) {
            break;
        } else {
            dominatorSequence.push_back(dominator);
        }
    }
    std::ranges::reverse(dominatorSequence);
#if 0
    cout << "Dominator sequence:";
    for(const vertex_descriptor v: dominatorSequence) {
        cout << " " << anchorIdToString(graph[v].anchorId);
    }
    cout << endl;
#endif
    SHASTA_ASSERT(dominatorSequence.size() >= 2);
    SHASTA_ASSERT(dominatorSequence.front() == v0);
    SHASTA_ASSERT(dominatorSequence.back() == v1);



    // Each pair of adjacent vertices (vA, vB) in the dominatorSequence
    // define a "segment" in the RestrictedAnchorGraph such that all paths
    // starting at v0 that go through vB must also go through vA first.
    // So vA and vB are "choke points" of the RestrictedAnchorGraph
    // and define a "segment" of the RestrictedAnchorGraph.
    for(uint64_t i=1; i<dominatorSequence.size(); i++) {
        const vertex_descriptor vA = dominatorSequence[i-1];
        const vertex_descriptor vB = dominatorSequence[i];

        // Skip trivial segments.
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(vA, vB, graph);
        if(edgeExists and (out_degree(vA, graph) == 1) and (in_degree(vB, graph) == 1)) {
            continue;
        }

        cout << "Working on non-trivial RestrictedAnchorGraph segment " <<
            anchorIdToString(graph[vA].anchorId) << " " <<
            anchorIdToString(graph[vB].anchorId) << endl;
    }
#endif
}



// Compute the dominator tree starting at a given vertex.
// Store the immediate dominator of each vertex in
// RestrictedAnchorGraphVertex::immediateDominator.
void RestrictedAnchorGraph::computeDominatorTree(vertex_descriptor v0)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    shasta::lengauer_tarjan_dominator_tree_general(graph, v0);
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



// Constructor from an AnchorId.
RestrictedAnchorGraph::RestrictedAnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    AnchorId anchorId,
    uint32_t distanceInJourney,
    ostream& html) :
    filteringPredicate(this),
    filteredGraph(*this, filteringPredicate, filteringPredicate)
{

    // Fill in the JourneyPortions by looking around this Anchor.
    const Anchor anchor = anchors[anchorId];
    for(const AnchorMarkerInfo& markerInfo: anchor) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const uint32_t positionInJourney = markerInfo.positionInJourney;
        const Journey journey = journeys[orientedReadId];

        // Find the journey portion of this OrientedReadId.
        const uint32_t begin = (distanceInJourney >= positionInJourney) ? 0 : (positionInJourney - distanceInJourney);
        uint32_t end = min(uint32_t(journey.size()), positionInJourney + distanceInJourney);
        journeyPortions.emplace_back(orientedReadId, begin, end);
    }



    // Create the RestrictedAnchorGraph using these JourneyPortions.
    create(anchors, journeys, html);
}

