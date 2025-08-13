// Shasta.
#include "RestrictedAnchorGraph.hpp"
#include "findReachableVertices.hpp"
#include "Journeys.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;

// Boost libraries.
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
    using OrientedReadInfo = TangleMatrix1::OrientedReadInfo;

    const vector<OrientedReadInfo>& entranceOrientedReadInfos = tangleMatrix1.entranceOrientedReadInfos[iEntrance];
    const vector<OrientedReadInfo>& exitOrientedReadInfos = tangleMatrix1.exitOrientedReadInfos[iExit];



    // Joint loop over the oriented reads of the entrance and exit.
    vector<JourneyPortion> journeyPortions;
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



    // Create the graph using these journey portions.
    create(anchors, journeys, journeyPortions, html);
}



void RestrictedAnchorGraph::create(
    const Anchors& anchors,
    const Journeys& journeys,
    const vector<JourneyPortion>& journeyPortions,
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

    // Loop over the JourneyPortions.
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
        dot << "\"" << anchorIdToString(anchorId) << "\"";
        if(find(highlightVertices.begin(), highlightVertices.end(), anchorId) != highlightVertices.end()) {
            dot << "[style=filled fillcolor=pink]";
        }
        dot << ";\n";
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
            " label=\"" << coverage << "\\n" << offset << "\""
            "];\n";
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
