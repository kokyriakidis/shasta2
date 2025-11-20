// Shasta.
#include "RestrictedAnchorGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "longestPath.hpp"
#include "Markers.hpp"
#include "TangleMatrix1.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



RestrictedAnchorGraph::RestrictedAnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    constructFromTangleMatrix1(anchors, journeys, tangleMatrix1, iEntrance, iExit, html);
}



void RestrictedAnchorGraph::constructFromTangleMatrix1(
    const Anchors& anchors,
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    fillJourneyPortions(journeys, tangleMatrix1, iEntrance, iExit, html);
    gatherAllAnchorIds(journeys);
    fillJourneyPortionsAnchorIndexes(journeys);
    gatherTransitions(html);

    // Create the CycleAvoider.
    cycleAvoider = new CycleAvoider<RestrictedAnchorGraph>(*this);


    // A vector<bool>, indexed by anchorIndex, to contain anchorIndexes
    // for vertices we already added.
    vector<bool> verticesAdded(allAnchorIds.size(), false);

    // Add the entrance and exit vertices.
    const AssemblyGraph::edge_descriptor entrance = tangleMatrix1.entrances[iEntrance];
    const AssemblyGraph::edge_descriptor exit = tangleMatrix1.exits[iExit];
    const AnchorId anchorId0 = tangleMatrix1.assemblyGraph[entrance].back().anchorPair.anchorIdB;
    const AnchorId anchorId1 = tangleMatrix1.assemblyGraph[exit].front().anchorPair.anchorIdA;
    SHASTA2_ASSERT(anchorId0 != anchorId1);
    const vertex_descriptor v0 = addVertex(anchorId0);
    const vertex_descriptor v1 = addVertex(anchorId1);
    sortVertexTable();
    const uint64_t anchorIndex0 = getAnchorIndex(anchorId0);
    const uint64_t anchorIndex1 = getAnchorIndex(anchorId1);
    verticesAdded[anchorIndex0] = true;
    verticesAdded[anchorIndex1] = true;

    /*
    cout << tangleMatrix1.assemblyGraph[entrance].id << " " << tangleMatrix1.assemblyGraph[exit].id <<
        " " << allAnchorIds.size() << endl;
    */
    if(html) {
        html << "<br>Total number of distinct AnchorIds on the JourneyPortions is " << allAnchorIds.size();
    }

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(allAnchorIds.size());
    vector<uint64_t> parent(allAnchorIds.size());
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<allAnchorIds.size(); i++) {
        disjointSets.make_set(i);
    }



    // Add edges in order of decreasing coverage.
    vector<uint64_t> anchorIndexesToAdd;
    for(uint64_t coverage=transitions.size()-1; coverage>0; coverage--) {
        const vector<Transition> transitionsToAdd = transitions[coverage];

        // Gather the anchorIndexes for the vertices that we need to add.
        anchorIndexesToAdd.clear();
        for(const Transition& transition: transitionsToAdd) {
            const uint64_t anchorIndex0 = transition.anchorIndex0;
            if(not verticesAdded[anchorIndex0]) {
                anchorIndexesToAdd.push_back(anchorIndex0);
                verticesAdded[anchorIndex0] = true;
            }
            const uint64_t anchorIndex1 = transition.anchorIndex1;
            if(not verticesAdded[anchorIndex1]) {
                anchorIndexesToAdd.push_back(anchorIndex1);
                verticesAdded[anchorIndex1] = true;
            }
        }

        // Add these vertices.
        for(const uint64_t anchorIndex: anchorIndexesToAdd) {
            const AnchorId anchorId = allAnchorIds[anchorIndex];
            addVertex(anchorId);
        }
        sortVertexTable();

        // Now we can add an edge for each of these transitions.
        for(const Transition& transition: transitionsToAdd) {
            const uint64_t anchorIndex0 = transition.anchorIndex0;
            const uint64_t anchorIndex1 = transition.anchorIndex1;
            const AnchorId anchorId0 = allAnchorIds[anchorIndex0];
            const AnchorId anchorId1 = allAnchorIds[anchorIndex1];
            const vertex_descriptor v0 = getExistingVertex(anchorId0);
            const vertex_descriptor v1 = getExistingVertex(anchorId1);

            edge_descriptor e;
            bool edgeWasAdded = false;
            tie(e, edgeWasAdded) = cycleAvoider->addEdge(v0, v1);
            if(edgeWasAdded) {
                RestrictedAnchorGraphEdge& edge = graph[e];
                edge.anchorPair.anchorIdA = anchorId0;
                edge.anchorPair.anchorIdB = anchorId1;
                disjointSets.union_set(anchorIndex0, anchorIndex1);
            }
        }

        // If the entrance and exit are different connected components,
        // continue the loop over coverage.
        if(disjointSets.find_set(anchorIndex0) != disjointSets.find_set(anchorIndex1)) {
            continue;
        }

        // If the exit is reachable from the entrance, stop the loop over coverage.
        if(isReachable(graph, v0, v1, 0)) {
            break;
        }

        // If still not reachable and coverage is 1, we can't reduce coverage further.
        if(coverage == 1) {
            delete cycleAvoider;
            cycleAvoider = 0;
            throw Unreachable();
        }
    }



    // Destroy the CycleAvoider.
    delete cycleAvoider;
    cycleAvoider = 0;



    // Now fill in the OrientedReadIds of each vertex and edge.
    for(uint64_t j=0; j<journeyPortions.size(); j++) {
        const JourneyPortion& journeyPortion = journeyPortions[j];
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const vector<uint64_t>& journeyPortionAnchorIndexes = journeyPortionsAnchorIndexes[j];

        // Vertices.
        for(const uint64_t anchorIndex: journeyPortionAnchorIndexes) {
            const uint64_t anchorId = allAnchorIds[anchorIndex];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                graph[v].orientedReadIds.push_back(orientedReadId);
            }
        }

        // Edges.
        for(uint64_t i1=1; i1<journeyPortionAnchorIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t anchorIndex0 = journeyPortionAnchorIndexes[i0];
            const uint64_t anchorIndex1 = journeyPortionAnchorIndexes[i1];
            const uint64_t anchorId0 = allAnchorIds[anchorIndex0];
            const vertex_descriptor v0 = getVertex(anchorId0);
            if(v0 == null_vertex()) {
                continue;
            }
            const uint64_t anchorId1 = allAnchorIds[anchorIndex1];
            const vertex_descriptor v1 = getVertex(anchorId1);
            if(v1 == null_vertex()) {
                continue;
            }
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, graph);
            if(edgeExists) {
                graph[e].anchorPair.orientedReadIds.push_back(orientedReadId);
            }
        }
    }

    // Fill in the offsets.
    BGL_FORALL_EDGES(e, graph, Graph) {
        RestrictedAnchorGraphEdge& edge = graph[e];
        edge.offset = edge.anchorPair.getAverageOffset(anchors);
    }

    // Remove unreachable portions.
    keepBetween(anchorId0, anchorId1);
}



// Gather all the distinct AnchorIds that appear in the JourneyPortions
// and store them sorted.
void RestrictedAnchorGraph::gatherAllAnchorIds(const Journeys& journeys)
{
    allAnchorIds.clear();
    for(const JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i=begin; i<end; i++) {
            const AnchorId anchorId = journey[i];
            allAnchorIds.push_back(anchorId);
        }
    }
    deduplicate(allAnchorIds);
}



void RestrictedAnchorGraph::fillJourneyPortionsAnchorIndexes(const Journeys& journeys)
{
    const uint64_t orientedReadCount = journeyPortions.size();
    journeyPortionsAnchorIndexes.clear();
    journeyPortionsAnchorIndexes.resize(orientedReadCount);

    // Loop over the oriented reads.
    for(uint64_t j=0; j<orientedReadCount; j++) {
        const JourneyPortion& journeyPortion = journeyPortions[j];
        vector<uint64_t>& journeyPortionAnchorIndexes =journeyPortionsAnchorIndexes[j];

        // Loop over the journey portion of this oriented read.
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        journeyPortionAnchorIndexes.reserve(end - begin);
        for(uint32_t i=begin; i<end; i++) {
            const AnchorId anchorId = journey[i];
            journeyPortionAnchorIndexes.push_back(getAnchorIndex(anchorId));
        }
    }
}



// Gather all transitions(anchorIndex0, anchorIndex1) for consecutive
// anchors in the journey portions. The number of times each
// transition appears in the journeys is its coverage.
// Store the transitions in a vector indexed by coverage.
void RestrictedAnchorGraph::gatherTransitions(ostream& html)
{
    // First, gather all transitions and store the anchorIndexes1
    // in a vector indexed by anchorIndex0.
    vector< vector<uint64_t> > anchorIndexes1(allAnchorIds.size());
    for(const vector<uint64_t>& journeyPortionAnchorIndexes: journeyPortionsAnchorIndexes) {
        for(uint64_t i1=1; i1<journeyPortionAnchorIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t anchorIndex0 = journeyPortionAnchorIndexes[i0];
            const uint64_t anchorIndex1 = journeyPortionAnchorIndexes[i1];
            SHASTA2_ASSERT(anchorIndex0 < allAnchorIds.size());
            SHASTA2_ASSERT(anchorIndex1 < allAnchorIds.size());
            anchorIndexes1[anchorIndex0].push_back(anchorIndex1);
        }
    }



    // Now for each anchorIndex0 deduplicate and count the anchorIndexes1.
    // Store its transitions according to its coverage.
    transitions.clear();
    vector<uint64_t> coverage;
    for(uint64_t anchorIndex0=0; anchorIndex0<allAnchorIds.size(); anchorIndex0++) {
        vector<uint64_t>& v = anchorIndexes1[anchorIndex0];
        deduplicateAndCount(v, coverage);

        for(uint64_t i=0; i<v.size(); i++) {
            const uint64_t c = coverage[i];
            if(c >= transitions.size()) {
                transitions.resize(c + 1);
            }
            transitions[c].push_back(Transition({anchorIndex0, v[i]}));
        }
    }

    if(html) {
        html << "<h3>Number of journey transitions by coverage</h3>"
            "<table><tr><th>Coverage<th>Frequency";
        for(uint64_t coverage=0; coverage<transitions.size(); coverage++) {
            const uint64_t frequency = transitions[coverage].size();
            if(frequency) {
                html << "<tr><td class=centered>" << coverage <<
                    "<td class=centered>" << frequency;
            }
        }
        html << "</table>";
    }

    if(transitions.empty()) {
        throw NoTransitions();
    }
}





// Fill the journey portions using a TangleMatrix1.
void RestrictedAnchorGraph::fillJourneyPortions(
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double drift = 0.2;

    using OrientedReadInfo = TangleMatrix1::OrientedReadInfo;

    const Anchors& anchors = tangleMatrix1.assemblyGraph.anchors;
    const Markers& markers = anchors.markers;

    const vector<OrientedReadInfo>& entranceOrientedReadInfos = tangleMatrix1.entranceOrientedReadInfos[iEntrance];
    const vector<OrientedReadInfo>& exitOrientedReadInfos = tangleMatrix1.exitOrientedReadInfos[iExit];



    // Loop over the common oriented reads to estimate an offset.
    auto itEntrance = entranceOrientedReadInfos.begin();
    auto itExit = exitOrientedReadInfos.begin();
    const auto itEntranceEnd = entranceOrientedReadInfos.end();
    const auto itExitEnd = exitOrientedReadInfos.end();
    uint64_t offsetSum = 0;
    uint64_t offsetCount = 0;
    while((itEntrance!=itEntranceEnd) and (itExit!=itExitEnd)) {
        if(itEntrance->orientedReadId < itExit->orientedReadId) {
            ++itEntrance;
            continue;
        }
        if(itExit->orientedReadId < itEntrance->orientedReadId) {
            ++itExit;
            continue;
        }
        const OrientedReadId orientedReadId = itEntrance->orientedReadId;
        SHASTA2_ASSERT(orientedReadId == itExit->orientedReadId);

        if(not tangleMatrix1.goesBackward(orientedReadId)) {
            const Journey journey = journeys[orientedReadId];
            const uint32_t entrancePositionInJourney = itEntrance->positionInJourney;
            const uint32_t exitPositionInJourney = itExit->positionInJourney;
            const AnchorId entranceAnchorId = journey[entrancePositionInJourney];
            const AnchorId exitAnchorId = journey[exitPositionInJourney];
            const uint32_t entranceOrdinal = anchors.getOrdinal(entranceAnchorId, orientedReadId);
            const uint32_t exitOrdinal = anchors.getOrdinal(exitAnchorId, orientedReadId);
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const uint64_t offset =
                orientedReadMarkers[exitOrdinal].position -
                orientedReadMarkers[entranceOrdinal].position;
            offsetSum += offset;
            ++offsetCount;
        }

        ++itEntrance;
        ++itExit;
    }
    const double averageOffset = double(offsetSum) / double(offsetCount);
    const uint32_t maxOffset = uint32_t(std::round(averageOffset * (1. + drift)));



    // Joint loop over the oriented reads of the entrance and exit.
    itEntrance = entranceOrientedReadInfos.begin();
    itExit = exitOrientedReadInfos.begin();
    while(true) {

        // If both iterators are at their end, we are done.
        if((itEntrance==itEntranceEnd) and (itExit==itExitEnd)) {
            break;
        }



        // Case where the OrientedReadId appears only in the entrance.
        // In this case we use a journey portion that begins at the
        // journey position stored in the OrientedReadInfo and
        // ends at a distance maxOffset after that.
        if(
            (itEntrance != itEntranceEnd) and
            (itExit==itExitEnd or (itEntrance->orientedReadId < itExit->orientedReadId))) {
            const OrientedReadId orientedReadId = itEntrance->orientedReadId;
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Journey journey = journeys[orientedReadId];
                const uint32_t begin = itEntrance->positionInJourney;
                uint32_t end = uint32_t(journey.size());

                // End the JourneyPortion at a distance maxOffset after begin.
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                uint32_t beginPosition = invalid<uint32_t>;
                for(uint32_t positionInJourney=begin; positionInJourney<end; positionInJourney++) {
                    const AnchorId anchorId = journey[positionInJourney];
                    const uint32_t ordinal = anchors.getOrdinal(anchorId, orientedReadId);
                    const uint32_t position = orientedReadMarkers[ordinal].position;
                    if(positionInJourney == begin) {
                        beginPosition = position;
                    } else {
                        const uint32_t offset = position - beginPosition;
                        if(offset >maxOffset) {
                            end = positionInJourney;
                            break;
                        }
                    }
                }

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
                uint32_t begin = 0;
                const uint32_t end = itExit->positionInJourney + 1;

                // Begin the JourneyPortion at a distance maxOffset before begin.
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                uint32_t endPosition = invalid<uint32_t>;
                for(uint32_t positionInJourney=end-1; /* Check later */ ; positionInJourney--) {
                    const AnchorId anchorId = journey[positionInJourney];
                    const uint32_t ordinal = anchors.getOrdinal(anchorId, orientedReadId);
                    const uint32_t position = orientedReadMarkers[ordinal].position;
                    if(positionInJourney == end - 1) {
                        endPosition = position;
                    } else {
                        const uint32_t offset = endPosition - position;
                        if(offset > maxOffset) {
                            begin = positionInJourney + 1;
                            break;
                        }
                    }

                    if(positionInJourney == 0) {
                        break;
                    }
                }

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
            SHASTA2_ASSERT(0);
        }
    }

    if(html) {
        html << "<h4>Journey portions used in this RestrictedAnchorGraph</h4>"
            "<table><tr>"
            "<th>OrientedReadId<th>Begin<th>End<th>Length";
        for(const JourneyPortion journeyPortion: journeyPortions) {
            html << "<tr>"
                "<td class=centered>" << journeyPortion.orientedReadId <<
                "<td class=centered>" << journeyPortion.begin <<
                "<td class=centered>" << journeyPortion.end <<
                "<td class=centered>" << journeyPortion.end - journeyPortion.begin;
        }
        html << "</table>";
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


    // Gather the AnchorIds that appear in the JourneyPortions.
    // Each distinct AnchorId will generate a vertex.
    vertexTable.clear();
    for(const JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i=begin; i<end; i++) {
            const AnchorId anchorId = journey[i];
            vertexTable.push_back(make_pair(anchorId, null_vertex()));
        }
    }
    deduplicate(vertexTable);

    // Each distinct AnchorId generates a vertex.
    for(pair<AnchorId, vertex_descriptor>& p: vertexTable) {
        const AnchorId anchorId = p.first;
        const vertex_descriptor v = add_vertex(RestrictedAnchorGraphVertex(anchorId), *this);
        p.second = v;
    }

    // Loop over the JourneyPortions to add oriented reads to the vertices.
    for(const JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i=begin; i<end; i++) {
            const AnchorId anchorId = journey[i];
            const vertex_descriptor v = getExistingVertex(anchorId);
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
            const vertex_descriptor v0 = getExistingVertex(anchorId0);
            const vertex_descriptor v1 = getExistingVertex(anchorId1);

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



// Create a new vertex and add it to the vertexTable, without resorting
// the vertexTable. This invalidates the vertexTable.
// If this is called with the AnchorId of an existing vertex,
// the call succeeds but the subsequent call to sortVertexTable will assert.
RestrictedAnchorGraph::vertex_descriptor RestrictedAnchorGraph::addVertex(AnchorId anchorId)
{

    // Add the vertex and update the vertexTable.
    // The vertexTable is no longer sorted and so it becomes invalid.
    vertex_descriptor v = null_vertex();
    if(cycleAvoider) {
        v = cycleAvoider->addVertex();
        (*this)[v].anchorId = anchorId;
    }
    else {
        v = add_vertex(RestrictedAnchorGraphVertex(anchorId), *this);
    }

    vertexTable.push_back(make_pair(anchorId, v));
    vertexTableIsValid = false;

    return v;
}



void RestrictedAnchorGraph::sortVertexTable()
{
    sort(vertexTable.begin(), vertexTable.end());

    // Check that we don't have any duplicates.
    for(uint64_t i=1; i<vertexTable.size(); i++) {
        SHASTA2_ASSERT(vertexTable[i-1].first != vertexTable[i].first);
    }

    vertexTableIsValid = true;
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

    const vertex_descriptor v0 = getExistingVertex(anchorId0);
    const vertex_descriptor v1 = getExistingVertex(anchorId1);

    std::set<vertex_descriptor> reachableVertices;
    vector<vertex_descriptor> verticesToBeRemoved;

    // Remove vertices that are not forward reachable from v0.
    // To permit future optimizations, we don't really remove
    // the vertices - we only disconnect them for the rest of the graph.
    findReachableVertices(graph, v0, 0, reachableVertices);
    SHASTA2_ASSERT(reachableVertices.contains(v1));
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
    SHASTA2_ASSERT(reachableVertices.contains(v0));
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
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
    vector<edge_descriptor> edgesSortedByDecreasingCoverage;
    for(const auto& p: edgeTable) {
        edgesSortedByDecreasingCoverage.push_back(p.first);
    }

    // Do the approximate topological sort.
    shasta::approximateTopologicalSort(graph, edgesSortedByDecreasingCoverage);
}



// Find the optimal assembly path.
void RestrictedAnchorGraph::findOptimalPath(
    AnchorId anchorId0,
    AnchorId anchorId1,
    vector<edge_descriptor>& optimalPath)
{
    using Graph = RestrictedAnchorGraph;
    Graph& graph = *this;

    // Find the vertices corresponding to anchorId0 and anchorId1.
    const vertex_descriptor v0 = getExistingVertex(anchorId0);
    const vertex_descriptor v1 = getExistingVertex(anchorId1);

    // Compute the longest path.
    shasta::longestPath(graph, optimalPath);
    if(
        optimalPath.empty() or
        (source(optimalPath.front(), graph) != v0) or
        (target(optimalPath.back(), graph) != v1)) {
        ofstream dot("RestrictedAnchorGraphFailure.dot");
        writeGraphviz(dot, vector<AnchorId>{anchorId0, anchorId1});
    }
    SHASTA2_ASSERT(not optimalPath.empty());
    SHASTA2_ASSERT(source(optimalPath.front(), graph) == v0);
    SHASTA2_ASSERT(target(optimalPath.back(), graph) == v1);


    // Set the isOptimalPathEdge flags.
    for(const edge_descriptor e: optimalPath) {
        graph[e].isOptimalPathEdge = true;
    }


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
    ostream& html)
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
