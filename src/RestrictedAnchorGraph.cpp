// Shasta.
#include "RestrictedAnchorGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "longestPath.hpp"
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
    ostream& html) :
    filteringPredicate(this),
    filteredGraph(*this, filteringPredicate, filteringPredicate)
{
    constructFromTangleMatrix1(anchors, journeys, tangleMatrix1, iEntrance, iExit, html);
}



// Original version.
void RestrictedAnchorGraph::constructFromTangleMatrix(
    const Anchors& anchors,
    const Journeys& journeys,
    const TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    fillJourneyPortions(journeys, tangleMatrix1, iEntrance, iExit);
    create(anchors, journeys, html);

    const AssemblyGraph::edge_descriptor entrance = tangleMatrix1.entrances[iEntrance];
    const AssemblyGraph::edge_descriptor exit = tangleMatrix1.exits[iExit];
    const AnchorId anchorId0 = tangleMatrix1.assemblyGraph[entrance].back().anchorPair.anchorIdB;
    const AnchorId anchorId1 = tangleMatrix1.assemblyGraph[exit].front().anchorPair.anchorIdA;

    removeLowCoverageEdges(anchorId0, anchorId1);
    keepBetween(anchorId0, anchorId1);
    removeCycles();
    keepBetween(anchorId0, anchorId1);
}



// More efficient version
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

    fillJourneyPortions(journeys, tangleMatrix1, iEntrance, iExit);
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
            SHASTA2_ASSERT(0);
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
    SHASTA2_ASSERT(not optimalPath.empty());
    SHASTA2_ASSERT(source(optimalPath.front(), graph) == v0);
    SHASTA2_ASSERT(target(optimalPath.back(), graph) == v1);

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
    const vertex_descriptor v0 = getExistingVertex(anchorId0);
    const vertex_descriptor v1 = getExistingVertex(anchorId1);

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
                    SHASTA2_ASSERT(v != v0);
                    SHASTA2_ASSERT(v != v1);
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
    SHASTA2_ASSERT(dominatorSequence.size() >= 2);
    SHASTA2_ASSERT(dominatorSequence.front() == v0);
    SHASTA2_ASSERT(dominatorSequence.back() == v1);



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
