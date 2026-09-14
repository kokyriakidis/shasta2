// Shasta2.
#include "StrandSplitter.hpp"
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "graphvizToHtml.hpp"
#include "html.hpp"
#include "Options.hpp"
using namespace shasta2;

// Standard library.



// This takes as input a self-complementary tangle
// in an AssemblyGraph and attempts to split the strands.
// The last two arguments are only used for debug output.
StrandSplitter::StrandSplitter(
    AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
    uint64_t tangleId,
    const string& debugOutputBaseName) :
    assemblyGraph(assemblyGraph),
    tangleId(tangleId),
    debugOutputBaseName(debugOutputBaseName),
    tangle(assemblyGraph, tangleVertices)
{
    writeInitialDebugOutput();
    SHASTA2_ASSERT(tangle.isSelfComplementary());
    gatherSegments();
    writeLowCoverageSegmentPairs();
    findReadOccurrences();
    createStrandSeparationGraph();

    if(not separateStrands()) {
        return;
    }

    createConnectionGraph();

#if 0
    findHangingSegments();
    findCandidateConnections();
#endif
}



StrandSplitter::~StrandSplitter()
{
    if(debug) {
        writeHtmlEnd(html);
    }
}



void StrandSplitter::writeInitialDebugOutput()
{
    if(debug) {
        cout << "StrandSplitter begins for tangle " << tangleId << endl;
        html.open(debugOutputBaseName + "-StrandSplitter-Tangle-" + to_string(tangleId) + ".html");
        writeHtmlBegin(html, "Tangle " + to_string(tangleId));
        html << "<h1>Self-complementary tangle " << tangleId << "</h1>";
        tangle.writeHtml(html);
    }

}


void StrandSplitter::gatherSegments()
{
    // Fill in allTangleSegments.
    const auto inserter = back_inserter(allTangleSegments);
    std::ranges::copy(tangle.tangleEdges, inserter);
    std::ranges::copy(tangle.entrances, inserter);
    std::ranges::copy(tangle.exits, inserter);
    sort(allTangleSegments.begin(), allTangleSegments.end(), assemblyGraph.orderById);

    // Fill in lowCoverageSegmentPairs.
    for(const Segment segment: allTangleSegments) {
        if(assemblyGraph[segment].lengthWeightedAverageCoverage() > maxCoverage) {
            continue;
        }
        const Segment segmentRc = assemblyGraph[segment].eRc;
        SHASTA2_ASSERT(segmentRc != segment);
        if(id(segment) < id(segmentRc)) {
            lowCoverageSegmentPairs.push_back({segment, segmentRc});
        }
    }

    // Fill in the lowCoverageSegments vector.
    for(const auto&[segment, segmentRc]: lowCoverageSegmentPairs) {
        lowCoverageSegments.push_back(segment);
        lowCoverageSegments.push_back(segmentRc);
    }

}



uint64_t StrandSplitter::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandSplitter::writeLowCoverageSegmentPairs()
{
    if(debug) {
        html << "<h2>Reverse complemented pairs of low coverage segments</h2>"
            "Low coverage segments are tangle segments "
            " (internal, entrances, exits) with coverage "
            "no greater than " << maxCoverage << ". These are the ones "
            "that are considered reliably single-copy and are used "
            "for strand separation. "
            "They don't necessarily include entrances and exits, "
            "but they can."
            "<br><br><table><tr><th>Pair index<th>Index0<th>Index1<th>Segment0<th>Segment1"
            "<th>Length<th>Coverage<th>Entrance<br>or<br>Exit";
        for(uint64_t segmentPairIndex=0; segmentPairIndex<lowCoverageSegmentPairs.size(); segmentPairIndex++) {
            const auto&[segment, segmentRc] = lowCoverageSegmentPairs[segmentPairIndex];
            const uint64_t length = assemblyGraph[segment].length();
            SHASTA2_ASSERT(length == assemblyGraph[segmentRc].length());
            const double coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
            SHASTA2_ASSERT(coverage == assemblyGraph[segmentRc].lengthWeightedAverageCoverage());
            html << "<tr>"
                "<td class=centered>" << segmentPairIndex <<
                "<td class=centered>" << 2*segmentPairIndex <<
                "<td class=centered>" << 2*segmentPairIndex+1 <<
                "<td class=centered>" << id(segment) <<
                "<td class=centered>" << id(segmentRc) <<
                "<td class=centered>" << length <<
                "<td class=centered>" << coverage <<
                "<td class=centered>";
            if(isEntrance(segment) or isExit(segment)) {
                html << "&check;";
            }

        }
        html << "</table>";
    }

}



void StrandSplitter::findReadOccurrences()
{
    for(uint64_t segmentPairIndex=0; segmentPairIndex<lowCoverageSegmentPairs.size(); segmentPairIndex++) {
        const auto&[segment, ignore] = lowCoverageSegmentPairs[segmentPairIndex];
        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();
                readOccurrenceMap[readId].emplace_back(ReadOccurrence(segmentPairIndex, strand));
            }
        }
    }

    // Deduplicate and count the Occurrences for each ReadId.
    vector<uint64_t> count;
    for(auto&[readId, occurrences]: readOccurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }

}



void StrandSplitter::StrandSeparationGraph::addToEdge(
    uint64_t segmentIndex0,
    uint64_t segmentIndex1,
    uint64_t frequency
    )
{
    auto[e, edgeExists] = boost::edge(segmentIndex0, segmentIndex1, *this);
    if(edgeExists) {
        (*this)[e].frequency += frequency;
    } else {
        boost::add_edge(segmentIndex0, segmentIndex1, StrandSeparationEdge(frequency), *this);
    }
}



void StrandSplitter::createStrandSeparationGraph()
{
    for(const Segment segment: lowCoverageSegments) {
        boost::add_vertex(segment, strandSeparationGraph);
    }
    for(const auto&[ignore, occurrences]: readOccurrenceMap) {
        for(uint64_t i0=0; i0<occurrences.size(); i0++) {
            const ReadOccurrence& occurrence0 = occurrences[i0];
            const uint64_t segmentPairIndex0 = occurrence0.segmentPairIndex;
            const uint64_t strand0 = occurrence0.strand;
            for(uint64_t i1=i0+1; i1<occurrences.size(); i1++) {
                const ReadOccurrence& occurrence1 = occurrences[i1];
                const uint64_t segmentPairIndex1 = occurrence1.segmentPairIndex;
                const uint64_t strand1 = occurrence1.strand;
                const uint64_t frequency = occurrence0.frequency * occurrence1.frequency;
                if(strand0 == strand1) {
                    uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                    uint64_t segmentIndex1 = 2 * segmentPairIndex1;
                    strandSeparationGraph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                    ++segmentIndex0;
                    ++segmentIndex1;
                    strandSeparationGraph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                } else {
                    uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                    uint64_t segmentIndex1 = 2 * segmentPairIndex1 + 1;
                    strandSeparationGraph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                    ++segmentIndex0;
                    --segmentIndex1;
                    strandSeparationGraph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                }
            }
        }
    }

    strandSeparationGraph.findEdgePairs();

}



void StrandSplitter::StrandSeparationGraph::findEdgePairs()
{
    StrandSeparationGraph& strandSeparationGraph = *this;

    std::set<edge_descriptor> edgesFound;
    BGL_FORALL_EDGES(e, strandSeparationGraph, StrandSeparationGraph) {
        if(edgesFound.contains(e)) {
            continue;
        }
        const vertex_descriptor v0 = source(e, strandSeparationGraph);
        const vertex_descriptor v1 = target(e, strandSeparationGraph);
        const vertex_descriptor v0Rc = v0 ^ 1;
        const vertex_descriptor v1Rc = v1 ^ 1;
        auto[eRc, edgeExists] = boost::edge(v0Rc, v1Rc, strandSeparationGraph);
        SHASTA2_ASSERT(edgeExists);
        SHASTA2_ASSERT(strandSeparationGraph[eRc].frequency == strandSeparationGraph[e].frequency);

        edgesFound.insert(e);
        edgesFound.insert(eRc);

        edgePairs.emplace_back(EdgePair({e, eRc, strandSeparationGraph[e].frequency}));
    }
    sort(edgePairs.begin(), edgePairs.end());

}



bool StrandSplitter::separateStrands()
{
    // Do strand separation by adding edges in order of decreasing frequency.
    DisjointSets disjointSets(lowCoverageSegments.size());
    uint64_t crossStrandEdgeCount = 0;
    for(const auto& edgePair: strandSeparationGraph.edgePairs) {
        const StrandSeparationGraph::edge_descriptor eA = edgePair.e;
        const StrandSeparationGraph::edge_descriptor eB = edgePair.eRc;
        const uint64_t v0A = source(eA, strandSeparationGraph);
        const uint64_t v1A = target(eA, strandSeparationGraph);
        const uint64_t v0B = source(eB, strandSeparationGraph);
        const uint64_t v1B = target(eB, strandSeparationGraph);
        const uint64_t v0ARc = v0A ^ 1;
        const uint64_t v1ARc = v1A ^ 1;
        const uint64_t v0BRc = v0B ^ 1;
        const uint64_t v1BRc = v1B ^ 1;
        const bool strandViolationA = (disjointSets.findSet(v1A) == disjointSets.findSet(v0ARc));
        const bool strandViolationB = (disjointSets.findSet(v1B) == disjointSets.findSet(v0BRc));
        const bool strandViolationARc = (disjointSets.findSet(v0A) == disjointSets.findSet(v1ARc));
        const bool strandViolationBRc = (disjointSets.findSet(v0B) == disjointSets.findSet(v1BRc));
        const bool strandViolation = strandViolationA;
        SHASTA2_ASSERT(strandViolationB == strandViolation);
        SHASTA2_ASSERT(strandViolationARc == strandViolation);
        SHASTA2_ASSERT(strandViolationBRc == strandViolation);
        if(strandViolation) {
            crossStrandEdgeCount += 2;
            strandSeparationGraph[eA].isCrossStrandEdge = true;
            strandSeparationGraph[eB].isCrossStrandEdge = true;
        } else {
            disjointSets.unionSet(v0A, v1A);
            disjointSets.unionSet(v0B, v1B);
        }
    }

    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);

    // Store the component of each vertex.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        for(const uint64_t v: component) {
            strandSeparationGraph[v].component = componentId;
        }
    }


    // Write out the StrandSeparationGraph.
    if(debug) {
        const string dotFileName = debugOutputBaseName + "-StrandSplitter-StrandSeparationGraph-Tangle-" +
            to_string(tangleId) + ".dot";
        strandSeparationGraph.writeGraphviz(dotFileName, assemblyGraph);
        const double timeout = 30.;
        const string options = "-Nshape=point -Nwidth=0.2 -Gratio=expand -Gsize=10";
        html << "<h2>Strand separation graph</h2>" << dotFileName;
        try {
            graphvizToHtml(dotFileName, "sfdp", timeout, options, html, true);
        } catch (std::exception&) {
            html << "The strand separation graph took too long to display.";
        }
    }



    // If we don't have exactly two components, do nothing.
    if(components.size() != 2) {
        if(debug) {
            html << "<br>Strand separation is not successful. "
                "Expected exactly 2 components.";
        }
        return false;
    }
    SHASTA2_ASSERT(components[0].size() == components[1].size());

    // Each component corresponds to a strand.
    // Store their Segments.
    for(uint64_t strand=0; strand<2; strand++) {
        const vector<uint64_t>& component = components[strand];
        for(uint64_t segmentIndex: component) {
            strandSegments[strand].push_back(strandSeparationGraph[segmentIndex].segment);
        }
        sort(strandSegments[strand].begin(), strandSegments[strand].end(),
            assemblyGraph.orderById);
    }

    if(debug) {
        for(uint64_t strand=0; strand<2; strand++) {
            html << "<h2>Strand " << strand << " segments</h2>";
            for(uint64_t i=0; i<strandSegments[strand].size(); i++) {
                if(i != 0) {
                    html << ",<wbr>";
                }
                html << id(strandSegments[strand][i]);
            }
        }

    }

    return true;
}



void StrandSplitter::StrandSeparationGraph::writeGraphviz(
    const string& fileName,
    const AssemblyGraph& assemblyGraph) const
{
    const StrandSeparationGraph& strandSeparationGraph = *this;

    ofstream dot(fileName);

    dot << "graph StrandSeparationGaph {\n";
    BGL_FORALL_VERTICES(segmentIndex, strandSeparationGraph, StrandSeparationGraph) {
        const string color = randomHslColor(strandSeparationGraph[segmentIndex].component, 0.75, 0.5);
        dot << assemblyGraph.id(strandSeparationGraph[segmentIndex].segment) <<
            " [style=filled fillcolor=\"" << color << "\"]"
            "\n";
    }

    BGL_FORALL_EDGES(e, strandSeparationGraph, StrandSeparationGraph) {
        const uint64_t segmentIndex0 = source(e, strandSeparationGraph);
        const uint64_t segmentIndex1 = target(e, strandSeparationGraph);
        dot << assemblyGraph.id(strandSeparationGraph[segmentIndex0].segment) << "--";
        dot << assemblyGraph.id(strandSeparationGraph[segmentIndex1].segment);
        if(strandSeparationGraph[e].isCrossStrandEdge) {
            dot << "[color=red]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



void StrandSplitter::createConnectionGraph()
{
    // Create vertices of the ConnectionGraph.
    // There is a vertex for each strand 0 Segment plus
    // a vertex for each entrance or exit that is not also a
    // strand0 segment.
    for(const Segment segment: strandSegments[0]) {
        connectionGraph.addVertex(segment);
    }
    for(const Segment segment: tangle.entrances) {
        connectionGraph.addVertex(segment);
    }
    for(const Segment segment: tangle.exits) {
        connectionGraph.addVertex(segment);
    }



    // Add the edges that correspond to connections already present
    // in the AssemblyGraph.
    BGL_FORALL_VERTICES(v0, connectionGraph, ConnectionGraph) {
        const Segment segment0 = connectionGraph[v0].segment;
        const ConnectionGraph::vertex_descriptor v1 = target(segment0, assemblyGraph);
        BGL_FORALL_OUTEDGES(v1, segment1, assemblyGraph,AssemblyGraph) {
            const auto it1 = connectionGraph.vertexMap.find(segment1);
            if(it1== connectionGraph.vertexMap.end()) {
                continue;
            }
            const ConnectionGraph::vertex_descriptor v1 = it1->second;
            boost::add_edge(v0, v1, ConnectionEdge(), connectionGraph);
        }
    }


    // Now walk at the AssemblyGraph to find additional connections.
    // When walking the AssemblyGraph, we avoid  strand 1 Segments.
    const vector<Segment>& forbiddenSegments = strandSegments[1];
    vector<Segment> stopSegments;
    BGL_FORALL_VERTICES(v, connectionGraph, ConnectionGraph) {
        stopSegments.push_back(connectionGraph[v].segment);
    }
    sort(stopSegments.begin(), stopSegments.end(), assemblyGraph.orderById);

    const uint32_t representativeRegionStepCount = uint32_t(assemblyGraph.options.representativeRegionStepCount);
    ostream noOutput(0);
    vector<Segment> reachableStopSegments;
    BGL_FORALL_VERTICES(vA, connectionGraph, ConnectionGraph) {
        const Segment segmentA = connectionGraph[vA].segment;
        const AssemblyGraph::vertex_descriptor vA0 = source(segmentA, assemblyGraph);
        const AssemblyGraph::vertex_descriptor vA1 = target(segmentA, assemblyGraph);

        // Look forward.
        assemblyGraph.bfs(vA1, 0, forbiddenSegments, stopSegments, reachableStopSegments);
        for(const Segment segmentB: reachableStopSegments) {
            const auto itB = connectionGraph.vertexMap.find(segmentB);
            SHASTA2_ASSERT(itB != connectionGraph.vertexMap.end());
            const ConnectionGraph::vertex_descriptor vB = itB->second;
            const auto[ignore, edgeExists] = boost::edge(vA, vB, connectionGraph);
            if((not edgeExists) and assemblyGraph.canConnect(segmentA, segmentB, false)) {
                boost::add_edge(vA, vB, ConnectionEdge(
                    SegmentStepSupport::analyzeSegmentPair(noOutput, assemblyGraph, segmentA, segmentB,
                    representativeRegionStepCount)),
                    connectionGraph);
            }
        }

        // Look backward.
        assemblyGraph.bfs(vA0, 1, forbiddenSegments, stopSegments, reachableStopSegments);
        for(const Segment segmentB: reachableStopSegments) {
            const auto itB = connectionGraph.vertexMap.find(segmentB);
            SHASTA2_ASSERT(itB != connectionGraph.vertexMap.end());
            const ConnectionGraph::vertex_descriptor vB = itB->second;
            const auto[ignore, edgeExists] = boost::edge(vB, vA, connectionGraph);
            if((not edgeExists) and assemblyGraph.canConnect(segmentB, segmentA, false)) {
                boost::add_edge(vB, vA, ConnectionEdge(
                    SegmentStepSupport::analyzeSegmentPair(noOutput, assemblyGraph, segmentB, segmentA,
                    representativeRegionStepCount)),
                    connectionGraph);
            }
        }
    }



    if(debug) {
        const string dotFileName = debugOutputBaseName + "-StrandSplitterConnectionGraph-Tangle-" +
            to_string(tangleId) + ".dot";
        connectionGraph.writeGraphviz(dotFileName, assemblyGraph);
        const double timeout = 30.;
        const string options = "-Nshape=rectangle";
        html << "<h2>Connection separation graph</h2>" << dotFileName;
        try {
            graphvizToHtml(dotFileName, "dot", timeout, options, html, true);
        } catch (std::exception&) {
            html << "The connection graph took too long to display.";
        }
    }
}



void StrandSplitter::ConnectionGraph::writeGraphviz(
    const string& fileName,
    const AssemblyGraph& assemblyGraph) const
{
    const ConnectionGraph& connectionGraph = *this;

    ofstream dot(fileName);
    dot << "digraph ConnectionGraph {\n";

    BGL_FORALL_VERTICES(v, connectionGraph, ConnectionGraph) {
        const Segment segment = connectionGraph[v].segment;
        dot << assemblyGraph.id(segment) << ";\n";
    }

    BGL_FORALL_EDGES(e, connectionGraph, ConnectionGraph) {
        const ConnectionEdge& edge = connectionGraph[e];
        const vertex_descriptor v0 = source(e, connectionGraph);
        const vertex_descriptor v1 = target(e, connectionGraph);
        const Segment segment0 = connectionGraph[v0].segment;
        const Segment segment1 = connectionGraph[v1].segment;
        dot <<
            assemblyGraph.id(segment0) << "->" <<
            assemblyGraph.id(segment1) << " [";

        if(edge.isDirectConnection) {
            dot << "color=green";
        } else {
            dot <<
                "label=\"" << edge.segmentPairInformation.commonCount <<
            "/" << edge.segmentPairInformation.missing() << "\"";
        }
        dot << "];\n";
    }

    dot << "}\n";
}



void StrandSplitter::ConnectionGraph::addVertex(Segment segment)
{
    ConnectionGraph& connectionGraph = *this;

    if(not vertexMap.contains(segment)) {
        const vertex_descriptor v = boost::add_vertex(ConnectionVertex(segment), connectionGraph);
        vertexMap.insert({segment, v});
    }
}


bool StrandSplitter::isEntrance(Segment segment) const
{
    return std::binary_search(tangle.entrances.begin(), tangle.entrances.end(),
        segment, assemblyGraph.orderById);
}



bool StrandSplitter::isExit(Segment segment) const
{
    return std::binary_search(tangle.exits.begin(), tangle.exits.end(),
        segment, assemblyGraph.orderById);
}




bool StrandSplitter::isStrand0Segment(Segment segment) const
{
    return std::binary_search(strandSegments[0].begin(), strandSegments[0].end(),
        segment, assemblyGraph.orderById);

}



bool StrandSplitter::isStrand1Segment(Segment segment) const
{
    return std::binary_search(strandSegments[1].begin(), strandSegments[1].end(),
        segment, assemblyGraph.orderById);
}

