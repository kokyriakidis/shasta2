// Shasta.
#include "LocalAssembly7.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Reads.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <cmath>
#include "fstream.hpp"
#include <iomanip>
#include "iostream.hpp"
#include <map>



LocalAssembly7::LocalAssembly7(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    ostream& html,
    const vector<OrientedReadId>& orientedReadIds) :
    anchors(anchors),
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB),
    html(html)
{
    if(html) {
        html << "<h3>Local assembly between " <<
            anchorIdToString(anchorIdA) << " and " <<
            anchorIdToString(anchorIdB) << "</h3>";
    }

    gatherOrientedReads(orientedReadIds);
    estimateOffset();
    gatherSequences();
    writeOrientedReads();
    writeSequences();

    createGraph();
    graph.disconnectUnreachableVertices();

    // If the graph has cycles, declare failure and stop here.
    {
        std::map<vertex_descriptor, uint64_t> componentMap;
        if(boost::strong_components(graph, boost::make_assoc_property_map(componentMap)) != num_vertices(graph)) {
            if(html) {
                html << "<br>The De Bruijn graph contains cycles.";
            }
            return;
        }
    }


    graph.computeAssemblyPath();
    writeGraph();

}



void LocalAssembly7::gatherOrientedReads(
    const vector<OrientedReadId>& orientedReadIds)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    if(html) {
        html << "<h4>" << orientedReadIds.size() << " input oriented reads</h4><table>";
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }



    // Joint loop over the input OrientedReadIds and
    // the OrientedReadIds of the two Anchors.
    auto itA = anchorA.begin();
    const auto endA = anchorA.end();
    auto itB = anchorB.begin();
    const auto endB = anchorB.end();
    for(const OrientedReadId orientedReadId: orientedReadIds) {

        // Check if this OrientedReadId is on the two anchors.
        while((itA != endA) and (itA->orientedReadId < orientedReadId)) {
            ++itA;
        }
        while((itB != endB) and (itB->orientedReadId < orientedReadId)) {
            ++itB;
        }

        const bool isOnA = (itA != endA) and (itA->orientedReadId == orientedReadId);
        const bool isOnB = (itB != endB) and (itB->orientedReadId == orientedReadId);

        // If on neither anchor, this OrientedReadId cannot be used.
        if(not (isOnA or isOnB)) {
            continue;
        }

        // If on both anchors and negative offset, this OrientedReadId cannot be used.
        if(isOnA and isOnB and (itA->position > itB->position)) {
            continue;
        }

        // We can use this OrientedReadId for assembly.

        // Create an OrientedReadInfo.
        OrientedReadInfo& info = orientedReadInfos.emplace_back();
        info.orientedReadId = orientedReadId;
        if(isOnA) {
            info.positionA = itA->position;
        }
        if(isOnB) {
            info.positionB = itB->position;
        }

    }
}



void LocalAssembly7::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        if(orientedReadInfo.isOnBothAnchors()) {
            sum += orientedReadInfo.positionOffsetAB();
            ++n;
        }
    }

    SHASTA2_ASSERT(n > 0);
    offset = uint32_t(std::round(double(sum) / double(n)));
}



void LocalAssembly7::gatherSequences()
{
    // For reads fixed on one side only, we use a sequence length
    // equal to offset + aDrift * offset + bDrift.
    const uint32_t length = offset + uint32_t(std::round(aDrift * double(offset) + bDrift));


    // Fill in the beginPosition and endPosition of each OrientedReadInfo.
    // Those are the position ranges of the sequences that will be used for assembly.
    for(OrientedReadInfo& info: orientedReadInfos) {

        if(info.isOnBothAnchors()) {
            info.positionBegin = info.positionA;
            info.positionEnd  = info.positionB;
        } else if(info.isOnAnchorA()) {
            SHASTA2_ASSERT(not info.isOnAnchorB());
            const ReadId readId = info.orientedReadId.getReadId();
            const uint32_t readLength = uint32_t(anchors.reads.getReadSequenceLength(readId));
            info.positionBegin = info.positionA;
            info.positionEnd  = min(readLength, info.positionBegin + length);
        } else if(info.isOnAnchorB()) {
            SHASTA2_ASSERT(not info.isOnAnchorA());
            info.positionEnd = info.positionB;
            if(info.positionEnd > length) {
                info.positionBegin = info.positionEnd - length;
            } else {
                info.positionBegin  = 0;
            }
        } else {
            SHASTA2_ASSERT(0);
        }
    }



    // Now we can create the sequences.
    const Reads& reads = anchors.reads;
    vector<Base> sequence;
    for(OrientedReadInfo& info: orientedReadInfos) {
        const OrientedReadId orientedReadId = info.orientedReadId;

        // Create the sequence of this read (portion to be used for this local assembly).
        sequence.clear();
        for(uint32_t position=info.positionBegin; position!=info.positionEnd; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }

        // See if we already have this sequence in this table.
        bool found = false;
        for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
            SequenceInfo& sequenceInfo = sequences[sequenceId];
            if(
                (info.isOnAnchorA() == sequenceInfo.isOnAnchorA) and
                (info.isOnAnchorB() == sequenceInfo.isOnAnchorB) and
                (sequence == sequenceInfo.sequence)) {
                info.sequenceId = sequenceId;
                sequenceInfo.orientedReadIds.push_back(orientedReadId);
                found = true;
                break;
            }
        }
        if(not found) {
            info.sequenceId = sequences.size();
            sequences.emplace_back(k, info.isOnAnchorA(), info.isOnAnchorB(), orientedReadId, sequence);
        }
    }
}



void LocalAssembly7::writeOrientedReads() const
{
    if(not html) {
        return;
    }

    html << "<h4>" << orientedReadInfos.size() << " usable oriented reads</h4>";

    html << "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>On A"
        "<th>On B"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length<br>(bases)"
        "<th>PositionAB<br>offset"
        "<th>Assembly<br>position<br>begin"
        "<th>Assembly<br>position<br>end"
        "<th>Assembly<br>position<br>length"
        "<th>Sequence<br>id"
        ;


    uint64_t commonCount = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        if(info.isOnBothAnchors()) {
            ++commonCount;
        }
        html <<
            "<tr>"
            "<th class=centered>" << info.orientedReadId <<
            "<td class=centered>" << (info.isOnAnchorA() ? "&check;" : "") <<
            "<td class=centered>" << (info.isOnAnchorB() ? "&check;" : "") <<
            "<td class=centered>" << (info.isOnAnchorA() ? to_string(info.positionA) : "") <<
            "<td class=centered>" << (info.isOnAnchorB() ? to_string(info.positionB) : "") <<
            "<td class=centered>" << anchors.reads.getReadSequenceLength(info.orientedReadId.getReadId()) <<
            "<td class=centered>" << (info.isOnBothAnchors() ? to_string(info.positionOffsetAB()) : "") <<
            "<td class=centered>" << info.positionBegin <<
            "<td class=centered>" << info.positionEnd <<
            "<td class=centered>" << info.positionOffsetForAssembly() <<
            "<td class=centered>" << info.sequenceId
            ;
    }

    html << "</table>";

    html << "<br>Estimated offset using " << commonCount <<
        " oriented reads common to the left and right anchors is " << offset << " bases.";


}



void LocalAssembly7::writeSequences() const
{
    if(not html) {
        return;
    }

    html << "<h4>Oriented read sequences used for assembly</h4>"
        "<table><tr>"
        "<th>Sequence<br>id"
        "<th>On A"
        "<th>On B"
        "<th>Length"
        "<th>Coverage"
        "<th class=left>Sequence<br>RLE sequence<br>Repeat count";

    std::map<char, uint32_t> repeatCountLegend;
    bool highRepeatCountIsPresent = false;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        html <<
            "<tr>"
            "<td class=centered>" << sequenceId <<
            "<td class=centered>" << (sequenceInfo.isOnAnchorA ? "&check;" : "") <<
            "<td class=centered>" << (sequenceInfo.isOnAnchorB ? "&check;" : "") <<
            "<td class=centered>" << sequenceInfo.sequence.size() <<
            "<td class=centered>" << sequenceInfo.orientedReadIds.size() <<
            "<td class=left style='font-family:monospace;white-space: nowrap'>";
        std::ranges::copy(sequenceInfo.sequence, ostream_iterator<Base>(html));
        html << "<br>";
        std::ranges::copy(sequenceInfo.rleSequence, ostream_iterator<Base>(html));
        html << "<br>";
        for(uint32_t repeatCount: sequenceInfo.repeatCount) {
            if(repeatCount == 1) {
                html << "&nbsp;";
            } else if(repeatCount < 10) {
                html << char('0' + repeatCount);
            } else if(repeatCount < 36) {
                const char c = char('A' + (repeatCount - 10));
                html << c;
                repeatCountLegend.insert({c, repeatCount});
            } else {
                html << '*';
                highRepeatCountIsPresent = true;
            }
        }
    }

    html << "</table>";

    // Write a repeat count legend.
    html <<
        "<br>Repeat count legend"
        "<table>"
        "<tr>"
        "<tr><th>Character<th>Repeat count"
        "<tr><td class=centered>Blank<td class=centered>1"
        "<tr><td class=centered>2-9<td class=centered>As displayed";
    for(const auto& [c, repeatCount]: repeatCountLegend) {
        html << "<tr><td class=centered>" << c << "<td class=centered>" << repeatCount;
    }
    if(highRepeatCountIsPresent) {
        html << "<tr><td class=centered>*<td class=centered>&gt;35";
    }

    html << "</table>";
}



void LocalAssembly7::SequenceInfo::constructRleSequence()
{
    for(const Base b: sequence) {
        if((not rleSequence.empty()) and (b == rleSequence.back())) {
            ++repeatCount.back();
        } else {
            rleSequence.push_back(b);
            repeatCount.push_back(1);
        }
    }
}



void LocalAssembly7::SequenceInfo::constructDeBruijnRleSequence(uint64_t k)
{
    if(isOnAnchorA) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnRleSequence.push_back(Base::fromInteger(uint8_t(10)));
            deBruijnRepeatCount.push_back(1);
        }
    }
    std::ranges::copy(rleSequence, back_inserter(deBruijnRleSequence));
    std::ranges::copy(repeatCount, back_inserter(deBruijnRepeatCount));
    if(isOnAnchorB) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnRleSequence.push_back(Base::fromInteger(uint8_t(20)));
            deBruijnRepeatCount.push_back(1);
        }
    }
}



// This can be sped up if necessary.
void LocalAssembly7::createGraph()
{
    // The k-mers of each of the sequences.
    using Kmer = vector<Base>;
    vector< vector<Kmer> > kmers;
    Kmer kmer;
    for(const SequenceInfo& sequenceInfo: sequences) {
        const vector<Base>& sequence = sequenceInfo.deBruijnRleSequence;
        vector<Kmer>& sequenceKmers = kmers.emplace_back();

        for(uint64_t position=0; position+k<=sequence.size(); position++) {
            kmer.clear();
            for(uint64_t i=0; i<k; i++) {
                kmer.push_back(sequence[position + i]);
            }
            sequenceKmers.push_back(kmer);
        }
    }


    // Gather k-mer occurrences by k-mer.
    std::map<Kmer, vector<KmerOccurrence> > kmerOccurrences;
    for(uint64_t sequenceId=0; sequenceId<kmers.size(); sequenceId++) {
        vector<Kmer>& sequenceKmers = kmers[sequenceId];
        for(uint64_t position=0; position<sequenceKmers.size(); position++) {
            const Kmer& kmer = sequenceKmers[position];
            kmerOccurrences[kmer].push_back(KmerOccurrence(sequenceId, position));
        }
    }


    // Each distinct k-mer generates a vertex.
    std::map<Kmer, vertex_descriptor> vertexMap;
    for(const auto& [kmer, occurrences]: kmerOccurrences) {
        uint64_t coverage = 0;
        for(const KmerOccurrence& occurrence: occurrences) {
            coverage += sequences[occurrence.sequenceId].coverage();
        }
        const Graph::vertex_descriptor v = boost::add_vertex(Vertex(kmer, occurrences, coverage), graph);
        vertexMap.insert({kmer, v});
    }

    // Flag the A and B vertices.
    const vector<Base> kmerA(k, Base::fromInteger(uint8_t(10)));
    const vector<Base> kmerB(k, Base::fromInteger(uint8_t(20)));
    graph.vA = vertexMap.at(kmerA);
    graph.vB = vertexMap.at(kmerB);
    graph[graph.vA].isAVertex = true;
    graph[graph.vB].isBVertex = true;



    // Now create the edges.
    for(uint64_t sequenceId=0; sequenceId<kmers.size(); sequenceId++) {
        const uint64_t coverage = sequences[sequenceId].coverage();
        vector<Kmer>& sequenceKmers = kmers[sequenceId];
        for(uint64_t position1=1; position1<sequenceKmers.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            const Kmer kmer0 = sequenceKmers[position0];
            const Kmer kmer1 = sequenceKmers[position1];
            const vertex_descriptor v0 = vertexMap.at(kmer0);
            const vertex_descriptor v1 = vertexMap.at(kmer1);
            auto [e, edgeExists] = boost::edge(v0, v1, graph);
            if(edgeExists) {
                graph[e].coverage += coverage;
            } else {
                auto[e, ignore] = boost::add_edge(v0, v1, graph);
                graph[e].coverage = coverage;
            }
        }
    }

    // Compute edge weights.
    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        const double logP = logPCoefficient * double(edge.coverage);
        edge.weight = std::pow(10., -0.1 * logP);
    }

}



void LocalAssembly7::writeGraph()
{
    if(not html) {
        return;
    }

    html <<
        "<h4>De Bruijn graph</h4>"
        "The De Bruijn graph has " << graph.countNonIsolatedVertices() <<
        " non-isolated vertices and " << num_edges(graph) << " edges.";

    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    graph.writeGraphviz(dotFileName);

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle";
    html << "<br>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);
}



void LocalAssembly7::Graph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void LocalAssembly7::Graph::writeGraphviz(ostream& dot) const
{
    const Graph& graph = *this;

    dot << "digraph DeBruijn {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) == 0) and (out_degree(v, graph) == 0)) {
            continue;
        }
        dot << v << " [";
        dot << "label=\"" << graph[v].coverage << "\"";
        if(graph[v].isAVertex) {
            dot << " style=filled fillcolor=Pink";
        } else if(graph[v].isBVertex) {
            dot << " style=filled fillcolor=LightBlue";
        } else if(graph[v].isOnAssemblyPath) {
            dot << " style=filled fillcolor=LightGreen";
        }
        dot << "];\n";
    }

    dot << std::fixed << std::setprecision(1);
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << v0 << "->" << v1 <<
            "[penwidth=" << 0.2*double(graph[e].coverage) <<
            " tooltip=\"" << graph[e].coverage << "\"";

        if(graph[v0].isOnAssemblyPath and graph[v1].isOnAssemblyPath) {
            dot << " color=green";
        }

        dot <<
            "]"
            ";\n";
    }

    dot << "}\n";
}




void LocalAssembly7::Graph::disconnectUnreachableVertices()
{
    Graph& graph = *this;

    // Disconnect the vertices that are not reachable from vA moving forward.
    std::set<vertex_descriptor> reachableVertices;
    findReachableVertices(graph, vA, 0, reachableVertices);
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            boost::clear_vertex(v, graph);
        }
    }

    // Disconnect the vertices that are not reachable from vB moving backward.
    findReachableVertices(graph, vB, 1, reachableVertices);
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            boost::clear_vertex(v, graph);
        }
    }
}



uint64_t LocalAssembly7::Graph::countNonIsolatedVertices() const
{
    const Graph& graph = *this;

    uint64_t n = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) > 0) and (out_degree(v, graph) > 0)) {
            ++n;
        }
    }

    return n;
}



void LocalAssembly7::Graph::computeAssemblyPath()
{
    Graph& graph = *this;

    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
    dijkstra_shortest_paths(graph, vA,
       weight_map(boost::get(&Edge::weight, graph)).
       predecessor_map(boost::make_assoc_property_map(predecessorMap))
       );

    // Walk back the predecessorMap staring at vB.
    vertex_descriptor v = vB;
    while(true) {
        assemblyPath.push_back(v);
        if(v == vA) {
            break;
        }
        v = predecessorMap.at(v);
    }
    std::ranges::reverse(assemblyPath);

    for(const vertex_descriptor v: assemblyPath) {
        graph[v].isOnAssemblyPath = true;
    }
}
