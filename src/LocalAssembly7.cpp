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
#include <stack>



LocalAssembly7::LocalAssembly7(
    Method method,
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

    // Fill in the orientedReads and the sequenceInfos.
    gatherOrientedReads(orientedReadIds);
    removeOutliers();
    estimateOffset();
    gatherSequences();
    writeOrientedReads();
    writeSequences();


    switch(method) {
    case Method::Adaptive:
        runAdaptive();
        break;
    case Method::Abpoa:
        runAbpoa();
        break;
    case Method::Poasta:
        runPoasta();
        break;
    case Method::Theseus:
        runTheseus();
        break;
    case Method::DeBruijn:
        runDeBruijn();
        break;
    }
}




// Loop over DeBruijn graphs with increasing k.
void LocalAssembly7::runDeBruijn()
{
    uint64_t k = kStart;
    while(true) {
        if(html) {
            html << "<br>Using a De Bruijn graph with k = " << k << ".";
        }

        runDeBruijn(k);

        if(success) {
            break;
        } else {
            k *= 2;
            if(k > kMax) {
                if(html) {
                    html << "<br>Cannot increase k above " << kMax << ".";
                }
                break;
            }
        }
    }
}



// Construct the Dr Bruijn graph for a given k
// and, if successful, use it to assemble sequence.
void LocalAssembly7::runDeBruijn(uint64_t k)
{
    for(SequenceInfo& sequenceInfo: sequences) {
        sequenceInfo.constructDeBruijnSequence(k);
        sequenceInfo.constructKmers(k);
    }
    gatherKmers();

    Graph graph;
    createGraph(k, graph);
    graph.disconnectUnreachableVertices();


    // If the graph has cycles, give up, leaving success set to false;
    {
        std::map<vertex_descriptor, uint64_t> componentMap;
        if(boost::strong_components(graph, boost::make_assoc_property_map(componentMap)) != num_vertices(graph)) {
            if(html) {
                html << "<br>The De Bruijn graph contains cycles.";
                graph.writeVertices(kmers, "DeBruijnGraph-" + to_string(k) + ".csv");
                writeKmerOccurrences(graph, "DeBruijnGraph-KmerOccurrences-" + to_string(k) + ".csv");
                writeGraph(k, graph);
            }
            return;
        }
    }

    // Compute the optimal assembly path.
    graph.computeAssemblyPath();

    if(html) {
        graph.writeVertices(kmers, "DeBruijnGraph-" + to_string(k) + ".csv");
    }
    writeGraph(k, graph);

    // Assemble sequence.
    assemble(k, graph);
    writeSequence();
    success = true;
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

        // Create an OrientedRead.
        OrientedRead& orientedRead = orientedReads.emplace_back();
        orientedRead.orientedReadId = orientedReadId;
        if(isOnA) {
            orientedRead.positionA = itA->position;
        }
        if(isOnB) {
            orientedRead.positionB = itB->position;
        }

    }
}



void LocalAssembly7::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedRead& orientedRead: orientedReads) {
        if(orientedRead.isOnBothAnchors()) {
            sum += orientedRead.positionOffsetAB();
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
    const uint32_t length = offset + uint32_t(std::round(aExtend * double(offset) + bExtend));


    // Fill in the beginPosition and endPosition of each OrientedRead.
    // Those are the position ranges of the sequences that will be used for assembly.
    for(OrientedRead& orientedRead: orientedReads) {

        if(orientedRead.isOnBothAnchors()) {
            orientedRead.positionBegin = orientedRead.positionA;
            orientedRead.positionEnd  = orientedRead.positionB;
        } else if(orientedRead.isOnAnchorA()) {
            SHASTA2_ASSERT(not orientedRead.isOnAnchorB());
            const ReadId readId = orientedRead.orientedReadId.getReadId();
            const uint32_t readLength = uint32_t(anchors.reads.getReadSequenceLength(readId));
            orientedRead.positionBegin = orientedRead.positionA;
            orientedRead.positionEnd  = min(readLength, orientedRead.positionBegin + length);
        } else if(orientedRead.isOnAnchorB()) {
            SHASTA2_ASSERT(not orientedRead.isOnAnchorA());
            orientedRead.positionEnd = orientedRead.positionB;
            if(orientedRead.positionEnd > length) {
                orientedRead.positionBegin = orientedRead.positionEnd - length;
            } else {
                orientedRead.positionBegin  = 0;
            }
        } else {
            SHASTA2_ASSERT(0);
        }
    }



    // Now we can create the sequences.
    const Reads& reads = anchors.reads;
    vector<Base> sequence;
    for(OrientedRead& orientedRead: orientedReads) {
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;

        // Create the sequence of this read (portion to be used for this local assembly).
        sequence.clear();
        for(uint32_t position=orientedRead.positionBegin; position!=orientedRead.positionEnd; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }

        // See if we already have this sequence in this table.
        bool found = false;
        for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
            SequenceInfo& sequenceInfo = sequences[sequenceId];
            if(
                (orientedRead.isOnAnchorA() == sequenceInfo.isOnAnchorA) and
                (orientedRead.isOnAnchorB() == sequenceInfo.isOnAnchorB) and
                (sequence == sequenceInfo.sequence)) {
                orientedRead.sequenceId = sequenceId;
                sequenceInfo.orientedReadIds.push_back(orientedReadId);
                found = true;
                break;
            }
        }
        if(not found) {
            orientedRead.sequenceId = sequences.size();
            sequences.emplace_back(orientedRead.isOnAnchorA(), orientedRead.isOnAnchorB(), orientedReadId, sequence);
        }
    }
}



void LocalAssembly7::writeOrientedReads() const
{
    if(not html) {
        return;
    }

    html << "<h4>" << orientedReads.size() << " usable oriented reads</h4>";

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
    for(const OrientedRead& orientedRead: orientedReads) {
        if(orientedRead.isOnBothAnchors()) {
            ++commonCount;
        }
        html <<
            "<tr>"
            "<th class=centered>" << orientedRead.orientedReadId <<
            "<td class=centered>" << (orientedRead.isOnAnchorA() ? "&check;" : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorB() ? "&check;" : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorA() ? to_string(orientedRead.positionA) : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorB() ? to_string(orientedRead.positionB) : "") <<
            "<td class=centered>" << anchors.reads.getReadSequenceLength(orientedRead.orientedReadId.getReadId()) <<
            "<td class=centered>" << (orientedRead.isOnBothAnchors() ? to_string(orientedRead.positionOffsetAB()) : "") <<
            "<td class=centered>" << orientedRead.positionBegin <<
            "<td class=centered>" << orientedRead.positionEnd <<
            "<td class=centered>" << orientedRead.positionOffsetForAssembly() <<
            "<td class=centered>" << orientedRead.sequenceId
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
        "<th class=left>Sequence";

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
    }

    html << "</table>";
}



void LocalAssembly7::SequenceInfo::constructDeBruijnSequence(uint64_t k)
{
    deBruijnSequence.clear();

    if(isOnAnchorA) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnSequence.push_back(Base::fromInteger(uint8_t(10)));
        }
    }
    std::ranges::copy(sequence, back_inserter(deBruijnSequence));
    if(isOnAnchorB) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnSequence.push_back(Base::fromInteger(uint8_t(20)));
        }
    }
}



void LocalAssembly7::SequenceInfo::constructKmers(uint64_t k)
{
    kmers.clear();

    // This can be made faster by shifting the k-mers instead of constructing
    // all of them from scratch.
    Kmer kmer;
    kmer.reserve(k);
    for(uint64_t position=0; position+k<=deBruijnSequence.size(); position++) {
        kmer.clear();
        for(uint64_t i=0; i<k; i++) {
            kmer.push_back(deBruijnSequence[position + i]);
        }
        kmers.push_back(kmer);
    }
}



void LocalAssembly7::gatherKmers()
{
    std::map<Kmer, vector<KmerOccurrence> > kmersMap;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        vector<Kmer>& sequenceKmers = sequences[sequenceId].kmers;
        for(uint64_t position=0; position<sequenceKmers.size(); position++) {
            const Kmer& kmer = sequenceKmers[position];
            kmersMap[kmer].push_back(KmerOccurrence(sequenceId, position));
        }
    }

    kmers.clear();
    std::ranges::copy(kmersMap, back_inserter(kmers));
}



void LocalAssembly7::createGraph(uint64_t k, Graph& graph)
{

    // Vectors to store the vertex corresponding to each position of each sequence.
    vector< vector<vertex_descriptor> > vertexTable;
    for(const SequenceInfo& sequenceInfo: sequences) {
        vector<vertex_descriptor>& v = vertexTable.emplace_back();
        v.resize(sequenceInfo.deBruijnSequence.size() - k + 1, Graph::null_vertex());
    }



    // In a standard De Bruijn graph, each distinct k-mer generates a vertex.
    // Here we do the same, with one important exception:
    // if a k-mer appears more than once in one of the sequences,
    // each occurrence of that k-mer generate a separate vertex.
    // This helps avoid cycles that are short compared to read length.
    // However some merging of the vertices generated in this way will be necessary.
    vector<bool> appearsInSequence(sequences.size());
    for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
        const auto& [kmer, occurrences] = kmers[kmerId];

        // Figure out if it appears more than once in any of the sequences.
        bool appearsMoreThanOnce = false;
        std::ranges::fill(appearsInSequence, false);
        for(const KmerOccurrence& occurrence: occurrences) {
            const uint64_t sequenceId = occurrence.sequenceId;
            if(appearsInSequence[sequenceId]) {
                appearsMoreThanOnce = true;
                break;
            }
            appearsInSequence[sequenceId] = true;
        }

        // If it appears more than once, generate a vertex for each occurrence.
        if(appearsMoreThanOnce) {
            vector<KmerOccurrence> occurrenceVector(1);
            for(const KmerOccurrence& occurrence: occurrences) {
                occurrenceVector.front() = occurrence;
                const uint64_t coverage = sequences[occurrence.sequenceId].coverage();
                const Graph::vertex_descriptor v = boost::add_vertex(
                    Vertex(graph.nextVertexId++, kmerId, occurrenceVector, coverage), graph);
                vertexTable[occurrence.sequenceId][occurrence.position] = v;
            }

        } else {

            // Generate a single vertex.
            uint64_t coverage = 0;
            for(const KmerOccurrence& occurrence: occurrences) {
                coverage += sequences[occurrence.sequenceId].coverage();
            }
            const Graph::vertex_descriptor v = boost::add_vertex(
                Vertex(graph.nextVertexId++, kmerId, occurrences, coverage), graph);
            for(const KmerOccurrence& occurrence: occurrences) {
                vertexTable[occurrence.sequenceId][occurrence.position] = v;
            }
        }
    }



    // Flag the A and B vertices.
    // Use the first and last vertex of a sequence fixed on both sides.
    graph.vA = Graph::null_vertex();
    graph.vB = Graph::null_vertex();
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequence = sequences[sequenceId];
        if(sequence.isOnAnchorA and sequence.isOnAnchorB) {
            graph.vA = vertexTable[sequenceId].front();
            graph.vB = vertexTable[sequenceId].back();
        }
    }
    SHASTA2_ASSERT(graph.vA != Graph::null_vertex());
    SHASTA2_ASSERT(graph.vB != Graph::null_vertex());
    graph[graph.vA].isAVertex = true;
    graph[graph.vB].isBVertex = true;



    // Now create the edges.
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const uint64_t coverage = sequences[sequenceId].coverage();
        const vector<vertex_descriptor>& sequenceVertices = vertexTable[sequenceId];
        for(uint64_t position1=1; position1<sequenceVertices.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            const vertex_descriptor v0 = sequenceVertices[position0];
            const vertex_descriptor v1 = sequenceVertices[position1];
            auto [e, edgeExists] = boost::edge(v0, v1, graph);
            if(edgeExists) {
                graph[e].coverage += coverage;
            } else {
                auto[e, ignore] = boost::add_edge(v0, v1, graph);
                graph[e].coverage = coverage;
            }
        }
    }

    graph.writeGraphviz("DeBruijnGraph-BeforeMerge-" + to_string(k) + ".dot");
    graph.merge();

    // Compute edge weights.
    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        const double logP = logPCoefficient * double(edge.coverage);
        edge.weight = std::pow(10., -0.1 * logP);
    }

}



void LocalAssembly7::writeGraph(uint64_t k, const Graph& graph)
{
    if(not html) {
        return;
    }

    html <<
        "<h4>De Bruijn graph</h4>"
        "The De Bruijn graph has " << graph.countNonIsolatedVertices() <<
        " non-isolated vertices and " << num_edges(graph) << " edges.";

    // Write it in graphviz format.
    const string dotFileName = "DeBruijnGraph-" + to_string(k) + ".dot";
    graph.writeGraphviz(dotFileName);

    // Display it in html in svg format.
    const double timeout = 3.;
    const string options = "-Nshape=rectangle";
    html << "<br>";
    try {
        graphvizToHtml(dotFileName, "dot", timeout, options, html, true);
    } catch(std::exception&) {
        html << "Unable to display the graph." << endl;
    }
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
        const Vertex& vertex = graph[v];
        dot << vertex.vertexId << " [";
        dot << "label=\"v" << vertex.vertexId << "\\nk" << vertex.kmerId << "\\n" << vertex.coverage << "\"";
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
        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];
        dot << vertex0.vertexId << "->" << vertex1.vertexId <<
            "[penwidth=" << 0.2*double(graph[e].coverage) <<
            " tooltip=\"" << graph[e].coverage << "\"";

        if(vertex0.isOnAssemblyPath and vertex1.isOnAssemblyPath) {
            dot << " color=green";
        }

        dot <<
            "]"
            ";\n";
    }

    dot << "}\n";
}



void LocalAssembly7::Graph::writeVertices(
    const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
    const string& fileName) const
{
    ofstream csv(fileName);
    writeVertices(kmers, csv);
}



void LocalAssembly7::Graph::writeVertices(
    const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
    ostream& csv) const
{
    const Graph& graph = *this;
    csv << "VertexId,KmerId,coverage,Kmer,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) == 0) and (out_degree(v, graph) == 0)) {
            continue;
        }
        const Vertex& vertex = graph[v];
        const Kmer& kmer = kmers[vertex.kmerId].first;
        csv << vertex.vertexId << ",";
        csv << vertex.kmerId << ",";
        csv << vertex.coverage << ",";
        std::ranges::copy(kmer, ostream_iterator<Base>(csv));
        csv << ",";
        csv << "\n";
    }
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




// This uses the assembly path to assemble sequence.
// The assembly path begins at vA, which contains invalid sequence (Base::fromInteger(10)),
// and ends at vB, which also contains invalid sequence (Base::fromInteger(20)).
// If the assembly path has N vertices, the base offset between vA and vB is N-1.
// The true assembled sequence begins k bases after vA and ends at the base before vB,
// so its length is n = N-1-k.

// For each assembled base position, we can choose among k vertices the one
// from which we get the base, which will be the same for all of the k vertices.
void LocalAssembly7::assemble(uint64_t k, Graph& graph)
{

    const uint64_t N = graph.assemblyPath.size();
    const uint64_t n = N - 1 - k;
    sequence.clear();
    sequence.resize(n);
    if(html) {
        html << "<br>The assembly path contains " << N << " vertices.";
        html << "<br>Assembled sequence is " << n << " bases long.";
    }



    // Look over all positions of assembled sequence.
    for(uint64_t position=0; position<n; position++) {

        // Loop over the k vertices we can use to get the
        // base at this position.
        Base b;
        for(uint64_t i=0; i<k; i++){
            const uint64_t positionInPath = position + k - i;
            const uint64_t positionInVertexKmer = i;
            const vertex_descriptor v = graph.assemblyPath[positionInPath];
            const Vertex& vertex = graph[v];
            const Kmer& kmer = kmers[vertex.kmerId].first;
            if(i == 0) {
                b = kmer[positionInVertexKmer];
            } else {
                SHASTA2_ASSERT(b == kmer[positionInVertexKmer]);
            }
        }
        sequence[position] = b;
    }
}



void LocalAssembly7::writeSequence() const
{
    if(not html) {
        return;
    }

    // Write the consensus.
    html <<
        "<h4>Assembled sequence</h4>"
        "<table>"
        "<tr><th class=left>Length<td class=left>" << sequence.size() <<
        "<tr><th class=left>Sequence<td class=left style='font-family:monospace;white-space:nowrap'>";
    std::ranges::copy(sequence, ostream_iterator<Base>(html));
    html << "</table>";

}



void LocalAssembly7::removeOutliers()
{
    // Gather the offsets.
    vector<pair<uint64_t, uint64_t > > offsetTable; // (index in orientedReadInfos, offset).
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const OrientedRead& orientedRead = orientedReads[i];
        if(orientedRead.isOnBothAnchors()) {
            const uint64_t offset = orientedRead.positionOffsetAB();
            offsetTable.push_back({i, offset});
        }
    }
    sort(offsetTable.begin(), offsetTable.end(), OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Find places where there is an unreasonably jump in the offset.
    vector<uint64_t> violations(1, 0);
    for(uint64_t i1=1; i1<offsetTable.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t offset0 = offsetTable[i0].second;
        const uint64_t offset1 = offsetTable[i1].second;
        if(not checkOffsets(offset0, offset1)) {
            violations.push_back(i1);
            // cout << "Violation " << offset0 << " " << offset1 << " " << i1 << endl;
        }
    }
    violations.push_back(offsetTable.size());

    // If no violations were found, keep all the OrientedReadInfos.
    // This is the most common case.
    if(violations.size() == 2) {
        return;
    }

#if 0
    cout << "violations vector ";
    std::ranges::copy(violations, ostream_iterator<uint64_t>(cout, " "));
    cout << endl;
#endif

    // Find the largest interval between violations.
    uint64_t keepBegin = 0;
    uint64_t keepEnd = 0;
    for(uint64_t i1=1; i1<violations.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t violation0 = violations[i0];
        const uint64_t violation1 = violations[i1];
        if(violation1 - violation0 > keepEnd - keepBegin) {
            keepBegin = violation0;
            keepEnd = violation1;
        }
    }
    // cout << "keepBegin " << keepBegin << ", keepEnd " << keepEnd << endl;

    // Only keep OrientedReads that are at positions [keepBegin, keepEnd)
    // in the offset table.
    std::set<uint64_t> discard;
    for(uint64_t i=0; i<keepBegin; i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReads[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }
    for(uint64_t i=keepEnd; i<offsetTable.size(); i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReads[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }

    vector<OrientedRead> newOrientedReads;
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        if(not discard.contains(i)) {
            newOrientedReads.push_back(orientedReads[i]);
        }
    }
    orientedReads.swap(newOrientedReads);
}



bool LocalAssembly7::checkOffsets(uint64_t offset0, uint64_t offset1)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double aDrift = 0.02;
    const double bDrift = 100.;

    if(offset1 == offset0) {
        return true;
    }

    SHASTA2_ASSERT(offset1 > offset0);

    const double average = 0.5 * double(offset0 + offset1);
    const uint64_t difference = offset1 - offset0;

    const double acceptableDifference = aDrift * average + bDrift;

    return difference < uint64_t(std::round(acceptableDifference));

}



void LocalAssembly7::writeKmerOccurrences(const Graph& graph, const string& fileName) const
{
    ofstream csv(fileName);
    writeKmerOccurrences(graph, csv);
}



void LocalAssembly7::writeKmerOccurrences(const Graph& graph, ostream& csv) const
{
    csv << "VertexId,Coverage,SequenceId,Position in sequence,Offset from left,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        for(const KmerOccurrence& kmerOccurrence: vertex.occurrences) {
            const SequenceInfo& sequence = sequences[kmerOccurrence.sequenceId];
            const uint64_t positionInSequence = kmerOccurrence.position;

            // Estimate the offset relative to the left anchor.
            int64_t offsetFromLeft;
            if(sequence.isOnAnchorA) {
                offsetFromLeft = positionInSequence;
            } else {
                const uint64_t sequenceLength = sequence.sequence.size();
                offsetFromLeft = int64_t(positionInSequence) + (int64_t(offset) - int64_t(sequenceLength));
            }

            csv << v << ",";
            csv << vertex.coverage << ",";
            csv << kmerOccurrence.sequenceId << ",";
            csv << positionInSequence << ",";
            csv << offsetFromLeft << ",";
            csv << "\n";
        }
    }
}



void LocalAssembly7::Graph::merge()
{
    Graph& graph = *this;
    const bool debug = false;

    // Check that all vertices have in-degree and out-degree
    // no greater than their number of occurrences.
    // This also implies that all vertices with one occurrence
    // cannot have more than one parent or children.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const uint64_t n = graph[v].occurrences.size();
        SHASTA2_ASSERT(in_degree(v, graph) <= n);
        SHASTA2_ASSERT(out_degree(v, graph) <= n);
    }


    // A stack of vertices that have one occurrence and that have
    // two or more children with one occurrence with the same kmerId.
    std::stack<vertex_descriptor> s;
    vector< vector<vertex_descriptor> > groups;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        findMergeableChildrenGroups(v, groups);
        if(not groups.empty()) {
            s.push(v);
            if(debug) {
                cout << "Added " << v << " to merge stack." << endl;
            }
        }
    }



    // Merge until the stack is empty.
    vector< vector<vertex_descriptor> > newGroups;
    while(not s.empty()) {

        // Dequeue a vertex that has children to be merged.
        const vertex_descriptor v = s.top();
        s.pop();
        if(debug) {
            cout << "Dequeued " << v << endl;
        }

        // Find the groups of children that should be merged.
        findMergeableChildrenGroups(v, groups);
        if(debug) {
            cout << "Found " << groups.size() <<  " groups." << endl;
        }

        // Merge each group.
        for(const vector<vertex_descriptor>& group: groups) {
            SHASTA2_ASSERT(group.size() > 1);
            const vertex_descriptor vNew = mergeGroup(group);
            if(debug) {
                cout << "Merged a group of " << group.size() << " vertices into new vertex " << vNew << endl;

                cout << "In-edges of new vertex:" << endl;
                BGL_FORALL_INEDGES(vNew, e, graph, Graph) {
                    cout << "From " << source(e, graph) << ", coverage " << graph[e].coverage << endl;
                }
                cout << "Out-edges of new vertex:" << endl;
                BGL_FORALL_OUTEDGES(vNew, e, graph, Graph) {
                    cout << "To " << target(e, graph) << ", coverage " << graph[e].coverage << endl;
                }
            }

            // See if the new vertex should be added to the stack.
            findMergeableChildrenGroups(vNew, newGroups);
            if(not newGroups.empty()) {
                s.push(vNew);
                if(debug) {
                    cout << "Enqueued " << vNew << ", stack size " << s.size() << endl;
                }
            }
        }
    }
}



void LocalAssembly7::Graph::findMergeableChildrenGroups(
    vertex_descriptor v0,
    vector< vector<vertex_descriptor> >& groups
    ) const
{
    const Graph& graph = *this;
    groups.clear();

    if(graph[v0].occurrences.size() < 2) {
        return;
    }
    if(out_degree(v0, graph) < 2) {
        return;
    }

    // Map with key = kmerId, value = children with that kmerId.
    // This can be made faster.
    std::map<uint64_t, vector<vertex_descriptor> > m;
    BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
        const vertex_descriptor v1 = target(e, graph);
        const Vertex& vertex1 = graph[v1];
        if(vertex1.occurrences.size() == 1) {
            m[vertex1.kmerId].push_back(v1);
        }
    }

    for(const auto& [kmerId, v]: m) {
        if(v.size() > 1) {
            groups.push_back(v);
        }
    }
}



LocalAssembly7::Graph::vertex_descriptor
    LocalAssembly7::Graph::mergeGroup(const vector<vertex_descriptor>& group)
{
    Graph& graph = *this;

    // Check that they all have the same kmerId,
    // one occurrence, and
    // their in-degree and out-degree are not greater than 1.
    uint64_t kmerId = invalid<uint64_t>;
    vector<Base> kmer;
    uint64_t coverage = 0;
    vector<KmerOccurrence> occurrences;
    std::map<vertex_descriptor, uint64_t> parentsWithEdgeCoverage;
    std::map<vertex_descriptor, uint64_t> childrenWithEdgeCoverage;
    for(const vertex_descriptor v: group) {
        const Vertex& vertex = graph[v];
        SHASTA2_ASSERT(vertex.occurrences.size() == 1);
        occurrences.push_back(vertex.occurrences.front());
        if(kmerId == invalid<uint64_t>) {
            kmerId = vertex.kmerId;
        } else {
            SHASTA2_ASSERT(kmerId == vertex.kmerId);
        }
        SHASTA2_ASSERT(in_degree(v, graph) < 2);
        SHASTA2_ASSERT(out_degree(v, graph) < 2);
        coverage += vertex.coverage;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const vertex_descriptor child = target(e, graph);
            const uint64_t coverage = graph[e].coverage;
            const auto it = childrenWithEdgeCoverage.find(child);
            if(it == childrenWithEdgeCoverage.end()) {
                childrenWithEdgeCoverage.insert(make_pair(child, coverage));
            } else {
                it->second += coverage;
            }
        }
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            const vertex_descriptor parent = source(e, graph);
            const uint64_t coverage = graph[e].coverage;
            const auto it = parentsWithEdgeCoverage.find(parent);
            if(it == parentsWithEdgeCoverage.end()) {
                parentsWithEdgeCoverage.insert(make_pair(parent, coverage));
            } else {
                it->second += coverage;
            }
        }
    }

    const vertex_descriptor vNew = add_vertex(Vertex(graph.nextVertexId++, kmerId, occurrences, coverage), graph);
    for(const auto& [child, coverage]: childrenWithEdgeCoverage) {
        auto[e, ignore] = boost::add_edge(vNew, child, graph);
        graph[e].coverage = coverage;
    }
    for(const auto& [parent, coverage]: parentsWithEdgeCoverage) {
        auto[e, ignore] = boost::add_edge(parent, vNew, graph);
        graph[e].coverage = coverage;
    }


    // Disconnect all the vertices in the group we merged.
    for(const vertex_descriptor v: group) {
        boost::clear_vertex(v, graph);
    }

    return vNew;
}



void LocalAssembly7::runAbpoa()
{
    throw runtime_error("LocalAssembly7::runAbpoa not implemented.");
}



void LocalAssembly7::runPoasta()
{
    throw runtime_error("LocalAssembly7::runPoasta not implemented.");
}



void LocalAssembly7::runTheseus()
{
    throw runtime_error("LocalAssembly7::runTheseus not implemented.");
}



void LocalAssembly7::runAdaptive()
{
    throw runtime_error("LocalAssembly7::runAdaptive not implemented.");
}
