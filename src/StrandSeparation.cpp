// Shasta2.
#include "StrandSeparation.hpp"
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "html.hpp"
using namespace shasta2;
using namespace StrandSeparation;

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"



StrandContact::StrandContact(
    AssemblyGraph& assemblyGraph,
    const vector<Segment>& allSegments, // All Segments,sorted by id.
    const string& debugOutputBaseName,  // Only used for debug output.
    uint64_t strandContactId            // Only used for debug output.
    ) :
    assemblyGraph(assemblyGraph),
    debugOutputBaseName(debugOutputBaseName),
    strandContactId(strandContactId)
{
    const bool debug = true;
    if(debug) {
        html.open(debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + ".html");
        cout << "Working on strand contact " << strandContactId << endl;
        writeHtmlBegin(html, "Strand contact " + to_string(strandContactId));
        html << "<h1>Strand contact " << strandContactId << "</h1>";
        writeAllSegments(allSegments);
    }

    gatherSegments(allSegments);
    countReadOccurrences();
    createBipartiteGraph();

    vector< vector<uint64_t> > componentsIndexes;
    bipartiteGraph.strandSeparation(componentsIndexes);
    writeBipartiteGraph(componentsIndexes);

    if(debug) {
        writeHtmlEnd(html);
        cout << "Done working on strand contact " << strandContactId << endl;
    }
}



uint64_t StrandContact::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandContact::writeAllSegments(const vector<Segment>& allSegments)
{
    if(not html) {
        return;
    }

    html << "<h2>All segments, sorted by id</h2>";
    for(uint64_t i=0; i<allSegments.size(); i++) {
        if(i != 0) {
            html << ",<wbr>";
        }
        html << id(allSegments[i]);
    }
}



void StrandContact::gatherSegments(const vector<Segment>& allSegments)
{
    // Sanity check: allSegments must be sorted by id.
    SHASTA2_ASSERT(std::is_sorted(allSegments.begin(), allSegments.end(), assemblyGraph.orderById));

    // Sanity check: if allSegments contains segment0, it must
    // also contain its reverse complement, segment1.
    for(const Segment segment0: allSegments) {
        const Segment segment1 = assemblyGraph[segment0].eRc;
        SHASTA2_ASSERT(segment0 != segment1);
        SHASTA2_ASSERT(std::binary_search(allSegments.begin(), allSegments.end(),
            segment1, assemblyGraph.orderById));
    }

    // Gather the SegmentPairs.
    for(const Segment segment0: allSegments) {
        const Segment segment1 = assemblyGraph[segment0].eRc;
        if(id(segment0) < id(segment1)) {
            allSegmentPairs.push_back({segment0, segment1});
        }
    }

    // Gather the lowCoverageSegmentPairs and the highCoverageSegmentPairs.
    for(const SegmentPair& segmentPair: allSegmentPairs) {
        const AssemblyGraphEdge& edge0 = assemblyGraph[segmentPair.segment0];
        const AssemblyGraphEdge& edge1 = assemblyGraph[segmentPair.segment1];
        const double coverage0 = edge0.lengthWeightedAverageCoverage();
        const double coverage1 = edge1.lengthWeightedAverageCoverage();
        SHASTA2_ASSERT(coverage0 == coverage1);
        if(coverage0 <= coverageThreshold) {
            lowCoverageSegmentPairs.push_back(segmentPair);
        } else {
            highCoverageSegmentPairs.push_back(segmentPair);
        }
    }

    // Fill in the lowCoverageSegments vector.
    for(const SegmentPair& segmentPair: lowCoverageSegmentPairs) {
        lowCoverageSegments.push_back(segmentPair.segment0);
        lowCoverageSegments.push_back(segmentPair.segment1);
    }

    writeSegments();
}



void StrandContact::writeSegments()
{
    if(not html) {
        return;
    }

    html << "<h2>All segment pairs</h2>";
    for(uint64_t i=0; i<allSegmentPairs.size(); i++) {
        const SegmentPair& segmentPair = allSegmentPairs[i];
        if(i != 0) {
            html << ",<wbr>";
        }
        html << id(segmentPair.segment0) << ",<wbr>" << id(segmentPair.segment1);
    }

    html << "<h2>Low coverage segment pairs</h2>";
    for(uint64_t i=0; i<lowCoverageSegmentPairs.size(); i++) {
        const SegmentPair& segmentPair = lowCoverageSegmentPairs[i];
        if(i != 0) {
            html << ",<wbr>";
        }
        html << id(segmentPair.segment0) << ",<wbr>" << id(segmentPair.segment1);
    }

    html << "<h2>High coverage segment pairs</h2>";
    for(uint64_t i=0; i<highCoverageSegmentPairs.size(); i++) {
        const SegmentPair& segmentPair = highCoverageSegmentPairs[i];
        if(i != 0) {
            html << ",<wbr>";
        }
        html << id(segmentPair.segment0) << ",<wbr>" << id(segmentPair.segment1);
    }
}



void StrandContact::countReadOccurrences()
{
    // Count occurrences of reads in the first Segment of each low coverage Segment pair.
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

    // Remove from the map reads that occur in just one Segment.
    std::map<ReadId, vector<ReadOccurrence> > newReadOccurrenceMap;
    for(const auto&p: readOccurrenceMap) {
        if(p.second.size() > 1) {
            newReadOccurrenceMap.insert(p);
        }
    }
    newReadOccurrenceMap.swap(readOccurrenceMap);

    writeReadOccurrences();
}



void StrandContact::writeReadOccurrences()
{
    if(not html) {
        return;
    }

    const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-ReadOccurrences.csv";
    ofstream csv(fileName);
    csv << "ReadId,Strand,Segment,Frequency\n";
    for(const auto&[readId, occurrences]: readOccurrenceMap) {
        for(const auto& occurrence: occurrences) {
            const Segment segment = lowCoverageSegmentPairs[occurrence.segmentPairIndex].segment0;
            csv << readId << ",";
            csv << occurrence.strand << ",";
            csv << id(segment) << ",";
            csv << occurrence.frequency << "\n";
        }
    }

    html << "<h2>Read occurrences</h2>"
        "For details of occurrences of reads in the low coverage segments, see "
        "<a href='"<< fileName << "'>" << fileName << "</a>.";
}



void StrandContact::createBipartiteGraph()
{
    // Generate the vertices corresponding to low coverage segments.
    bipartiteGraph.segmentIndexToVertexMap.resize(lowCoverageSegments.size());
    for(uint64_t segmentIndex=0; segmentIndex<lowCoverageSegments.size(); segmentIndex++) {
        const BipartiteGraph::vertex_descriptor v =
            boost::add_vertex(BipartiteGraphVertex(segmentIndex), bipartiteGraph);
        bipartiteGraph.segmentIndexToVertexMap[segmentIndex] = v;
    }

    // Generate the vertices corresponding to OrientedReadIds.
    for(const auto&[readId, ignore]: readOccurrenceMap) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const BipartiteGraph::vertex_descriptor v =
                boost::add_vertex(BipartiteGraphVertex(orientedReadId), bipartiteGraph);
            bipartiteGraph.orientedReadIdToVertexMap.insert({orientedReadId, v});
        }
    }

    // Map vertices t indexes and vice versa.
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        bipartiteGraph.vertexIndexMap.insert({v, vertexIndex++});
        bipartiteGraph.vertexTable.push_back(v);
    }



    // Now generate the edges.
    // Each ReadOccurrence generates a pair of reverse complemented edges.
    for(const auto&[readId, occurrences]: readOccurrenceMap) {
        for(const auto& occurrence: occurrences) {

            const OrientedReadId orientedReadIdA(readId, occurrence.strand);;
            const BipartiteGraph::vertex_descriptor vOrientedReadA =
                bipartiteGraph.orientedReadIdToVertexMap.at(orientedReadIdA);
            const uint64_t segmentIndexA = 2 * occurrence.segmentPairIndex;
            const BipartiteGraph::vertex_descriptor vSegmentA =
                bipartiteGraph.segmentIndexToVertexMap[segmentIndexA];
            auto[e, ignore] = boost::add_edge(vOrientedReadA, vSegmentA,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);


            const OrientedReadId orientedReadIdB(readId, 1 - occurrence.strand);;
            const BipartiteGraph::vertex_descriptor vOrientedReadB =
                bipartiteGraph.orientedReadIdToVertexMap.at(orientedReadIdB);
            const uint64_t segmentIndexB = 2 * occurrence.segmentPairIndex + 1;
            const BipartiteGraph::vertex_descriptor vSegmentB =
                bipartiteGraph.segmentIndexToVertexMap[segmentIndexB];
            auto [eRc, ignoreRc] = boost::add_edge(vOrientedReadB, vSegmentB,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);

            bipartiteGraph.edgePairs.push_back({e, eRc, occurrence.frequency});
        }
    }
    sort(bipartiteGraph.edgePairs.begin(), bipartiteGraph.edgePairs.end());
}



void BipartiteGraph::writeGraphviz(
    const string& fileName,
    const vector<Segment>& lowCoverageSegments,
    const vector< vector<uint64_t> >& componentsIndexes,
    const AssemblyGraph& assemblyGraph) const
{
    const BipartiteGraph& bipartiteGraph = *this;

    // Find the component that each vertex belongs to.
    vector<uint64_t> componentTable(vertexTable.size(), invalid<uint64_t>);
    for(uint64_t componentId=0; componentId<componentsIndexes.size(); componentId++) {
        const vector<uint64_t>& componentIndexes = componentsIndexes[componentId];
        for(const uint64_t vertexIndex: componentIndexes) {
            componentTable[vertexIndex] = componentId;
        }
    }

    ofstream dot(fileName);

    dot << "graph BipartiteGraph {\n";

    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex = bipartiteGraph[v];
        const Segment segment = lowCoverageSegments[vertex.segmentIndex];
        const uint64_t vertexIndex = vertexIndexMap.at(v);
        const uint64_t componentId = componentTable[vertexIndex];
        const string color = randomHslColor(componentId, 0.75, 0.5);
        SHASTA2_ASSERT(componentId != invalid<uint64_t>);
        if(vertex.isSegment) {
            dot << assemblyGraph.id(segment);
            dot << "[width=0.1";
        } else {
            dot << "\"" << vertex.orientedReadId << "\"";
            dot << "[width=0.02";
        }
        dot << " color=\"" << color << "\"";
        dot << "]";
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const vertex_descriptor v0 = source(e, bipartiteGraph);
        const vertex_descriptor v1 = target(e, bipartiteGraph);
        const BipartiteGraphVertex& vertex0 = bipartiteGraph[v0];
        const BipartiteGraphVertex& vertex1 = bipartiteGraph[v1];

        if(vertex0.isSegment) {
            const Segment segment0 = lowCoverageSegments[vertex0.segmentIndex];
            dot << assemblyGraph.id(segment0);
        } else {
            dot << "\"" << vertex0.orientedReadId << "\"";
        }
        dot << "--";

        if(vertex1.isSegment) {
            const Segment segment1 = lowCoverageSegments[vertex1.segmentIndex];
            dot << assemblyGraph.id(segment1);
        } else {
            dot << "\"" << vertex1.orientedReadId << "\"";
        }

        dot << "[";
        dot << "tooltip=\"";
        if(vertex0.isSegment) {
            const Segment segment0 = lowCoverageSegments[vertex0.segmentIndex];
            dot << assemblyGraph.id(segment0);
        } else {
            dot << vertex0.orientedReadId;
        }
        dot << " ";
        if(vertex1.isSegment) {
            const Segment segment1 = lowCoverageSegments[vertex1.segmentIndex];
            dot << assemblyGraph.id(segment1);
        } else {
            dot << vertex1.orientedReadId;
        }
        dot << " " << bipartiteGraph[e].frequency;
        dot << "\"";

        dot << "]";

        dot << ";\n";

    }

    dot << "}\n";
}



void StrandContact::writeBipartiteGraph(const vector< vector<uint64_t> >& componentsIndexes)
{
    if(not html) {
        return;
    }

    const string dotFileName = debugOutputBaseName + "-StrandContact-" +
        to_string(strandContactId) + ".dot";
    bipartiteGraph.writeGraphviz(dotFileName, lowCoverageSegments, componentsIndexes, assemblyGraph);

    const double timeout = 30.;
    const string options = "-Nshape=point -Epenwidth=0.2 -Gratio=expand -Gsize=15";
    html << "<h2>Bipartite graph</h2>"
        "<br>In the bipartite graph, each vertex represents a low coverage segment or "
        "an oriented read. Oriented reads are displayed as small dots. "
        "<br><br>" << dotFileName << "<br>";

    try {
        graphvizToHtml(dotFileName, "sfdp", timeout, options, html, true);
    } catch (std::exception&) {
        html << "The bipartite graph is too complex to display.";
    }
}



BipartiteGraph::vertex_descriptor BipartiteGraph::reverseComplement(vertex_descriptor v) const
{
    const BipartiteGraph& bipartiteGraph = *this;
    const auto& vertex = bipartiteGraph[v];

    if(vertex.isSegment) {
        const uint64_t segmentIndex = vertex.segmentIndex;
        const uint64_t segmentIndexRc = segmentIndex ^ 1;
        return segmentIndexToVertexMap[segmentIndexRc];
    } else {
        OrientedReadId orientedReadId = vertex.orientedReadId;
        orientedReadId.flipStrand();
        return orientedReadIdToVertexMap.at(orientedReadId);
    }
}




// One attempt at strand separation.
// This add edge pairs in the given order to a disjoint sets
// data structure. However an edge pair is discarded if it
// would cause contacts between strands.
// This leaves the BipartiteGraph unchanged and returns
// the connected components computed in this way,
// sorted by decreasing size. Reverse complemented
// connected components are not guaranteed to be consecutive.
// Vertex indexes in each component are sorted.
// They can be converted to vertex_descriptors
// via the vertexTable.
void BipartiteGraph::strandSeparation(
    const vector<uint64_t>& edgePairsIndexes,
    vector< vector<uint64_t> >& componentsIndexes) const
{
    const BipartiteGraph& bipartiteGraph = *this;

    DisjointSets disjointSets(vertexIndexMap.size());

    for(const uint64_t edgePairIndex: edgePairsIndexes) {
        const auto& edgePair = edgePairs[edgePairIndex];
        const auto eA = edgePair.e;
        const auto eB = edgePair.eRc;

        const auto v0A = source(eA, bipartiteGraph);
        const auto v1A = target(eA, bipartiteGraph);
        const auto v0B = source(eB, bipartiteGraph);
        const auto v1B = target(eB, bipartiteGraph);

        const auto v0ARc = reverseComplement(v0A);
        const auto v1ARc = reverseComplement(v1A);
        const auto v0BRc = reverseComplement(v0B);
        const auto v1BRc = reverseComplement(v1B);

        const uint64_t i0A = vertexIndexMap.at(v0A);
        const uint64_t i1A = vertexIndexMap.at(v1A);
        const uint64_t i0B = vertexIndexMap.at(v0B);
        const uint64_t i1B = vertexIndexMap.at(v1B);

        const uint64_t i0ARc = vertexIndexMap.at(v0ARc);
        const uint64_t i1ARc = vertexIndexMap.at(v1ARc);
        const uint64_t i0BRc = vertexIndexMap.at(v0BRc);
        const uint64_t i1BRc = vertexIndexMap.at(v1BRc);

        const bool strandViolationA = (disjointSets.findSet(i1A) == disjointSets.findSet(i0ARc));
        const bool strandViolationB = (disjointSets.findSet(i1B) == disjointSets.findSet(i0BRc));
        const bool strandViolationARc = (disjointSets.findSet(i0A) == disjointSets.findSet(i1ARc));
        const bool strandViolationBRc = (disjointSets.findSet(i0B) == disjointSets.findSet(i1BRc));

        const bool strandViolation = strandViolationA;
        SHASTA2_ASSERT(strandViolationB == strandViolation);
        SHASTA2_ASSERT(strandViolationARc == strandViolation);
        SHASTA2_ASSERT(strandViolationBRc == strandViolation);

        if(not strandViolation) {
            disjointSets.unionSet(i0A, i1A);
            disjointSets.unionSet(i0B, i1B);
        }
    }


    // Gather the components.
    disjointSets.gatherComponents(1, componentsIndexes);
}



void BipartiteGraph::strandSeparation(vector< vector<uint64_t> >& componentsIndexes)
{
    // For now just do one attempt, using the EdgePairs in
    // the order in which they are stored, that is,
    // by decreasing frequency.
    vector<uint64_t> edgePairIndexes(edgePairs.size());
    std::ranges::iota(edgePairIndexes, 0);

    strandSeparation(edgePairIndexes, componentsIndexes);


    // Reorder the componentsIndexes so pairs of
    // reverse complemented components are numbered consecutively.
    // Because of the way vertices are created, this can be done simply by
    // sorting by the first index of each component index.
    class SortByFirstElement {
    public:
    public:
         bool operator()(const vector<uint64_t>& x, const vector<uint64_t>& y) const
        {
             SHASTA2_ASSERT(not x.empty());
             SHASTA2_ASSERT(not y.empty());
             return x.front() < y.front();
        }
    };
    sort(componentsIndexes.begin(), componentsIndexes.end(), SortByFirstElement());
}
