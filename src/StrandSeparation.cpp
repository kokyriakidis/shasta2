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
#include <iomanip>
#include "iostream.hpp"



StrandContact::StrandContact(
    AssemblyGraph& assemblyGraph,
    const vector<Segment>& allSegmentsById, // All Segments,sorted by id.
    const string& debugOutputBaseName,  // Only used for debug output.
    uint64_t strandContactId            // Only used for debug output.
    ) :
    assemblyGraph(assemblyGraph),
    allSegmentsById(allSegmentsById),
    debugOutputBaseName(debugOutputBaseName),
    strandContactId(strandContactId)
{
    const bool debug = true;
    if(debug) {
        html.open(debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + ".html");
        cout << "Working on strand contact " << strandContactId << endl;
        writeHtmlBegin(html, "Strand contact " + to_string(strandContactId));
        html << "<h1>Strand contact " << strandContactId << "</h1>";
        writeAllSegmentsById();
    }

    gatherSegments();
    countReadOccurrences();
    createBipartiteGraph();

    vector< vector<uint64_t> > componentsIndexes;
    bipartiteGraph.strandSeparation(componentsIndexes);
    writeBipartiteGraph(componentsIndexes);
    writeComponents(componentsIndexes);
    SHASTA2_ASSERT(not componentsIndexes.empty());
    gatherStrands(componentsIndexes);
    classifySegments();

    if(debug) {
        writeHtmlEnd(html);
        cout << "Done working on strand contact " << strandContactId << endl;
    }
}



uint64_t StrandContact::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandContact::writeAllSegmentsById()
{
    if(not html) {
        return;
    }

    html << "<h2>All segments, sorted by id</h2>";
    for(uint64_t i=0; i<allSegmentsById.size(); i++) {
        if(i != 0) {
            html << ",<wbr>";
        }
        html << id(allSegmentsById[i]);
    }
}



void StrandContact::gatherSegments()
{
    SHASTA2_ASSERT(std::is_sorted(allSegmentsById.begin(), allSegmentsById.end(), assemblyGraph.orderById));

    // Sanity check: if allSegments contains segment0, it must
    // also contain its reverse complement, segment1.
    for(const Segment segment0: allSegmentsById) {
        const Segment segment1 = assemblyGraph[segment0].eRc;
        SHASTA2_ASSERT(segment0 != segment1);
        SHASTA2_ASSERT(std::binary_search(allSegmentsById.begin(), allSegmentsById.end(),
            segment1, assemblyGraph.orderById));
    }

    // Gather the SegmentPairs.
    for(const Segment segment0: allSegmentsById) {
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
        if(coverage0 != coverage1) {
            throw runtime_error("Coverage check failed at segments " +
                to_string(edge0.id) + " " + to_string(edge1.id) + ".");
        }
        if(coverage0 <= coverageThreshold) {
            lowCoverageSegmentPairs.push_back(segmentPair);
        } else {
            highCoverageSegmentPairs.push_back(segmentPair);
        }
    }

    // Fill in the allSegments vector.
    for(const SegmentPair& segmentPair: allSegmentPairs) {
        allSegments.push_back(segmentPair.segment0);
        allSegments.push_back(segmentPair.segment1);
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

    /*
    // Remove from the map reads that occur in just one Segment.
    std::map<ReadId, vector<ReadOccurrence> > newReadOccurrenceMap;
    for(const auto&p: readOccurrenceMap) {
        if(p.second.size() > 1) {
            newReadOccurrenceMap.insert(p);
        }
    }
    newReadOccurrenceMap.swap(readOccurrenceMap);
    */

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

        string color;
        if(componentId == 0) {
            color = StrandContact::color(StrandContact::SegmentClassification::LowCoverageStrand0);
        } else if(componentId == 1) {
            color = StrandContact::color(StrandContact::SegmentClassification::LowCoverageStrand1);
        } else {
            color = randomHslColor(componentId, 0.75, 0.5);
        }

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
        const uint64_t vertexIndex0 = vertexIndexMap.at(v0);
        const uint64_t vertexIndex1 = vertexIndexMap.at(v1);
        const BipartiteGraphVertex& vertex0 = bipartiteGraph[v0];
        const BipartiteGraphVertex& vertex1 = bipartiteGraph[v1];

        const uint64_t componentId0 = componentTable[vertexIndex0];
        const uint64_t componentId1 = componentTable[vertexIndex1];

        const uint64_t frequency = bipartiteGraph[e].frequency;
        const double thickness = 0.3 * (1. + std::log10(frequency));

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

        dot << " penwidth=\"" << thickness << "\"";

        if(componentId0 != componentId1) {
            dot << " color=red";
        }

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



void StrandContact::writeComponents(const vector< vector<uint64_t> >& componentsIndexes)
{
    if(not html) {
        return;
    }

    for(uint64_t componentId=0; componentId<componentsIndexes.size(); componentId++) {
        const vector<uint64_t>& componentIndexes = componentsIndexes[componentId];

        html << "<h2>Component "<< componentId << "</h2>";

        // Segments.
        bool isFirstTime = true;
        for(uint64_t i=0; i<componentIndexes.size(); i++) {
            const uint64_t vertexIndex = componentIndexes[i];
            const BipartiteGraph::vertex_descriptor v = bipartiteGraph.vertexTable[vertexIndex];
            const BipartiteGraphVertex& vertex = bipartiteGraph[v];
            if(vertex.isSegment) {
                const Segment segment = lowCoverageSegments[vertex.segmentIndex];
                if(isFirstTime) {
                    isFirstTime = false;
                } else {
                    html << ",<wbr>";
                }
                html << id(segment);
            }
        }

        // Oriented reads.
        html << "<br><br>";
        isFirstTime = true;
        for(uint64_t i=0; i<componentIndexes.size(); i++) {
            const uint64_t vertexIndex = componentIndexes[i];
            const BipartiteGraph::vertex_descriptor v = bipartiteGraph.vertexTable[vertexIndex];
            const BipartiteGraphVertex& vertex = bipartiteGraph[v];
            if(not vertex.isSegment) {
                if(isFirstTime) {
                    isFirstTime = false;
                } else {
                    html << ",<wbr>";
                }
                html << vertex.orientedReadId;
            }
        }
    }
}



void StrandContact::gatherStrands(const vector< vector<uint64_t> >& componentsIndexes)
{
    SHASTA2_ASSERT(componentsIndexes.size() >= 2);

    for(uint64_t strand=0; strand<2; strand++) {
        const vector<uint64_t>& componentIndexes = componentsIndexes[strand];
        for(const uint64_t vertexIndex: componentIndexes) {
            const BipartiteGraph::vertex_descriptor v = bipartiteGraph.vertexTable[vertexIndex];
            const BipartiteGraphVertex& vertex = bipartiteGraph[v];
            if(vertex.isSegment) {
                const Segment segment = lowCoverageSegments[vertex.segmentIndex];
                strandSegments[strand].push_back(segment);
            } else {
                strandOrientedReadIds[strand].push_back(vertex.orientedReadId);
            }
        }
    }



    if(html) {
        for(uint64_t strand=0; strand<2; strand++) {
            html << "<h2>Strand " << strand << " low coverage segments</h2>";
            for(uint64_t i=0; i<strandSegments[strand].size(); i++) {
                if(i != 0) {
                    html << ",<wbr>";
                }
                html << id(strandSegments[strand][i]);
            }

            html << "<h2>Strand " << strand << " oriented reads</h2>";
            for(uint64_t i=0; i<strandOrientedReadIds[strand].size(); i++) {
                if(i != 0) {
                    html << ",<wbr>";
                }
                html << strandOrientedReadIds[strand][i];
            }

            sort(
                strandSegments[strand].begin(),
                strandSegments[strand].end(),
                assemblyGraph.orderById);
            sort(
                strandOrientedReadIds[strand].begin(),
                strandOrientedReadIds[strand].end());
        }
    }

}



void StrandContact::classifySegments()
{

    segmentClassifications.resize(allSegments.size(), SegmentClassification::Invalid);

    // Classify low coverage segments.
    for(const Segment segment: lowCoverageSegments) {
        const bool isStrand0 = std::binary_search(
            strandSegments[0].begin(),
            strandSegments[0].end(),
            segment, assemblyGraph.orderById);
        const bool isStrand1 = std::binary_search(
            strandSegments[1].begin(),
            strandSegments[1].end(),
            segment, assemblyGraph.orderById);
        const uint64_t indexInAllSegments = std::lower_bound(
            allSegmentsById.begin(),
            allSegmentsById.end(),
            segment,
            assemblyGraph.orderById) - allSegmentsById.begin();
        if(isStrand0) {
            SHASTA2_ASSERT(not isStrand1);
            segmentClassifications[indexInAllSegments] = SegmentClassification::LowCoverageStrand0;
        } else if(isStrand1) {
            segmentClassifications[indexInAllSegments] = SegmentClassification::LowCoverageStrand1;
        } else {
            segmentClassifications[indexInAllSegments] = SegmentClassification::LowCoverageUnclassified;
        }

    }



    // Classify high coverage segments.
    if(html) {
        html << "<h2>Classifying high coverage segments</h2>"
            "<table>"
            "<tr><th>Segment0<th>Segment1"
            "<th>n0<th>n1"
            "<th>Fraction0<th>Fraction1" <<
            std::fixed << std::setprecision(2);
    }

    for(const SegmentPair& segmentPair: highCoverageSegmentPairs) {
        const Segment segment0 = segmentPair.segment0;
        const Segment segment1 = segmentPair.segment1;
        const AssemblyGraphEdge& edge = assemblyGraph[segment0];
        uint64_t n0 = 0;
        uint64_t n1 = 0;
        for(const AssemblyGraphEdgeStep& step: edge) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                if(std::binary_search(
                    strandOrientedReadIds[0].begin(), strandOrientedReadIds[0].end(), orientedReadId)) {
                    ++n0;
                }
                if(std::binary_search(
                    strandOrientedReadIds[1].begin(), strandOrientedReadIds[1].end(), orientedReadId)) {
                    ++n1;
                }
            }
        }
        const double fraction0 = double(n0) / double(n0 + n1);
        const double fraction1 = 1. - fraction0;

        const bool isStrand0 = (fraction0 >= strandFractionThreshold);
        const bool isStrand1 = (fraction1 >= strandFractionThreshold);

        const uint64_t indexInAllSegments0 = std::lower_bound(
            allSegmentsById.begin(),
            allSegmentsById.end(),
            segment0,
            assemblyGraph.orderById) - allSegmentsById.begin();
        const uint64_t indexInAllSegments1 = std::lower_bound(
            allSegmentsById.begin(),
            allSegmentsById.end(),
            segment1,
            assemblyGraph.orderById) - allSegmentsById.begin();
        if(isStrand0) {
            SHASTA2_ASSERT(not isStrand1);
            segmentClassifications[indexInAllSegments0] = SegmentClassification::HighCoverageStrand0;
            segmentClassifications[indexInAllSegments1] = SegmentClassification::HighCoverageStrand1;
        } else if(isStrand1) {
            segmentClassifications[indexInAllSegments0] = SegmentClassification::HighCoverageStrand1;
            segmentClassifications[indexInAllSegments1] = SegmentClassification::HighCoverageStrand0;
        } else {
            segmentClassifications[indexInAllSegments0] = SegmentClassification::HighCoverageAmbiguous;
            segmentClassifications[indexInAllSegments1] = SegmentClassification::HighCoverageAmbiguous;
        }


        if(html) {
            html << "<tr>";

            html <<"<td class=centered";
            if(isStrand0) {
                html << " style='background-color:" << color(SegmentClassification::HighCoverageStrand0) << "'";
            }
            if(isStrand1) {
                html << " style='background-color:" << color(SegmentClassification::HighCoverageStrand1) << "'";
            }
            html << ">" << id(segment0);

            html << "<td class=centered";
            if(isStrand1) {
                html << " style='background-color:" << color(SegmentClassification::HighCoverageStrand0) << "'";
            }
            if(isStrand0) {
                html << " style='background-color:" << color(SegmentClassification::HighCoverageStrand1) << "'";
            }
            html << ">" << id(segment1);

            html <<
                "<td class=centered>" << n0 <<
                "<td class=centered>" << n1 <<
                "<td class=centered>" << fraction0 <<
                "<td class=centered>" << fraction1;
        }
    }
    html << "</table>";

    // Check that all segments have a valid classification.
    for(const auto& segmentClassification: segmentClassifications) {
        SHASTA2_ASSERT(segmentClassification != SegmentClassification::Invalid);
    }



    // Write a csv file containing the file of each Segment which can be loaded in Bandage.
    if(html) {
        const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-Bandage.csv";
        ofstream csv(fileName);
        csv << "Segment,Classification,Color\n";
        for(uint64_t i=0; i<allSegmentsById.size(); i++) {
            const auto segmentClassification = segmentClassifications[i];

            string type = "Nothing";

            csv << id(allSegmentsById[i]) << ",";
            csv << type << ",";
            csv << color(segmentClassification) << "\n";
        }

    }
}



string StrandContact::color(SegmentClassification segmentClassification)
{
    // The colorTable is filled in at the first call.
    static array<string, static_cast<uint64_t>(SegmentClassification::MaxValue)> colorTable;

    static bool isFirstTime = true;
    if(isFirstTime) {
        isFirstTime = false;


        // Strand 0 uses hue=0.6 (blue).
        colorTable[static_cast<uint64_t>(SegmentClassification::LowCoverageStrand0)] =
            hslToRgbString(0.6, 1., 0.5);
        colorTable[static_cast<uint64_t>(SegmentClassification::HighCoverageStrand0)] =
            hslToRgbString(0.6, 1, 0.75);

        // Strand 1 uses hue=0.1 (orange).
        colorTable[static_cast<uint64_t>(SegmentClassification::LowCoverageStrand1)] =
            hslToRgbString(0.1, 1., 0.5);
        colorTable[static_cast<uint64_t>(SegmentClassification::HighCoverageStrand1)] =
            hslToRgbString(0.1, 1, 0.75);

        // Unclassified/ambiguous uses hue=0.9 (purple).
        colorTable[static_cast<uint64_t>(SegmentClassification::LowCoverageUnclassified)] =
            hslToRgbString(0.9, 1., 0.5);
        colorTable[static_cast<uint64_t>(SegmentClassification::HighCoverageAmbiguous)] =
            hslToRgbString(0.9, 1, 0.75);

    }

    const uint64_t index = static_cast<uint64_t>(segmentClassification);
    SHASTA2_ASSERT(index < colorTable.size());
    return colorTable[index];
}
