// Shasta2.
#include "StrandSeparation1.hpp"
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "graphvizToHtml.hpp"
#include "html.hpp"
using namespace shasta2;
using namespace StrandSeparation1;

// Standard library.
#include <iomanip>


// StrandContact constructor.
// The strandContactVerticesmust be sorted by id.
// The debugOutputBaseName and strandContactId are only used for debug output.
StrandContact::StrandContact(
    AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraph::vertex_descriptor>& strandContactVertices,
    const string& debugOutputBaseName,
    uint64_t strandContactId
    ) :
    assemblyGraph(assemblyGraph),
    strandContactVertices(strandContactVertices),
    debugOutputBaseName(debugOutputBaseName),
    strandContactId(strandContactId)
{
    const bool debug = true;
    if(debug) {
        html.open(debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + ".html");
        cout << "Working on strand contact " << strandContactId << endl;
        writeHtmlBegin(html, "Strand contact " + to_string(strandContactId));
        writeMakeAllTablesCopyable(html);
        html <<
            "</head>"
            "<body onload='makeAllTablesCopyable()'>";
        html << "<h1>Strand contact " << strandContactId << "</h1>";
    }

    gatherSegmentPairs();
    countReadOccurrences();
    createBipartiteGraph();
    bipartiteGraph.strandSeparation();
    writeBipartiteGraph();
    evaluateStrandSeparation(true);

    if(debug) {
        writeHtmlEnd(html);
        cout << "Done working on strand contact " << strandContactId << endl;
    }
}


uint64_t StrandContact::id(AssemblyGraphBaseClass::vertex_descriptor v) const
{
    return assemblyGraph.id(v);
}



uint64_t StrandContact::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandContact::gatherSegmentPairs()
{
    SHASTA2_ASSERT(std::ranges::is_sorted(strandContactVertices, assemblyGraph.orderById));

    // Out-edges of the strandContactVertices give us Segments
    // internal to the StrandContact plus the exits.
    for(const AssemblyGraph::vertex_descriptor v0: strandContactVertices) {
        BGL_FORALL_OUTEDGES(v0, segment, assemblyGraph, AssemblyGraph) {
            const Segment segmentRc = assemblyGraph[segment].eRc;
            SHASTA2_ASSERT(segmentRc != segment);

            if(id(segment) < id(segmentRc)) {
                const AssemblyGraph::vertex_descriptor v1 = target(segment, assemblyGraph);
                const bool isExit = not std::ranges::binary_search(strandContactVertices, v1, assemblyGraph.orderById);
                SegmentPair& segmentPair = segmentPairs.emplace_back();

                SegmentInfo& segmentInfo0 = segmentPair.segmentInfos[0];
                SegmentInfo& segmentInfo1 = segmentPair.segmentInfos[1];

                segmentInfo0.segment = segment;
                segmentInfo0.id = id(segment);
                segmentInfo0.isExit = isExit;
                segmentInfo1.segment = segmentRc;
                segmentInfo1.id = id(segmentRc);
                segmentInfo1.isEntrance = isExit;

                segmentPair.length = assemblyGraph[segment].length();
                SHASTA2_ASSERT(segmentPair.length == assemblyGraph[segmentRc].length());

                segmentPair.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
                SHASTA2_ASSERT(segmentPair.coverage == assemblyGraph[segmentRc].lengthWeightedAverageCoverage());
            }
        }
    }



    // In-edges of the strandContactVertices give us the entrances.
    for(const AssemblyGraph::vertex_descriptor v0: strandContactVertices) {
        BGL_FORALL_INEDGES(v0, segment, assemblyGraph, AssemblyGraph) {
            const Segment segmentRc = assemblyGraph[segment].eRc;
            SHASTA2_ASSERT(segmentRc != segment);

            if(id(segment) < id(segmentRc)) {
                const AssemblyGraph::vertex_descriptor v1 = source(segment, assemblyGraph);
                const bool isEntrance = not std::ranges::binary_search(strandContactVertices, v1, assemblyGraph.orderById);

                if(isEntrance) {
                    SegmentPair& segmentPair = segmentPairs.emplace_back();

                    SegmentInfo& segmentInfo0 = segmentPair.segmentInfos[0];
                    SegmentInfo& segmentInfo1 = segmentPair.segmentInfos[1];

                    segmentInfo0.segment = segment;
                    segmentInfo0.id = id(segment);
                    segmentInfo0.isEntrance = true;
                    segmentInfo1.segment = segmentRc;
                    segmentInfo1.id = id(segmentRc);
                    segmentInfo1.isExit = true;

                    segmentPair.length = assemblyGraph[segment].length();
                    SHASTA2_ASSERT(segmentPair.length == assemblyGraph[segmentRc].length());

                    segmentPair.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
                    SHASTA2_ASSERT(segmentPair.coverage == assemblyGraph[segmentRc].lengthWeightedAverageCoverage());
                }
            }
        }
    }

    std::ranges::sort(segmentPairs, {}, &SegmentPair::id0);
    writeSegmentPairs();
}



void StrandContact::writeSegmentPairs()
{
    if(not html) {
        return;
    }

    html << std::fixed << std::setprecision(1);
    html << "<h2>Segment pairs</h2>"
        "<table>"
        "<tr>"
        "<th>Segment<br>pair id"
        "<th>Segment 0"
        "<th>Segment 1"
        "<th>Length"
        "<th>Coverage"
        "<th>Segment 0<br>is entrance"
        "<th>Segment 0<br>is exit"
        "<th>Segment 1<br>is entrance"
        "<th>Segment 1<br>is exit";

    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        const SegmentPair& segmentPair = segmentPairs[segmentPairId];
        html <<
            "<tr>"
            "<td class=centered>" << segmentPairId <<
            "<td class=centered>" << segmentPair.segmentInfos[0].id <<
            "<td class=centered>" << segmentPair.segmentInfos[1].id <<
            "<td class=centered>" << segmentPair.length <<
            "<td class=centered>" << segmentPair.coverage;

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[0].isEntrance) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[0].isExit) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[1].isEntrance) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[1].isExit) {
            html << "&check;";
        }
    }

    html << "</table>";



    // Write a csv file that cna be loaded in Bandage to show this StrandContact
    // with its entrances ane exits.
    const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-Bandage.csv";
    ofstream csv(fileName);
    csv << "Segment,Classification,Color\n";
    for(const SegmentPair& segmentPair: segmentPairs) {
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            csv << segmentInfo.id << ",";
            if(segmentInfo.isEntrance) {
                SHASTA2_ASSERT(not segmentInfo.isExit);
                csv << "Entrance,";
                csv << hslToRgbString(0.333, 0.5, 0.6) << ",";  // Green
            } else if(segmentInfo.isExit) {
                SHASTA2_ASSERT(not segmentInfo.isEntrance);
                csv << "Exit,";
                csv << hslToRgbString(0., 0.5, .6) << ",";      // Red
            } else {
                csv << "Internal,";
                csv << hslToRgbString(0.6, 0.5, .6) << ",";     // Blue
            }
            csv << endl;
        }
    }

}



void StrandContact::countReadOccurrences()
{
    // Count occurrences of reads in the first Segment of each SegmentPair.
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        const SegmentPair& segmentPair = segmentPairs[segmentPairId];

        // Use the first Segment of the SegmentPair.
        const Segment segment = segmentPair.segmentInfos[0].segment;

        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();
                readOccurrenceMap[readId].emplace_back(ReadOccurrence(segmentPairId, strand));
            }
        }
    }

    // Deduplicate and count the occurrences for each ReadId.
    vector<uint64_t> count;
    for(auto&[readId, occurrences]: readOccurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }

    // Remove from the map reads that occur in just one Segment.
    for(auto it=readOccurrenceMap.begin(); it!=readOccurrenceMap.end(); ++it) {
        auto itNext = it;
        ++itNext;
        if(it->second.size() == 1) {
            readOccurrenceMap.erase(it);
        }
        it = itNext;
    }

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
            const Segment segment = segmentPairs[occurrence.segmentPairId].segmentInfos[0].segment;
            csv << readId << ",";
            csv << occurrence.strand << ",";
            csv << id(segment) << ",";
            csv << occurrence.frequency << "\n";
        }
    }

    html << "<h2>Read occurrences</h2>"
        "For details of occurrences of reads in the segments, see "
        "<a href='"<< fileName << "'>" << fileName << "</a>.";
}



void StrandContact::createBipartiteGraph()
{
    // Create the vertices corresponding to Segments.
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        SegmentPair& segmentPair = segmentPairs[segmentPairId];
        for(uint64_t segmentIndexInPair=0; segmentIndexInPair<2; segmentIndexInPair++) {
            SegmentInfo& segmentInfo = segmentPair.segmentInfos[segmentIndexInPair];
            segmentInfo.v = bipartiteGraph.addVertex(segmentPairId, segmentIndexInPair);
        }
    }

    // Create the vertices corresponding to OrientedReadIds.
    for(const auto&[readId, ignore]: readOccurrenceMap) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const BipartiteGraph::vertex_descriptor v = bipartiteGraph.addVertex(orientedReadId);
            orientedReadIdVertexMap.insert({orientedReadId, v});
        }
    }



    // Now generate the edges.
    // Each ReadOccurrence generates a pair of reverse complemented edges.
    for(const auto&[readId, occurrences]: readOccurrenceMap) {


        // Now loop over the occurrences of this ReadId.
        for(const auto& occurrence: occurrences) {

            // Get the two OrientedReadIds of this ReadId.
            const OrientedReadId orientedReadId0(readId, occurrence.strand);;
            const OrientedReadId orientedReadId1(readId, 1 - occurrence.strand);;

            // Get corresponding BipartiteGraph vertices
            const BipartiteGraph::vertex_descriptor vOrientedRead0 =
                orientedReadIdVertexMap.at(orientedReadId0);
            const BipartiteGraph::vertex_descriptor vOrientedRead1 =
                orientedReadIdVertexMap.at(orientedReadId1);

            // The SegmentPair contains the two BipartiteGraph vertices
            // for this SegmentPair.
            const SegmentPair& segmentPair = segmentPairs[occurrence.segmentPairId];

            const BipartiteGraph::vertex_descriptor vSegment0 = segmentPair.segmentInfos[0].v;
            auto[e, ignore] = boost::add_edge(vOrientedRead0, vSegment0,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);

            const BipartiteGraph::vertex_descriptor vSegment1 = segmentPair.segmentInfos[1].v;
            auto[eRc, ignoreRc] = boost::add_edge(vOrientedRead1, vSegment1,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);

            bipartiteGraph.edgePairs.push_back({e, eRc, occurrence.frequency});

        }
    }

    // Sort the EdgePairs by decreasing frequency.
    std::ranges::sort(
        bipartiteGraph.edgePairs,
        std::greater<uint64_t>(),
        &BipartiteGraph::EdgePair::frequency);
}



// Add a vertex representing a Segment.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(uint64_t segmentPairId, uint64_t segmentIndexInPair)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(segmentPairId, segmentIndexInPair), bipartiteGraph);
}




void StrandContact::writeBipartiteGraph()
{
    if(not html) {
        return;
    }

    const string dotFileName = debugOutputBaseName + "-StrandContact-" +
        to_string(strandContactId) + ".dot";
    bipartiteGraph.writeGraphviz(dotFileName, segmentPairs, assemblyGraph);

    const double timeout = 30.;
    const string options = "-Nshape=point -Epenwidth=0.2 -Gratio=expand -Gsize=15";
    html << "<h2>Bipartite graph</h2>"
        "<br>In the bipartite graph, each vertex represents a segment or "
        "an oriented read. Oriented reads are displayed as small dots. "
        "<br><br>" << dotFileName << "<br>";

    try {
        graphvizToHtml(dotFileName, "sfdp", timeout, options, html, true);
    } catch (std::exception&) {
        html << "The bipartite graph is too complex to display.";
    }
}



void BipartiteGraph::writeGraphviz(
    const string& fileName,
    const vector<SegmentPair>& segmentPairs,
    const AssemblyGraph& assemblyGraph) const
{
    const BipartiteGraph& bipartiteGraph = *this;

    ofstream dot(fileName);

    dot << "graph BipartiteGraph {\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex = bipartiteGraph[v];
        const string color = randomHslColor(vertex.componentId, 0.5, 0.6);

        if(vertex.isSegment) {
            const uint64_t segmentPairId = vertex.segmentPairId;
            const uint64_t segmentIndexInPair = vertex.segmentIndexInPair;
            const SegmentPair& segmentPair = segmentPairs[segmentPairId];
            const SegmentInfo& segmentInfo = segmentPair.segmentInfos[segmentIndexInPair];
            const Segment segment = segmentInfo.segment;
            dot << assemblyGraph.id(segment);
            dot << " [width=0.1";
        } else {
            dot << "\"" << vertex.orientedReadId << "\"";
            dot << " [width=0.02";
        }
        dot << " color=\"" << color << "\"";
        dot << "]";
        dot << ";\n";
    }



    // Edges.
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const vertex_descriptor v0 = source(e, bipartiteGraph);
        const vertex_descriptor v1 = target(e, bipartiteGraph);
        const BipartiteGraphVertex& vertex0 = bipartiteGraph[v0];
        const BipartiteGraphVertex& vertex1 = bipartiteGraph[v1];

        const uint64_t frequency = bipartiteGraph[e].frequency;
        const double thickness = 0.1 * (1. + std::log10(frequency));

        if(vertex0.isSegment) {
            const uint64_t segmentPairId0 = vertex0.segmentPairId;
            const uint64_t segmentIndexInPair0 = vertex0.segmentIndexInPair;
            const SegmentPair& segmentPair0 = segmentPairs[segmentPairId0];
            const SegmentInfo& segmentInfo0 = segmentPair0.segmentInfos[segmentIndexInPair0];
            const Segment segment0 = segmentInfo0.segment;
            dot << assemblyGraph.id(segment0);
        } else {
            dot << "\"" << vertex0.orientedReadId << "\"";
        }

        dot << "--";

        if(vertex1.isSegment) {
            const uint64_t segmentPairId1 = vertex1.segmentPairId;
            const uint64_t segmentIndexInPair1 = vertex1.segmentIndexInPair;
            const SegmentPair& segmentPair1 = segmentPairs[segmentPairId1];
            const SegmentInfo& segmentInfo1 = segmentPair1.segmentInfos[segmentIndexInPair1];
            const Segment segment1 = segmentInfo1.segment;
            dot << assemblyGraph.id(segment1);
        } else {
            dot << "\"" << vertex1.orientedReadId << "\"";
        }

        dot << "[";
        dot << "penwidth=\"" << thickness << "\"";

        if(vertex0.componentId != vertex1.componentId) {
            dot << " color=red";
        }

        dot << "]";

        dot << ";\n";

    }

    dot << "}\n";
}

// Add a vertex representing an OrientedReadId.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(OrientedReadId orientedReadId)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(orientedReadId), bipartiteGraph);

}



void BipartiteGraph::strandSeparation()
{
    vector<uint64_t> edgePairsIndexes(edgePairs.size());
    std::ranges::iota(edgePairsIndexes, 0);
    strandSeparation(edgePairsIndexes);
}



// Return the reverse complement of a vertex.
// Because vertices are added in reverse complemented pairs,
// pairs of reverse complemented vertices have consecutive vertex_descriptors.
BipartiteGraph::vertex_descriptor BipartiteGraph::reverseComplement(vertex_descriptor v) const
{
    return v ^ 1;
}



// This processes the EdgePairs in the order described by
// the edgePairsIndexes.
// It stores the components in the components vector
// and also fills in the componentIndex in all the vertices.
// Reverse complemented components are numbered consecutively.
void BipartiteGraph::strandSeparation(const vector<uint64_t>& edgePairsIndexes)
{
    BipartiteGraph& bipartiteGraph = *this;

    DisjointSets disjointSets(num_vertices(bipartiteGraph));

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

        const bool strandViolationA = (disjointSets.findSet(v1A) == disjointSets.findSet(v0ARc));
        const bool strandViolationB = (disjointSets.findSet(v1B) == disjointSets.findSet(v0BRc));
        const bool strandViolationARc = (disjointSets.findSet(v0A) == disjointSets.findSet(v1ARc));
        const bool strandViolationBRc = (disjointSets.findSet(v0B) == disjointSets.findSet(v1BRc));

        const bool strandViolation = strandViolationA;
        SHASTA2_ASSERT(strandViolationB == strandViolation);
        SHASTA2_ASSERT(strandViolationARc == strandViolation);
        SHASTA2_ASSERT(strandViolationBRc == strandViolation);

        if(not strandViolation) {
            disjointSets.unionSet(v0A, v1A);
            disjointSets.unionSet(v0B, v1B);
        }
    }

    // Gather the components.
    disjointSets.gatherComponents(1, components);

    // Reorder the components so pairs of
    // reverse complemented components are numbered consecutively.
    // Because of the way vertices are created, this can be done simply by
    // sorting by the first index of each component.
    class SortByFirstElement {
    public:
    public:
         bool operator()(const vector<vertex_descriptor>& v0, const vector<vertex_descriptor>& v1) const
        {
             SHASTA2_ASSERT(not v0.empty());
             SHASTA2_ASSERT(not v1.empty());
             return v0.front() < v1.front();
        }
    };
    sort(components.begin(), components.end(), SortByFirstElement());

    // Store the componentId of the vertices.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<vertex_descriptor>& component = components[componentId];
        for(const vertex_descriptor v: component) {
            bipartiteGraph[v].componentId = componentId;
        }
    }
}



void StrandContact::evaluateStrandSeparation(bool writeCsvFile)
{

    ofstream csv;
    if(writeCsvFile) {
        csv.open(debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-EvaluateStrandSeparation.csv");
    }

    csv <<
        "Vertex,Segment,OrientedReadId,"
        "Total edge count,Cross-strand edge count,Cross-strand edge ratio,"
        "Total edge frequency,Cross-strand edge frequency,Cross-strand edge frequency ratio,\n";

    BGL_FORALL_VERTICES(v0, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex0 = bipartiteGraph[v0];
        uint64_t totalEdgeCount = 0;
        uint64_t crossStrandEdgeCount = 0;
        uint64_t totalEdgeFrequency = 0;
        uint64_t crossStrandEdgeFrequency = 0;
        BGL_FORALL_OUTEDGES(v0, e, bipartiteGraph, BipartiteGraph) {
            const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
            const BipartiteGraphVertex& vertex1 = bipartiteGraph[v1];
            const uint64_t frequency = bipartiteGraph[e].frequency;
            ++totalEdgeCount;
            totalEdgeFrequency += frequency;
            if(vertex0.componentId != vertex1.componentId) {
                ++crossStrandEdgeCount;
                crossStrandEdgeFrequency += frequency;
            }
        }
        const double crossStrandEdgeRatio = double(crossStrandEdgeCount) / double(totalEdgeCount);
        const double crossStrandEdgeFrequencyRatio = double(crossStrandEdgeFrequency) / double(totalEdgeFrequency);

        if(writeCsvFile) {
            csv << v0 << ",";

            if(vertex0.isSegment) {
                csv << segmentPairs[vertex0.segmentPairId].segmentInfos[vertex0.segmentIndexInPair].id;
            }
            csv << ",";

            if(not vertex0.isSegment) {
                csv << vertex0.orientedReadId;
            }
            csv << ",";

            csv << totalEdgeCount << ",";

            if(crossStrandEdgeCount) {
                csv << crossStrandEdgeCount;
            }
            csv << ",";

            if(crossStrandEdgeCount) {
                csv << crossStrandEdgeRatio;
            }
            csv << ",";

            csv << totalEdgeFrequency << ",";

            if(crossStrandEdgeFrequency) {
                csv << crossStrandEdgeFrequency;
            }
            csv << ",";

            if(crossStrandEdgeFrequency) {
                csv << crossStrandEdgeFrequencyRatio;
            }
            csv << ",";

            csv << "\n";
        }
    }

    html << "</table>";
}
