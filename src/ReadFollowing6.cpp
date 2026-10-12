// Shasta2.
#include "ReadFollowing6.hpp"
#include "AssemblyGraph.hpp"
#include "computeLayout.hpp"
#include "deduplicate.hpp"
#include "color.hpp"
using namespace shasta2;
using namespace ReadFollowing6;

// Standard library.
#include <iomanip>
#include "fstream.hpp"



ReadFollowing6::Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
    ostream& html,
    uint64_t tangleId,
    const string& debugOutputBaseName) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVertices),
    html(html),
    tangleId(tangleId),
    debugOutputBaseName(debugOutputBaseName)
{
    if(html) {
        html << "<h1>ReadFollowing6</h1>";
    }

    checkTangleVertices();

    gatherSegments();
    writeSegments();

    countReadOccurrences();
    writeReadOccurrences();

    createBipartiteGraph();
    writeBipartiteGraph();

}


// The tangleVertices must be sorted by id.
// If the Tangle includes a vertex, it cannot also include its reverse complement.
void ReadFollowing6::Tangle::checkTangleVertices()
{
    SHASTA2_ASSERT(std::ranges::is_sorted(tangleVertices, assemblyGraph.orderById));

    for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
        const AssemblyGraph::vertex_descriptor vRc = assemblyGraph[v].vRc;
        SHASTA2_ASSERT(not std::ranges::binary_search(tangleVertices, vRc, assemblyGraph.orderById));
    }
}



void ReadFollowing6::Tangle::gatherSegments()
{

    // Out-edges of the tangleVertices give us Segments
    // internal to the Tangle plus the exits.
    for(const AssemblyGraph::vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, segment, assemblyGraph, AssemblyGraph) {

            // If the target of this Segment is not a Tangle vertex, this is an exit.
            const AssemblyGraph::vertex_descriptor v1 = target(segment, assemblyGraph);
            const bool isExit = not std::ranges::binary_search(tangleVertices, v1, assemblyGraph.orderById);

            SegmentInfo& segmentInfo = segmentInfos.emplace_back();

            segmentInfo.segment = segment;
            segmentInfo.id = id(segment);
            segmentInfo.isExit = isExit;
            segmentInfo.length = assemblyGraph[segment].length();
            segmentInfo.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
        }
    }



    // In-edges of the tangleVertices give us the entrances.
    for(const AssemblyGraph::vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_INEDGES(v0, segment, assemblyGraph, AssemblyGraph) {

            // If the source of this Segment is not a Tangle vertex, this is an entrance.
            const AssemblyGraph::vertex_descriptor v1 = source(segment, assemblyGraph);
            const bool isEntrance = not std::ranges::binary_search(tangleVertices, v1, assemblyGraph.orderById);

            if(isEntrance) {
                SegmentInfo& segmentInfo = segmentInfos.emplace_back();

                segmentInfo.segment = segment;
                segmentInfo.id = id(segment);
                segmentInfo.isEntrance = isEntrance;
                segmentInfo.length = assemblyGraph[segment].length();
                segmentInfo.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
            }
        }
    }

    std::ranges::sort(segmentInfos, {}, &SegmentInfo::id);

}



void ReadFollowing6::Tangle::writeSegments()
{
    if(not html) {
        return;
    }

    html << std::fixed << std::setprecision(1);
    html << "<h2>Segments</h2>"
        "<table>"
        "<tr>"
        "<th>Segment<br>index"
        "<th>Segment"
        "<th>Length"
        "<th>Coverage"
        "<th>Is<br>entrance"
        "<th>Is>br>exit";

    for(uint64_t segmentIndex=0; segmentIndex<segmentInfos.size(); segmentIndex++) {
        const SegmentInfo& segmentInfo = segmentInfos[segmentIndex];
        html <<
            "<tr>"
            "<td class=centered>" << segmentIndex <<
            "<td class=centered>" << segmentInfo.id <<
            "<td class=centered>" << segmentInfo.length <<
            "<td class=centered>" << segmentInfo.coverage;

        html << "<td class=centered>";
        if(segmentInfo.isEntrance) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentInfo.isExit) {
            html << "&check;";
        }
    }
    html << "</table>";



    // Write a csv file that can be loaded in Bandage to show this Tangle
    // with its entrances and exits.
    const string fileName = debugOutputBaseName + "-Tangle-" + to_string(tangleId) + "-Bandage.csv";
    ofstream csv(fileName);
    csv << "Segment,Classification,Color\n";
    for(const SegmentInfo& segmentInfo: segmentInfos) {
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



uint64_t ReadFollowing6::Tangle::id(AssemblyGraphBaseClass::vertex_descriptor v) const
{
    return assemblyGraph.id(v);
}



uint64_t ReadFollowing6::Tangle::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void ReadFollowing6::Tangle::countReadOccurrences()
{
    // The readOccurrenceMap is a filtered version of the completeReadOccurrenceMap.
    std::map<OrientedReadId, vector<ReadOccurrence> > completeReadOccurrenceMap;

    // Count occurrences of oriented reads in the Segments.
    for(uint64_t segmentIndex=0; segmentIndex<segmentInfos.size(); segmentIndex++) {
        const SegmentInfo& segmentInfo = segmentInfos[segmentIndex];
        const Segment segment = segmentInfo.segment;

        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                completeReadOccurrenceMap[orientedReadId].emplace_back(ReadOccurrence(segmentIndex));
            }
        }
    }



    // Deduplicate and count the occurrences for each OrientedReadId.
    vector<uint64_t> count;
    for(auto&[orientedReadId, occurrences]: completeReadOccurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }



    // The readOccurrenceMap is a filtered version of the completeReadOccurrenceMap.
    for(auto&[orientedReadId, completeOccurrences]: completeReadOccurrenceMap) {

        // Find the entrances and exits in these occurrences.
        vector<uint64_t> entranceSegmentIndexes;
        vector<uint64_t> exitSegmentIndexes;
        for(const auto& occurrence: completeOccurrences) {
            const uint64_t segmentIndex = occurrence.segmentIndex;
            const SegmentInfo& segmentInfo = segmentInfos[segmentIndex];
            if(segmentInfo.isEntrance) {
                entranceSegmentIndexes.push_back(segmentIndex);
            }
            if(segmentInfo.isExit) {
                exitSegmentIndexes.push_back(segmentIndex);
            }
        }
        deduplicate(entranceSegmentIndexes);
        deduplicate(exitSegmentIndexes);


        // If there is more than one entrance, we remove occurrences
        // on all entrances.
        // If there is more than one exit, we remove occurrences
        // on all entrances.
        vector<ReadOccurrence> occurrences;
        if((entranceSegmentIndexes.size() < 2) and (exitSegmentIndexes.size() < 2)) {
            occurrences = completeOccurrences;
        } else {
            for(const auto& occurrence: completeOccurrences) {
                const uint64_t segmentIndex = occurrence.segmentIndex;
                if((entranceSegmentIndexes.size() >= 2) and (std::ranges::binary_search(entranceSegmentIndexes, segmentIndex))) {
                    continue;
                }
                if((exitSegmentIndexes.size() >= 2) and (std::ranges::binary_search(exitSegmentIndexes, segmentIndex))) {
                    continue;
                }
            }
        }

        // If there is only one occurrence, skip this OrientedReadId.
        if(occurrences.size() < 2) {
            continue;
        }

        readOccurrenceMap.insert({orientedReadId, occurrences});
    }


}



void ReadFollowing6::Tangle::writeReadOccurrences()
{
    if(not html) {
        return;
    }

    const string fileName = debugOutputBaseName + "-Tangle-" + to_string(tangleId) + "-ReadOccurrences.csv";
    ofstream csv(fileName);
    csv << "OrientedReadId,Segment,Frequency\n";
    for(const auto&[orientedReadId, occurrences]: readOccurrenceMap) {
        for(const auto& occurrence: occurrences) {
            const Segment segment = segmentInfos[occurrence.segmentIndex].segment;
            csv << orientedReadId << ",";
            csv << id(segment) << ",";
            csv << occurrence.frequency << "\n";
        }
    }

    html << "<h2>Read occurrences</h2>"
        "For details of occurrences of oriented reads in the segments, see "
        "<a href='"<< fileName << "'>" << fileName << "</a>.";

}



void ReadFollowing6::Tangle::createBipartiteGraph()
{
    // Create the vertices corresponding to Segments.
    for(uint64_t segmentIndex=0; segmentIndex<segmentInfos.size(); segmentIndex++) {
        SegmentInfo& segmentInfo = segmentInfos[segmentIndex];
        segmentInfo.v = bipartiteGraph.addVertex(segmentIndex);
    }

    // Create the vertices corresponding to OrientedReadIds.
    for(const auto&[orientedReadId, _]: readOccurrenceMap) {
        const BipartiteGraph::vertex_descriptor v = bipartiteGraph.addVertex(orientedReadId);
        orientedReadIdVertexMap.insert({orientedReadId, v});
    }



    // Now generate the edges.
    for(const auto&[orientedReadId, occurrences]: readOccurrenceMap) {

        // Compute the total frequency for this OrientedReadId.
        uint64_t totalFrequency = 0;
        for(const auto& occurrence: occurrences) {
            totalFrequency += occurrence.frequency;
        }

        // Get corresponding BipartiteGraph vertex.
        const BipartiteGraph::vertex_descriptor vOrientedRead =
            orientedReadIdVertexMap.at(orientedReadId);

        // Now loop over the occurrences of this OrientedReadId.
        for(const auto& occurrence: occurrences) {
            SHASTA2_ASSERT(occurrences.size() > 1);

            // Get the BipartiteGraph vertex for this Segment.
            const SegmentInfo& segmentInfo = segmentInfos[occurrence.segmentIndex];
            const BipartiteGraph::vertex_descriptor vSegment = segmentInfo.v;

            const uint64_t frequency = occurrence.frequency;
            const double weight = double(frequency) / double(totalFrequency);

            boost::add_edge(vOrientedRead, vSegment,
                BipartiteGraphEdge({frequency, weight}), bipartiteGraph);
        }
    }


}


// Add a vertex representing a Segment.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(uint64_t segmentIndex)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(segmentIndex), bipartiteGraph);
}



// Add a vertex representing an OrientedReadId.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(OrientedReadId orientedReadId)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(orientedReadId), bipartiteGraph);
}




void ReadFollowing6::Tangle::writeBipartiteGraph()
{
    if(not html) {
        return;
    }

    html << "<h3>Bipartite graph</h3>"
        "The bipartite graph has " << num_vertices(bipartiteGraph) <<
        " vertices and " << num_edges(bipartiteGraph)<< " edges.";

    // If command "customLayout" is not available, don't do anything.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA2_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    if(returnCode != 0) {
        html << "<br>The bipartite graph cannot be displayes because command "
            "customLayout is not available.";
        return;
    }

    // Create a map containing the desired length for each edge.
    std::map<BipartiteGraph::edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        // const BipartiteGraph::vertex_descriptor v0 = source(e, bipartiteGraph);
        // const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
        // const uint64_t frequency = bipartiteGraph[e].frequency;
        // const double length = std::pow(frequency, -0.1);
        const double length = 1.1 - bipartiteGraph[e].weight;
        edgeLengthMap.insert({e, length});
    }

    // Compute the layout.
    std::map<BipartiteGraph::vertex_descriptor, array<double, 2> > positionMap;
    const int quality = 2;
    const double timeout = 30.;
    const auto layoutReturnCode = computeLayoutCustom(bipartiteGraph, edgeLengthMap, positionMap, quality, timeout);
    if(layoutReturnCode != ComputeLayoutReturnCode::Success) {
        html << "<br>The custom layout of the bipartite graph cannot be displayed.";
        return;
    }

    // Compute the bounding box of the layout.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    for(const auto& p: positionMap) {
        const array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
    }

    // Enlarge the bounding box a bit.
    const double extend = 0.05 * max(xMax-xMin, yMax-yMin);
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;

    // Make it square,
    if((xMax - xMin) > (yMax - yMin)) {
        const double delta = ((xMax - xMin) - (yMax - yMin)) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = ((yMax - yMin) - (xMax - xMin)) / 2.;;
        xMin -= delta;
        xMax += delta;
    }


    // Begin the svg.
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    const uint64_t sizePixels = 900;
    html << std::scientific;
    const string svgId = "BipartiteGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << xMin << " " << yMin << " " <<
        xMax-xMin << " " <<
        yMax-yMin << "'"
        " style='background-color:#f0f0f0'"
        ">\n";



    // Write the edges first so they don't obscure the vertices.
    const double edgeThicknessFactor = (xMax - xMin) * 0.0003;
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraph::vertex_descriptor v0 = source(e, bipartiteGraph);
        const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
        const auto&[x0, y0] = positionMap.at(v0);
        const auto&[x1, y1] = positionMap.at(v1);
        const string color = "Black";

        // const uint64_t frequency = bipartiteGraph[e].frequency;
        // const double thickness = edgeThicknessFactor * std::pow(frequency, 0.5);
        const double thickness = edgeThicknessFactor * (1. + 10. * bipartiteGraph[e].weight);

        html <<
            "<line x1='" << x0 << "' y1='" << y0 <<
            "' x2='" << x1 << "' y2='" << y1 <<
            "' stroke='" << color <<
            "' stroke-width='" << thickness <<
            "' />";
    }



    // Write the vertices.
    const double segmentRadius = (xMax - xMin) * 0.01;
    const double orientedReadRadius = (xMax - xMin) * 0.004;
    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex = bipartiteGraph[v];
        const string color = (vertex.isSegment ? "Green" : "Red");
        const double radius = (bipartiteGraph[v].isSegment ? segmentRadius : orientedReadRadius);
        const auto&[x, y] = positionMap.at(v);

        html << "<g><title>";
        if(vertex.isSegment) {
            html << id(segmentInfos[vertex.segmentIndex].segment);
        } else {
            html << vertex.orientedReadId;
        }
        html << "</title>";

        html << "<circle cx='" << x << "' cy='" << y <<
            "' fill='" << color <<
            "' r='" << radius <<
            "' /></g>";
    }


    // Finish the svg.
    html << "</svg></div>";
}
