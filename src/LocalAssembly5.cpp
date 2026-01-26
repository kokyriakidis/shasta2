// Shasta.
#include "LocalAssembly5.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "graphvizToHtml.hpp"
#include "Kmer.hpp"
#include "MarkerKmers.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "tmpDirectory.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "algorithm.hpp"
#include <cmath>
#include "fstream.hpp"



LocalAssembly5::LocalAssembly5(
    const Anchors& anchors,
    [[maybe_unused]] uint64_t abpoaMaxLength,
    ostream& html,
    [[maybe_unused]] bool debug,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    anchors(anchors),
    html(html),
    leftAnchorId(anchorPair.anchorIdA),
    rightAnchorId(anchorPair.anchorIdB)
{
    SHASTA2_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));
    if(html) {
        writeInput(anchorPair, additionalOrientedReadIds);
    }

    gatherAllOrientedReads(anchorPair, additionalOrientedReadIds);
    if(html) {
        writeAllOrientedReadIds();
    }

    fillMarkerInfos();
    if(html) {
        writeMarkerInfos();
    }

    gatherOrientedReads();
    estimateOffset();
    gatherKmers();
    fillOrdinalsForAssembly();
    if(html) {
        writeOrientedReads();
    }

    createGraph();
    if(html) {
        writeGraph();
    }



    SHASTA2_ASSERT(0);
}



void LocalAssembly5::gatherAllOrientedReads(
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
{
    std::ranges::set_union(anchorPair.orientedReadIds, additionalOrientedReadIds,
        back_inserter(allOrientedReadIds));

}



void LocalAssembly5::fillMarkerInfos()
{
    // Get the Kmers.
    leftKmer = anchors.anchorKmer(leftAnchorId);
    rightKmer = anchors.anchorKmer(rightAnchorId);

    // Get the MarkerInfos.
    anchors.markerKmers.get(leftKmer, leftMarkerInfos);
    anchors.markerKmers.get(rightKmer, rightMarkerInfos);
}



// This gathers the oriented reads that will be used used in this local assembly.
// These are all OrientedReadIds that are in allOrientedReadIds and
// also on the left or right Marker.
// This can be done more efficiently but the code would be less readable.
void LocalAssembly5::gatherOrientedReads()
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    // Comparator for MarkerInfos that takes into account just the OrientedReadId.
    // This is safe because the left and right Markers are Anchors,
    // and so an OrientedReadId can appear more than once.
    class MarkerInfoComparator {
    public:
        bool operator() (const MarkerInfo& x, const MarkerInfo& y) const
        {
            return x.orientedReadId < y.orientedReadId;
        }
    };
    MarkerInfoComparator markerInfoComparator;



    // Loop over allOrientedReadIds.
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {

        // Look for left and right MarkerInfos containing this oriented read.
        const auto itLeft = std::lower_bound(leftMarkerInfos.begin(), leftMarkerInfos.end(),
            MarkerInfo(orientedReadId, 0), markerInfoComparator);
        const auto itRight = std::lower_bound(rightMarkerInfos.begin(), rightMarkerInfos.end(),
            MarkerInfo(orientedReadId, 0), markerInfoComparator);

        const bool isOnLeftMarker =
            (itLeft != leftMarkerInfos.end()) and
            (itLeft->orientedReadId == orientedReadId);
        const bool isOnRightMarker =
            (itRight != rightMarkerInfos.end()) and
            (itRight->orientedReadId == orientedReadId);

        // We can use it if it is on at least one of the left and right markers.
        // If on both, we also have to check that the right ordinal is greater than the left ordinal
        // (order violation).
        if(isOnLeftMarker or isOnRightMarker) {

            // Check for order violation.
            const bool isOnBothMarkers = isOnLeftMarker and isOnRightMarker;
            const bool orderViolation = isOnBothMarkers and (itRight->ordinal <= itLeft->ordinal);
            if(not orderViolation) {

                // Ok, we can use this oriented read.
                const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
                OrientedReadInfo& orientedReadInfo = orientedReadInfos.emplace_back();
                orientedReadInfo.orientedReadId = orientedReadId;
                orientedReadInfo.isOnLeftMarker = isOnLeftMarker;
                orientedReadInfo.isOnRightMarker = isOnRightMarker;
                if(isOnLeftMarker) {
                    orientedReadInfo.leftOrdinal = itLeft->ordinal;
                    orientedReadInfo.leftPosition = orientedReadMarkers[itLeft->ordinal].position + kHalf;
                }
                if(isOnRightMarker) {
                    orientedReadInfo.rightOrdinal = itRight->ordinal;
                    orientedReadInfo.rightPosition = orientedReadMarkers[itRight->ordinal].position + kHalf;
                }
            }
        }
    }
}



void LocalAssembly5::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        if(orientedReadInfo.isOnBothMarkers()) {
            sum += orientedReadInfo.positionOffset();
            ++n;
        }
    }

    SHASTA2_ASSERT(n > 0);
    offset = uint32_t(std::round(double(sum) / double(n)));

    if(html) {
        html << "<br>Estimated offset using " << n <<
            " oriented reads on the left and right markers is " << offset << " bases.";
    }
}


void LocalAssembly5::gatherKmers()
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);
    const uint32_t length = uint32_t(std::round((1. + drift) * double(offset)));

    vector<uint64_t> count;
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        orientedReadInfo.gatherKmers(
            kHalf,
            length,
            anchors.markers[orientedReadId.getValue()]);

        // Fill in allKmers.
        for(uint32_t ordinal=orientedReadInfo.firstOrdinalForAssembly;
            ordinal<=orientedReadInfo.lastOrdinalForAssembly; ordinal++) {
            orientedReadInfo.allKmers.push_back(anchors.markers.getKmer(orientedReadId, ordinal));
        }

        // Fill in nonUniqueKmers.
        orientedReadInfo.nonUniqueKmers = orientedReadInfo.allKmers;
        deduplicateAndCountWithThreshold(orientedReadInfo.nonUniqueKmers, count, 2UL);
    }
}



void LocalAssembly5::OrientedReadInfo::gatherKmers(
    uint32_t kHalf,
    uint32_t length,
    const span<const Marker>& orientedReadMarkers)
{
    if(isOnBothMarkers()) {
        firstOrdinalForAssembly = leftOrdinal;
        lastOrdinalForAssembly = rightOrdinal;
    }

    else if(isOnLeftMarker) {
        firstOrdinalForAssembly = leftOrdinal;

        if(firstOrdinalForAssembly ==  uint32_t(orientedReadMarkers.size() - 1)) {
            lastOrdinalForAssembly = uint32_t(orientedReadMarkers.size() - 1);
        } else {
            // Move right by at least length bases.
            for(lastOrdinalForAssembly=firstOrdinalForAssembly+1;
                lastOrdinalForAssembly<orientedReadMarkers.size()-1;
                lastOrdinalForAssembly++) {
                if(orientedReadMarkers[lastOrdinalForAssembly].position - leftPosition >= length) {
                    break;
                }
            }
        }

    }

    else if(isOnRightMarker) {
        lastOrdinalForAssembly = rightOrdinal;

        if(lastOrdinalForAssembly == 0) {
            firstOrdinalForAssembly = 0;
        } else {

            // Move left by at least length bases.
            for(firstOrdinalForAssembly=lastOrdinalForAssembly-1;
                /* Check later */ ; firstOrdinalForAssembly--) {
                if(rightPosition - orientedReadMarkers[firstOrdinalForAssembly].position >= length) {
                    break;
                }
                if(firstOrdinalForAssembly == 0) {
                    break;
                }
            }
        }

    }

    else {
        SHASTA2_ASSERT(0);
    }
    SHASTA2_ASSERT(firstOrdinalForAssembly < orientedReadMarkers.size());
    SHASTA2_ASSERT(lastOrdinalForAssembly < orientedReadMarkers.size());


    firstPositionForAssembly = orientedReadMarkers[firstOrdinalForAssembly].position + kHalf;
    lastPositionForAssembly = orientedReadMarkers[lastOrdinalForAssembly].position + kHalf;
}



void LocalAssembly5::writeInput(
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << anchorIdToString(leftAnchorId) <<
        "<tr><th class=left>Right anchor<td class=centered>" << anchorIdToString(rightAnchorId) <<
        "</table>";
    html << "<h3>OrientedReadIds in the AnchorPair</h3><table>";
    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << anchorPair.orientedReadIds.size() << " oriented reads in the AnchorPair.";

    html << "<h3>Additional OrientedReadIds</h3>"
        "These are additional OrientedReadIds that can be used in this assembly."
        "<table>";
    for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << additionalOrientedReadIds.size() << " additional oriented reads.";

}



void LocalAssembly5::writeAllOrientedReadIds() const
{
    html << "<h3>All OrientedReadIds</h3>"
        "These are the union of the OrientedReadIds in the AnchorPair "
        "and the additional OrientedReadIds"
        "<table>";
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << allOrientedReadIds.size() << " oriented reads.";
}



void LocalAssembly5::writeMarkerInfos() const
{
    writeMarkerInfos("Left", leftMarkerInfos);
    writeMarkerInfos("Right", rightMarkerInfos);
}



void LocalAssembly5::writeMarkerInfos(
    const string& side,
    const vector<MarkerInfo>& markerInfos) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    html << "<h3>" << side << " marker</h3>"
        "<table>"
        "<tr><th>OrientedReadId<th>Ordinal<th>Position";
    for(const MarkerInfo& markerInfo: markerInfos) {
        const auto orientedReadMarkers = anchors.markers[markerInfo.orientedReadId.getValue()];
        html <<
            "<tr>"
            "<td class=centered>" << markerInfo.orientedReadId <<
            "<td class=centered>" << markerInfo.ordinal <<
            "<td class=centered>" << orientedReadMarkers[markerInfo.ordinal].position + kHalf;
    }
    html << "</table>" << markerInfos.size() << " oriented reads.";

}



void LocalAssembly5::writeOrientedReads() const
{
    html <<
        "<h3>Oriented reads used in this local assembly</h3>"
        "<table><tr>"
        "<th>OrientedReadId"
        "<th>Marker<br>Count"
        "<th>Left<br>ordinal<th>Right<br>ordinal<th>Ordinal<br>offset"
        "<th>First<br>ordinal<br>for assembly<th>Last<br>ordinal<br>for assembly"
        "<th>Ordinal<br>offset<br>for assembly"
        "<th>Number<br>of markers<br>for assembly"
        "<th>Number<br>of non-unique<br>markers<br>for assembly"
        "<th>Initial<br>number<br>of markers<br>for assembly"
        "<th>Length"
        "<th>Left<br>position<th>Right<br>position<th>Position<br>offset"
        "<th>First<br>position<br>for assembly<th>Last<br>position<br>for assembly"
        "<th>Position<br>offset<br>for assembly";

    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const uint64_t markerCount = anchors.markers[orientedReadInfo.orientedReadId.getValue()].size();
        const uint64_t baseCount = anchors.reads.getRead(orientedReadInfo.orientedReadId.getReadId()).baseCount;


        html <<
            "<tr><td class=centered>" << orientedReadInfo.orientedReadId <<
            "<td class=centered>" << markerCount;



        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker) {
            html << orientedReadInfo.leftOrdinal;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightOrdinal;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker and orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightOrdinal - orientedReadInfo.leftOrdinal;
        }

        html << "<td class=centered>";
        html << orientedReadInfo.firstOrdinalForAssembly;
        html << "<td class=centered>";
        html << orientedReadInfo.lastOrdinalForAssembly;
        html << "<td class=centered>";
        html << orientedReadInfo.lastOrdinalForAssembly - orientedReadInfo.firstOrdinalForAssembly;
        html << "<td class=centered>";
        html << orientedReadInfo.allKmers.size();
        html << "<td class=centered>";
        html << orientedReadInfo.nonUniqueKmers.size();
        html << "<td class=centered>";
        html << orientedReadInfo.ordinalsForAssembly.size();

        html << "<td class=centered>" << baseCount;

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker) {
            html << orientedReadInfo.leftPosition;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightPosition;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker and orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightPosition - orientedReadInfo.leftPosition;
        }

        html << "<td class=centered>";
        html << orientedReadInfo.firstPositionForAssembly;
        html << "<td class=centered>";
        html << orientedReadInfo.lastPositionForAssembly;
        html << "<td class=centered>";
        html << orientedReadInfo.lastPositionForAssembly - orientedReadInfo.firstPositionForAssembly;
    }

    html << "</table>" << orientedReadInfos.size() << " oriented reads.";

}



void LocalAssembly5::fillOrdinalsForAssembly()
{
    // Gather all Kmers that are non unique in one or more reads.
    vector<Kmer> allNonUniqueKmers;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        std::ranges::copy(orientedReadInfo.nonUniqueKmers, back_inserter(allNonUniqueKmers));
    }
    deduplicate(allNonUniqueKmers);

    if(html) {
        html << "<br>Found " << allNonUniqueKmers.size() <<
            " marker k-mers that appear more than once in the assembly region of one or more oriented reads. "
            "These marker k-mers will not be used in graph generation.";
    }



    // Fill in the initial ordinalsForAssembly.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint32_t ordinal=orientedReadInfo.firstOrdinalForAssembly;
            ordinal<=orientedReadInfo.lastOrdinalForAssembly; ordinal++) {
            const Kmer& kmer = orientedReadInfo.allKmers[ordinal - orientedReadInfo.firstOrdinalForAssembly];
            if(not binary_search(allNonUniqueKmers.begin(), allNonUniqueKmers.end(), kmer)) {
                orientedReadInfo.ordinalsForAssembly.push_back(ordinal);
            }
        }
    }
}



// Create the graph using the current ordinalsForAssembly.
void LocalAssembly5::createGraph()
{
    LocalAssembly5& graph = *this;

    // Clear first, to permit iteration.
    graph.clear();

    // Gather all the Kmers.
    vector<Kmer> kmers;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint32_t ordinal: orientedReadInfo.ordinalsForAssembly) {
            kmers.push_back(orientedReadInfo.getKmer(ordinal));
        }
    }
    deduplicate(kmers);



    // Generate vertices.
    // Each Kmer generates a vertex.
    // The vertex_descriptor is the index in the kmers vector.
    vLeft = null_vertex();
    vRight = null_vertex();
    for(const Kmer& kmer: kmers) {
        if(kmer == leftKmer) {
            vLeft = num_vertices(graph);
        }
        if(kmer == rightKmer) {
            vRight = num_vertices(graph);
        }
        add_vertex(graph);
    }
    SHASTA2_ASSERT(vLeft != null_vertex());
    SHASTA2_ASSERT(vRight != null_vertex());



    // Fill in the infos of the vertices.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];
        for(uint32_t ordinal: orientedReadInfo.ordinalsForAssembly) {
            const Kmer& kmer = orientedReadInfo.getKmer(ordinal);
            const auto it = std::lower_bound(kmers.begin(), kmers.end(), kmer);
            SHASTA2_ASSERT(it != kmers.end());
            SHASTA2_ASSERT(*it == kmer);
            const vertex_descriptor v = it - kmers.begin();
            LocalAssembly5Vertex& vertex = graph[v];
            auto& info = vertex.infos.emplace_back();
            info.orientedReadIndex = orientedReadIndex;
            info.ordinal = ordinal;
        }
    }


    // Generate the edges.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];
        for(uint64_t i1=1; i1<orientedReadInfo.ordinalsForAssembly.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint32_t ordinal0 = orientedReadInfo.ordinalsForAssembly[i0];
            const uint32_t ordinal1 = orientedReadInfo.ordinalsForAssembly[i1];
            const Kmer& kmer0 = orientedReadInfo.getKmer(ordinal0);
            const Kmer& kmer1 = orientedReadInfo.getKmer(ordinal1);
            const auto it0 = std::lower_bound(kmers.begin(), kmers.end(), kmer0);
            SHASTA2_ASSERT(it0 != kmers.end());
            SHASTA2_ASSERT(*it0 == kmer0);
            const auto it1 = std::lower_bound(kmers.begin(), kmers.end(), kmer1);
            SHASTA2_ASSERT(it1 != kmers.end());
            SHASTA2_ASSERT(*it1 == kmer1);
            const vertex_descriptor v0 = it0 - kmers.begin();
            const vertex_descriptor v1 = it1 - kmers.begin();
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                SHASTA2_ASSERT(edgeExists);
            }
            LocalAssembly5Edge& edge = graph[e];
            auto& info = edge.infos.emplace_back();
            info.orientedReadIndex = orientedReadIndex;
            info.ordinal0 = ordinal0;
            info.ordinal1 = ordinal1;
        }
    }



    if(html) {
        html << "<br>The LocalAssembly5 graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
    }

}



void LocalAssembly5::writeGraphviz(const string& fileName)
{
    ofstream dot(fileName);
    writeGraphviz(dot);

}



void LocalAssembly5::writeGraphviz(ostream& dot)
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;

    // To make the output more readable, we output vertices in
    // approximate topological order.
    vector< pair<edge_descriptor, uint64_t> > edgesWithCoverage;
    BGL_FORALL_EDGES(e, graph, Graph) {
        edgesWithCoverage.push_back(make_pair(e, graph[e].infos.size()));
    }
    sort(edgesWithCoverage.begin(), edgesWithCoverage.end(),
        OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
    vector<edge_descriptor> edgesOrderedByDecreasingCoverage;
    for(const auto& [e, ignore]: edgesWithCoverage) {
        edgesOrderedByDecreasingCoverage.push_back(e);
    }
    approximateTopologicalSort(graph, edgesOrderedByDecreasingCoverage);
    vector< pair<vertex_descriptor, uint64_t> > verticesOrderedByRank;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        verticesOrderedByRank.push_back(make_pair(v, graph[v].rank));
    }
    sort(verticesOrderedByRank.begin(), verticesOrderedByRank.end(),
        OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());


    dot << "digraph LocalAssembly5 {\n";

    for(const auto& [v, ignore]: verticesOrderedByRank) {
        const LocalAssembly5Vertex& vertex = graph[v];
        dot << v;

        // Begin vertex attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << vertex.infos.size() << "\"";

        // Color.
        if(v== vLeft) {
            dot << " style=filled fillcolor=LightGreen";
        }
        if(v== vRight) {
            dot << " style=filled fillcolor=LightPink";
        }

        // End vertex attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        const LocalAssembly5Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << v0 << "->" << v1;

        // Begin edge attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << edge.infos.size() << "\"";

        // Thickness.
        dot << " penwidth=" << std::fixed << std::setprecision(2) << 0.5 * double(edge.infos.size());

        // End edge attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }

    dot << "}\n";
}



void LocalAssembly5::writeGraph()
{
    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName);


    // Display it in html in svg format.
    const double timeout = 120.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    try {
        graphvizToHtml(dotFileName, "dot", timeout, options, html);
    } catch(const std::exception& exception) {
        html << "<br>Error during graph layout: " << exception.what();
    }

}
