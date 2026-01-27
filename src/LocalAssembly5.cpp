// Shasta.
#include "LocalAssembly5.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Kmer.hpp"
#include "longestPath.hpp"
#include "MarkerKmers.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "tmpDirectory.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "algorithm.hpp"
#include <cmath>
#include "fstream.hpp"
#include <queue>



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
    // Sanity check.
    SHASTA2_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));

    // Summarize the input for this local assembly.
    if(html) {
        writeInput(anchorPair, additionalOrientedReadIds);
    }

    // Gather the oriented reads that could possibly participate in this local assembly.
    gatherAllOrientedReads(anchorPair, additionalOrientedReadIds);
    if(false) {
        writeAllOrientedReadIds();
    }

    // Gather the MarkerInfos on the left and right marker
    // for this local assembly.
    fillMarkerInfos();
    if(false) {
        writeMarkerInfos();
    }

    // Gather the oriented reads that will actually be used in this local assembly.
    gatherOrientedReads();

    // Estimate the base offset for this local assembly using the oriented reads
    // that appear on both the left and right markers.
    estimateOffset();

    // Fill in the LocalRegions of all oriented reads.
    // These are the oriented read regions that will be used in this local assembly.
    fillLocalRegions();

    // Gather the Kmers that will be used in this local assembly.
    gatherKmers();
    if(html) {
        writeKmers();
        writeOrientedReads();
    }

    // Create the initial graph.
    createGraph();
    removeInaccessibleVertices();
    computeDominatorTree();
    if(html) {
        html << "<h3>Initial LocalAssembly5 graph</h3>";
        writeGraph();
    }

    // Graph cleanup.
    removeLowCoverageEdges();
    removeInaccessibleVertices();
    if(html) {
        html << "<h3>LocalAssembly5 graph after cleanup</h3>";
        writeGraph();
    }

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



void LocalAssembly5::fillLocalRegions()
{
    const uint32_t length = uint32_t(std::round((1. + drift) * double(offset)));

    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        orientedReadInfo.fillLocalRegion(anchors, length);
    }
}



void LocalAssembly5::OrientedReadInfo::fillLocalRegion(
    const Anchors& anchors,
    uint32_t length)
{
    const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    if(isOnBothMarkers()) {
        localRegion.firstOrdinal = leftOrdinal;
        localRegion.lastOrdinal = rightOrdinal;
    }

    else if(isOnLeftMarker) {
        localRegion.firstOrdinal = leftOrdinal;

        if(localRegion.firstOrdinal ==  uint32_t(orientedReadMarkers.size() - 1)) {
            localRegion.lastOrdinal = uint32_t(orientedReadMarkers.size() - 1);
        } else {
            // Move right by at least length bases.
            for(localRegion.lastOrdinal = localRegion.firstOrdinal + 1;
                localRegion.lastOrdinal<orientedReadMarkers.size();
                localRegion.lastOrdinal++) {
                if(orientedReadMarkers[localRegion.lastOrdinal].position - leftPosition >= length) {
                    break;
                }
            }
        }

    }

    else if(isOnRightMarker) {
        localRegion.lastOrdinal = rightOrdinal;

        if(localRegion.lastOrdinal == 0) {
            localRegion.firstOrdinal = 0;
        } else {

            // Move left by at least length bases.
            for(localRegion.firstOrdinal = localRegion.lastOrdinal - 1;
                /* Check later */ ; localRegion.firstOrdinal--) {
                if(rightPosition - orientedReadMarkers[localRegion.firstOrdinal].position >= length) {
                    break;
                }
                if(localRegion.firstOrdinal == 0) {
                    break;
                }
            }
        }

    }

    else {
        SHASTA2_ASSERT(0);
    }
    SHASTA2_ASSERT(localRegion.firstOrdinal < orientedReadMarkers.size());
    SHASTA2_ASSERT(localRegion.lastOrdinal < orientedReadMarkers.size());


    localRegion.firstPosition = orientedReadMarkers[localRegion.firstOrdinal].position + kHalf;
    localRegion.lastPosition = orientedReadMarkers[localRegion.lastOrdinal].position + kHalf;


    // Fill in the LocalRegion Kmers.
    for(uint32_t ordinal = localRegion.firstOrdinal;
        ordinal <= localRegion.lastOrdinal; ordinal++) {
        localRegion.kmers.push_back(anchors.markers.getKmer(orientedReadId, ordinal));
    }

    // Fill in nonUniqueKmers.
    localRegion.computeNonUniqueKmers();
}



void LocalAssembly5::OrientedReadInfo::LocalRegion::computeNonUniqueKmers()
{
    nonUniqueKmers = kmers;
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(nonUniqueKmers, count, 2UL);
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
        "<th>Left<br>marker<br>ordinal"
        "<th>Right<br>marker<br>ordinal"
        "<th>Left<br>to right<br>ordinal<br>offset"
        "<th>Local<br>region<br>first<br>ordinal"
        "<th>Local<br>region<br>last<br>ordinal"
        "<th>Local<br>region<br>ordinal<br>offset"
        "<th>Local<br>region<br>number<br>of markers"
        "<th>Number<br>of markers<br>for assembly"
        "<th>Read<br>length"
        "<th>Left<br>marker<br>position"
        "<th>Right<br>marker<br>position"
        "<th>Left<br>to right<br>position<br>offset"
        "<th>Local<br>region<br>first<br>position"
        "<th>Local<br>region<br>last<br>position"
        "<th>Local<br>region<br>position<br>offset<br>(length)";

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
        html << orientedReadInfo.localRegion.firstOrdinal;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.lastOrdinal;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.ordinalOffset();
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.kmers.size();
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
        html << orientedReadInfo.localRegion.firstPosition;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.lastPosition;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.positionOffset();
    }

    html << "</table>" << orientedReadInfos.size() << " oriented reads.";

}



void LocalAssembly5::gatherKmers()
{
    // Gather all Kmers that are not unique
    // in the LocalRegion of one or more reads.
    vector<Kmer> allNonUniqueKmers;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        std::ranges::copy(orientedReadInfo.localRegion.nonUniqueKmers, back_inserter(allNonUniqueKmers));
    }
    deduplicate(allNonUniqueKmers);

    if(html) {
        html << "<br>Found " << allNonUniqueKmers.size() <<
            " marker k-mers that appear more than once in the local rethis local assembly graph generation.";
    }


    // Gather the Kmers that are unique
    // in the LocalRegions of all oriented reads.
    // Fill in the ordinalsForAssembly in the OrientedReadInfos.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint32_t ordinal = orientedReadInfo.localRegion.firstOrdinal;
            ordinal <= orientedReadInfo.localRegion.lastOrdinal; ordinal++) {

            const Kmer& kmer = orientedReadInfo.localRegion.getKmer(ordinal);
            if(not binary_search(allNonUniqueKmers.begin(), allNonUniqueKmers.end(), kmer)) {
                orientedReadInfo.ordinalsForAssembly.push_back(ordinal);
                kmers.push_back(kmer);
            }
        }
    }
    deduplicate(kmers);
}



// Create the graph using the current ordinalsForAssembly.
void LocalAssembly5::createGraph()
{
    LocalAssembly5& graph = *this;

    // Clear first, just in case.
    graph.clear();

    // Generate vertices.
    // Each Kmer generates a vertex.
    vector<vertex_descriptor> vertexTable;
    vLeft = null_vertex();
    vRight = null_vertex();
    for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
        const Kmer& kmer = kmers[kmerId];
        const vertex_descriptor v = add_vertex(LocalAssembly5Vertex(kmerId), graph);
        vertexTable.push_back(v);
        if(kmer == leftKmer) {
            vLeft = v;
        }
        if(kmer == rightKmer) {
            vRight = v;
        }
    }
    SHASTA2_ASSERT(vLeft != null_vertex());
    SHASTA2_ASSERT(vRight != null_vertex());



    // Fill in the infos of the vertices.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];
        for(uint32_t ordinal: orientedReadInfo.ordinalsForAssembly) {
            const Kmer& kmer = orientedReadInfo.localRegion.getKmer(ordinal);
            const auto it = std::lower_bound(kmers.begin(), kmers.end(), kmer);
            SHASTA2_ASSERT(it != kmers.end());
            SHASTA2_ASSERT(*it == kmer);
            const vertex_descriptor v = vertexTable[it - kmers.begin()];
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
            const Kmer& kmer0 = orientedReadInfo.localRegion.getKmer(ordinal0);
            const Kmer& kmer1 = orientedReadInfo.localRegion.getKmer(ordinal1);
            const auto it0 = std::lower_bound(kmers.begin(), kmers.end(), kmer0);
            SHASTA2_ASSERT(it0 != kmers.end());
            SHASTA2_ASSERT(*it0 == kmer0);
            const auto it1 = std::lower_bound(kmers.begin(), kmers.end(), kmer1);
            SHASTA2_ASSERT(it1 != kmers.end());
            SHASTA2_ASSERT(*it1 == kmer1);
            const vertex_descriptor v0 = vertexTable[it0 - kmers.begin()];
            const vertex_descriptor v1 = vertexTable[it1 - kmers.begin()];
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
        dot << vertex.kmerId;

        // Begin vertex attributes.
        dot << "[";

        // Label.
        dot << "label=\"V" << vertex.kmerId << "\\n" << vertex.infos.size() << "\"";

        // Tooltip.
        dot << " tooltip=\"";
        kmers[vertex.kmerId].write(dot, anchors.k);
        dot << "\"";

        // Color.
        if(v== vLeft) {
            dot << " style=filled fillcolor=LightGreen";
        }
        else if(v== vRight) {
            dot << " style=filled fillcolor=LightPink";
        }
        else if(vertex.isOnDominatorTreePath) {
            dot << " style=filled fillcolor=LightBlue";
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
        dot << graph[v0].kmerId << "->" << graph[v1].kmerId;

        // Begin edge attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << edge.infos.size() << "\"";

        // Tooltip.
        dot << " tooltip=\"" << edge.infos.size() << "\"";

        // Thickness.
        dot << " penwidth=" << std::fixed << std::setprecision(2) << 0.3 * double(edge.infos.size());

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



// Remove vertices that are not forward accessible from vLeft
// and backward accessible from vRight.
void LocalAssembly5::removeInaccessibleVertices()
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;

    std::set<vertex_descriptor> leftAccessibleVertices;
    findReachableVertices(graph, vLeft, 0, leftAccessibleVertices);

    std::set<vertex_descriptor> rightAccessibleVertices;
    findReachableVertices(graph, vRight, 1, rightAccessibleVertices);

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const bool keep =  leftAccessibleVertices.contains(v) and rightAccessibleVertices.contains(v);
        if(not keep) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

}



// Compute a dominator tree starting at vLeft
// and the dominator tree path from vLeft to vRight.
void LocalAssembly5::computeDominatorTree()
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;

    // Compute the dominator tree.
    // This sets the dominator field of the vertices
    // reachable by leftAnchorVertex.
    lengauer_tarjan_dominator_tree_general(graph, vLeft);

    // To compute the dominator tree path from vLeft
    // to vRight we walk back the dominator tree
    // starting at vRight.
    dominatorTreePath.push_back(vRight);
    while(true) {
        const vertex_descriptor v = dominatorTreePath.back();
        const vertex_descriptor dominator = graph[v].dominator;
        if(dominator == null_vertex()) {
            break;
        } else {
            dominatorTreePath.push_back(dominator);
        }
    }
    std::ranges::reverse(dominatorTreePath);

    for(const vertex_descriptor v: dominatorTreePath) {
        graph[v].isOnDominatorTreePath = true;
    }
    SHASTA2_ASSERT(dominatorTreePath.front() == vLeft);
    SHASTA2_ASSERT(dominatorTreePath.back() == vRight);
}



// Find the longest path.
void LocalAssembly5::computeAssemblyPath()
{
    vector<edge_descriptor> longestPath;
    shasta2::longestPath(*this, longestPath);
}



// The assembly path is computed by joining together
// partial assembly paths between adjacent vertices
// of the dominator tree path.
void LocalAssembly5::removeLowCoverageEdges()
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;

    // We use color for BFSs.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        graph[v].color = 0;
    }

    // Remove low coverage edges, using a different threshold
    // for each leg of the dominator tree path.
    for(uint64_t i1=1; i1<dominatorTreePath.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const vertex_descriptor v0 = dominatorTreePath[i0];
        const vertex_descriptor v1 = dominatorTreePath[i1];
        removeLowCoverageEdges(v0, v1);
    }
}



// Given two vertices which are adjacent in the dominator tree path, vA and vB,
// remove low coverage edges in-between without destroying reachability
// of vB from vA.
void LocalAssembly5::removeLowCoverageEdges(
    vertex_descriptor vA,
    vertex_descriptor vB)
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;
    const bool debug = false;

    if(debug) {
        cout << "Working on assembly path portion between V" <<
            graph[vA].kmerId << " and V" << graph[vB].kmerId << endl;
    }

    // Find all the vertices "between" vA and vB.
    // Use a BFS that starts at vA and stops at vB.
    // This sets the color of these vertices to 1.
    vector<vertex_descriptor> inBetweenVertices;
    std::queue<vertex_descriptor> q;
    q.push(vA);
    graph[vA].color = 1;
    graph[vB].color = 1;
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);
            if(v1 != vB) {
                if(graph[v1].color == 0) {
                    inBetweenVertices.push_back(v1);
                    graph[v1].color = 1;
                    q.push(v1);
                }
            }
        }

    }

    if(debug) {
        cout << "This portion of the graph contains " << inBetweenVertices.size() << " vertices:" << endl;
        for(const vertex_descriptor v: inBetweenVertices) {
            cout << "V" << graph[v].kmerId << " ";
        }
        cout << endl;
    }



    // A vertex predicate that select only vertices with color>0.
    class VertexPredicate {
    public:
        bool operator()(const vertex_descriptor& v) const
        {
            return (*graph)[v].color > 0;
        }
        VertexPredicate(const Graph& graph) : graph(&graph) {}
        VertexPredicate() : graph(0) {}
        const Graph* graph;
    };
    const VertexPredicate vertexPredicate(graph);



    // An edge predicate that select only edges with coverage >= minCoverage.
    uint64_t minCoverage = 0;
    class EdgePredicate {
    public:
        bool operator()(const edge_descriptor& e) const
        {
            return (*graph)[e].infos.size() >= *minCoveragePointer;
        }
        EdgePredicate(const Graph& graph, uint64_t* minCoveragePointer) :
            graph(&graph),
            minCoveragePointer(minCoveragePointer) {}
        EdgePredicate() : graph(0) {}
        const Graph* graph;
        uint64_t* minCoveragePointer;
    };
    EdgePredicate edgePredicate(graph, &minCoverage);

    // A filtered graph that includes all vertices with color>0 and
    // all the edges with coverage >= minCoverage.
    using FilteredGraph = boost::filtered_graph<Graph, EdgePredicate, VertexPredicate>;
    FilteredGraph filteredGraph(graph, edgePredicate, vertexPredicate);



    // Find the highest edge coverage threshold that ensures reachability.
    for(minCoverage=1; ; ++minCoverage) {
        if(not isReachable(filteredGraph, vA, vB, 0)) {
            --minCoverage;
            break;
        }
    }
    if(debug) {
        cout << "Edge coverage threshold for reachability is " << minCoverage << endl;
    }

    // Remove edges with lower coverage.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_OUTEDGES(vA, e, graph, Graph) {
        if(graph[e].infos.size() < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const vertex_descriptor v: inBetweenVertices) {
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            if(graph[e].infos.size() < minCoverage) {
                edgesToBeRemoved.push_back(e);
            }
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }


    // Reset the colors.
    graph[vA].color = 0;
    graph[vB].color = 0;
    for(const vertex_descriptor v: inBetweenVertices) {
        graph[v].color = 0;
    }

}



void LocalAssembly5::removeIsolatedVertices()
{
    using Graph = LocalAssembly5;
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) == 0) and (out_degree(v, graph)==0)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, graph);
    }
}



void LocalAssembly5::writeKmers() const
{
    html <<
        "<h3>Kmers</h3>"
        "<table><tr><th>Vertex<th>Kmer";

    for(uint64_t vertexId=0; vertexId<kmers.size(); vertexId++) {
        html <<
            "<tr><td class=centered>" << vertexId <<
            "<td style='font-family:monospace'>";
        kmers[vertexId].write(html, anchors.k);
    }

    html << "</table>";
}
