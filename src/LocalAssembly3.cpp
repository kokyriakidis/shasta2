// Shasta.
#include "LocalAssembly3.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "graphvizToHtml.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>



// This assembles between anchorIdA and anchorIdB
// of the given AnchorPair. It uses all the OrientedReadIds
// stored in the AnchorPair, and which all appear in both
// anchorIdA and anchorIdB. In addition, it uses OrientedReadIds
// stored in additionalOrientedReadIds that:
// - Are not also in the AnchorPair.
// - Appear in at least one of anchorIdA and anchorIdB.
LocalAssembly3::LocalAssembly3(
    const Anchors& anchors,
    ostream& html,
    bool debug,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    leftAnchorId(anchorPair.anchorIdA),
    rightAnchorId(anchorPair.anchorIdB)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double drift = 0.1;

    SHASTA_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));
    if(html) {
        writeInput(html, debug, anchorPair, additionalOrientedReadIds);
    }

    gatherOrientedReads(anchors,anchorPair, additionalOrientedReadIds);
    estimateOffset();
    if(html) {
        html << "<p>Estimated offset is " << offset << " bases.";
    }

    fillFirstLastOrdinalForAssembly(anchors.markers, drift);
    if(html) {
        writeOrientedReads(anchors, html);
    }

    fillOrientedReadKmers(anchors.markers);
    gatherKmers(anchors);
    if(html) {
        html << "<p>Found " << kmers.size() << " distinct kmers." << endl;
        if(debug) {
            writeOrientedReadKmers(html);
        }
    }

    createVertices();
    createEdges();
    computeDominatorTree();
    if(html) {
        html << "<p>The graph for this local assembly has " << num_vertices(*this) <<
            " vertices and " << num_edges(*this) << " edges.";
        writeHtml(html, anchors.markers);
    }
}



void LocalAssembly3::writeInput(
    ostream& html,
    bool debug,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << anchorIdToString(leftAnchorId) <<
        "<tr><th class=left>Right anchor<td class=centered>" << anchorIdToString(rightAnchorId) <<
        "</table>";

    if(debug) {
        html << "<h3>OrientedReadIds in the AnchorPair</h3><table>";
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";

        html << "<h3>Additional OrientedReadIds</h3><table>";
        for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }

}



// Gather all the OrientedReadIds that we can consider for use in this assembly.
void LocalAssembly3::gatherOrientedReads(
    const Anchors& anchors,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
    {

    // Merge together the OrientedReadIds in the AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> inputOrientedReadIds = anchorPair.orientedReadIds;
    std::ranges::copy(additionalOrientedReadIds, back_inserter(inputOrientedReadIds));
    deduplicate(inputOrientedReadIds);

    // Access the two Anchors.
    const Anchor leftAnchor = anchors[leftAnchorId];
    const Anchor rightAnchor = anchors[rightAnchorId];

    // To gather the OrientedReadIds that we can use for assembly,
    // use a joint loop over the inputOrientedReadIds and the oriented reads in
    // the two anchors.
    auto itLeft = leftAnchor.begin();
    const auto endLeft = leftAnchor.end();
    auto itRight = rightAnchor.begin();
    const auto endRight = rightAnchor.end();
    for(auto it=inputOrientedReadIds.begin(); it!=inputOrientedReadIds.end(); ++it) {
        const OrientedReadId orientedReadId = *it;

        while((itLeft != endLeft) and (itLeft->orientedReadId < orientedReadId)) {
            ++itLeft;
        }
        while((itRight != endRight) and (itRight->orientedReadId < orientedReadId)) {
            ++itRight;
        }
        const bool isOnLeftAnchor = (itLeft != endLeft) and (itLeft->orientedReadId == orientedReadId);
        const bool isOnRightAnchor = (itRight != endRight) and (itRight->orientedReadId == orientedReadId);

        if(isOnLeftAnchor or isOnRightAnchor) {
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            OrientedReadInfo orientedReadInfo;
            orientedReadInfo.orientedReadId = orientedReadId;
            orientedReadInfo.isOnLeftAnchor = isOnLeftAnchor;
            orientedReadInfo.isOnRightAnchor = isOnRightAnchor;
            if(isOnLeftAnchor) {
                orientedReadInfo.leftOrdinal = itLeft->ordinal;
                orientedReadInfo.leftPosition = orientedReadMarkers[itLeft->ordinal].position;
            }
            if(isOnRightAnchor) {
                orientedReadInfo.rightOrdinal = itRight->ordinal;
                orientedReadInfo.rightPosition = orientedReadMarkers[itRight->ordinal].position;
            }
            orientedReadInfos.push_back(orientedReadInfo);
        }
    }
}



void LocalAssembly3::estimateOffset()
{
    uint64_t n = 0;
    uint64_t offsetSum = 0;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        if(orientedReadInfo.isOnBothAnchors()) {
            ++n;
            offsetSum += orientedReadInfo.positionOffset();
        }
    }
    SHASTA_ASSERT(n > 0);

    const double preciseOffset = double(offsetSum) / double(n);
    offset = uint32_t(std::round(preciseOffset));
}



// Use the estimated offset to fill in the firstOrdinalForAssembly and lastOrdinalForAssembly
// of each oriented read. These define the portion of this OrientedReadId
// sequence that will be used in this local assembly.
void LocalAssembly3::fillFirstLastOrdinalForAssembly(const Markers& markers, double drift)
{
    const uint32_t length = uint32_t(std::round((1. + drift) * double(offset)));

    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        orientedReadInfo.fillFirstLastOrdinalForAssembly(markers, length);
    }
}



void LocalAssembly3::OrientedReadInfo::fillFirstLastOrdinalForAssembly(
    const Markers& markers,
    uint32_t length)
{
    if(isOnLeftAnchor and isOnRightAnchor) {
        firstOrdinalForAssembly = leftOrdinal;
        lastOrdinalForAssembly = rightOrdinal;
    }

    else if(isOnLeftAnchor and not isOnRightAnchor){
        firstOrdinalForAssembly = leftOrdinal;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Move right by at least length bases.
        for(lastOrdinalForAssembly=firstOrdinalForAssembly;
            lastOrdinalForAssembly<orientedReadMarkers.size()-1; lastOrdinalForAssembly++) {
            if(orientedReadMarkers[lastOrdinalForAssembly].position - leftPosition >= length) {
                break;
            }
        }
    }

    else if(isOnRightAnchor and not isOnLeftAnchor){
        lastOrdinalForAssembly = rightOrdinal;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Move left by at least length bases.
        for(firstOrdinalForAssembly=lastOrdinalForAssembly; /* Check later */ ; firstOrdinalForAssembly--) {
            if(rightPosition - orientedReadMarkers[firstOrdinalForAssembly].position >= length) {
                break;
            }
            if(firstOrdinalForAssembly == 0) {
                break;
            }
        }
    }

    else {
        SHASTA_ASSERT(0);
    }
}



uint32_t LocalAssembly3::OrientedReadInfo::firstPositionForAssembly(const Markers& markers) const
{
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    return orientedReadMarkers[firstOrdinalForAssembly].position;
}



uint32_t LocalAssembly3::OrientedReadInfo::lastPositionForAssembly(const Markers& markers) const
{
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    return orientedReadMarkers[lastOrdinalForAssembly].position;
}



void LocalAssembly3::fillOrientedReadKmers(const Markers& markers)
{
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        orientedReadInfo.fillOrientedReadKmers(markers);
    }
}



void LocalAssembly3::OrientedReadInfo::fillOrientedReadKmers(const Markers& markers)
{
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    orientedReadKmerInfos.resize(lastOrdinalForAssembly - firstOrdinalForAssembly + 1);
    for(uint32_t ordinal=firstOrdinalForAssembly; ordinal<=lastOrdinalForAssembly; ordinal++) {
        orientedReadKmerInfos[ordinal - firstOrdinalForAssembly].kmer =
            markers.getKmer(orientedReadId, ordinal);
    }

}



void LocalAssembly3::gatherKmers(const Anchors& anchors)
{
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(const auto& kmerInfo: orientedReadInfo.orientedReadKmerInfos) {
            kmers.push_back(kmerInfo.kmer);
        }
    }
    deduplicate(kmers);

    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(auto& kmerInfo: orientedReadInfo.orientedReadKmerInfos) {
            kmerInfo.kmerIndex = getKmerIndex(kmerInfo.kmer);
        }
    }

    leftAnchorKmerIndex = getKmerIndex(anchors.anchorKmer(leftAnchorId));
    rightAnchorKmerIndex = getKmerIndex(anchors.anchorKmer(rightAnchorId));
}



// Return the index of the given Kmer in the kmers vector.
uint64_t LocalAssembly3::getKmerIndex(const Kmer& kmer) const
{
    const auto it = std::lower_bound(kmers.begin(), kmers.end(), kmer);
    SHASTA_ASSERT(it != kmers.end());
    SHASTA_ASSERT(*it == kmer);
    return it - kmers.begin();
}



void LocalAssembly3::writeOrientedReads(
    const Anchors& anchors,
    ostream& html) const
{

    html << "<h3>Oriented reads used in this local assembly</h3>"
        "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>Number<br>of<br>markers"
        "<th>Sequence<br>read<br>length"
        "<th>On left<br>anchor"
        "<th>On right<br>anchor"
        "<th>On both<br>anchors"
        "<th>Left<br>ordinal"
        "<th>Right<br>ordinal"
        "<th>Ordinal<br>offset"
        "<th>Left<br>position"
        "<th>Right<br>position"
        "<th>Position<br>offset"
        "<th>First<br>ordinal<br>for assembly"
        "<th>Last<br>ordinal<br>for assembly"
        "<th>Ordinal<br>offset<br>for assembly"
        "<th>First<br>position<br>for assembly"
        "<th>Last<br>position<br>for assembly"
        "<th>Position<br>offset<br>for assembly";

    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        html <<
            "<tr><td class=centered>" << orientedReadInfo.orientedReadId <<
            "<td class=centered>" <<
            anchors.markers.size(orientedReadInfo.orientedReadId.getValue()) <<
            "<td class=centered>" <<
            anchors.reads.getReadSequenceLength(orientedReadInfo.orientedReadId.getReadId());

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << orientedReadInfo.leftOrdinal;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << orientedReadInfo.rightOrdinal;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << orientedReadInfo.ordinalOffset();
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << orientedReadInfo.leftPosition;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << orientedReadInfo.rightPosition;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << orientedReadInfo.positionOffset();
        }

        html <<
            "<td class=centered>" << orientedReadInfo.firstOrdinalForAssembly <<
            "<td class=centered>" << orientedReadInfo.lastOrdinalForAssembly <<
            "<td class=centered>" <<
            orientedReadInfo.lastOrdinalForAssembly -
            orientedReadInfo.firstOrdinalForAssembly <<
            "<td class=centered>" << orientedReadInfo.firstPositionForAssembly(anchors.markers) <<
            "<td class=centered>" << orientedReadInfo.lastPositionForAssembly(anchors.markers) <<
        "<td class=centered>" <<
            orientedReadInfo.lastPositionForAssembly(anchors.markers) -
            orientedReadInfo.firstPositionForAssembly(anchors.markers);
    }
    html << "</table>";

}



void LocalAssembly3::writeOrientedReadKmers(ostream& html) const
{
    html <<
        "<h3>Sequences of k-mer indexes</h3>"
        "<p>Left anchor has k-mer index " << leftAnchorKmerIndex <<
        "<p>Right anchor has k-mer index " << rightAnchorKmerIndex <<
        "<table>";
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        html << "<tr><th>" << orientedReadInfo.orientedReadId;
        for(const auto& kmerInfo: orientedReadInfo.orientedReadKmerInfos) {
            html << "<td class=centered";
            if(kmerInfo.kmerIndex == leftAnchorKmerIndex) {
                html << " style='background-color:LightBlue'";
            }
            if(kmerInfo.kmerIndex == rightAnchorKmerIndex) {
                html << " style='background-color:LightPink'";
            }
            html << ">" << kmerInfo.kmerIndex;
        }
    }
    html << "</table>";
}



void LocalAssembly3::createVertices()
{
    using Graph = LocalAssembly3;
    Graph& graph = *this;

    vertexMap.clear();
    vertexMap.resize(kmers.size(), null_vertex());

    // Loop over oriented reads.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];

        // Loop over Kmers of this oriented read.
        for(uint32_t ordinal=orientedReadInfo.firstOrdinalForAssembly;
            ordinal<=orientedReadInfo.lastOrdinalForAssembly; ordinal++) {
            const auto& orientedReadKmerInfo =
                orientedReadInfo.orientedReadKmerInfos[ordinal - orientedReadInfo.firstOrdinalForAssembly];
            const uint64_t kmerIndex = orientedReadKmerInfo.kmerIndex;

            // Get the vertex, creating it if necessary.
            vertex_descriptor v = vertexMap[kmerIndex];
            if(v == null_vertex()) {
                v = boost::add_vertex(LocalAssembly3Vertex(kmerIndex), graph);
                vertexMap[kmerIndex] = v;
            }

            // Store this marker in the vertex.
            graph[v].data.emplace_back(LocalAssembly3Vertex::Data({orientedReadIndex, ordinal}));
        }
    }

    leftAnchorVertex = vertexMap[leftAnchorKmerIndex];
    rightAnchorVertex = vertexMap[rightAnchorKmerIndex];
}



void LocalAssembly3::createEdges()
{
    using Graph = LocalAssembly3;
    Graph& graph = *this;

    // Loop over oriented reads.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];

        // Loop over Kmers of this oriented read.
        for(uint32_t ordinal0=orientedReadInfo.firstOrdinalForAssembly;
            ordinal0<orientedReadInfo.lastOrdinalForAssembly; ordinal0++) {

            const auto& orientedReadKmerInfo0 =
                orientedReadInfo.orientedReadKmerInfos[ordinal0 - orientedReadInfo.firstOrdinalForAssembly];
            const uint64_t kmerIndex0 = orientedReadKmerInfo0.kmerIndex;

            const uint32_t ordinal1 = ordinal0 + 1;
            const auto& orientedReadKmerInfo1 =
                orientedReadInfo.orientedReadKmerInfos[ordinal1 - orientedReadInfo.firstOrdinalForAssembly];
            const uint64_t kmerIndex1 = orientedReadKmerInfo1.kmerIndex;

            vertex_descriptor v0 = vertexMap[kmerIndex0];
            vertex_descriptor v1 = vertexMap[kmerIndex1];

            // Get the edge, creating it if necessary.
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                SHASTA_ASSERT(edgeExists);
            }

            const LocalAssembly3Edge::Data data({orientedReadIndex, ordinal0, ordinal1});
            graph[e].data.emplace_back(LocalAssembly3Edge::Data({orientedReadIndex, ordinal0, ordinal1}));
        }
    }

}



void LocalAssembly3::writeGraphviz(const string& fileName, const Markers& markers) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, markers);
}



void LocalAssembly3::writeGraphviz(ostream& dot, const Markers& markers) const
{
    using Graph = LocalAssembly3;
    const Graph& graph = *this;

    dot << "digraph LocalAssembly3 {\n";



    BGL_FORALL_VERTICES(v, graph, Graph) {
        const LocalAssembly3Vertex& vertex = graph[v];
        dot << vertex.kmerIndex;

        // Begin vertex attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << vertex.kmerIndex << "\\n" <<
            vertex.coverage() << "\"";

        // Color.
        if(v == leftAnchorVertex) {
            dot << " style=filled fillcolor=LightGreen";
        } else if(v == rightAnchorVertex) {
            dot << " style=filled fillcolor=LightPink";
        } else if(vertex.isOnDominatorTreePath) {
            dot << " style=filled fillcolor=LightBlue";
        }

        // End vertex attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        const LocalAssembly3Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        dot << graph[v0].kmerIndex << "->" << graph[v1].kmerIndex;

        // Begin edge attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << edge.coverage() << "\\n" <<
            edgePositionOffset(e, markers) << "\"";

        // Thickness.
        dot << " penwidth=" << std::setprecision(2) << 0.5 * double(edge.coverage());

        // End edge attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }
    dot << "}\n";
}



void LocalAssembly3::writeHtml(ostream& html, const Markers& markers) const
{
    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, markers);


    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);

}



// Compute the average position offset for an edge.
uint32_t LocalAssembly3::edgePositionOffset(
    edge_descriptor e,
    const Markers& markers) const
{
    using Graph = LocalAssembly3;
    const Graph& graph = *this;

    const LocalAssembly3Edge& edge = graph[e];

    uint64_t positionOffsetSum = 0;
    for(const auto& data: edge.data) {
        const OrientedReadId orientedReadId = orientedReadInfos[data.orientedReadIndex].orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t positionOffset =
            orientedReadMarkers[data.ordinal1].position -
            orientedReadMarkers[data.ordinal0].position;
        positionOffsetSum += positionOffset;
    }

    const double preciseOffset = double(positionOffsetSum) / double(edge.data.size());
    return uint32_t(std::round(preciseOffset));
}



// Compute a dominator tree starting at leftAnchorVertex
// and the dominator tree path from leftAnchorVertex
// to rightAnchorVertex.
void LocalAssembly3::computeDominatorTree()
{
    using Graph = LocalAssembly3;
    Graph& graph = *this;

    // Compute the dominator tree.
    // This sets the dominator field of the vertices
    // reachable by leftAnchorVertex.
    lengauer_tarjan_dominator_tree_general(graph, leftAnchorVertex);

    // To compute the dominator tree path from leftAnchorVertex
    // to rightAnchorVertex we walk back the dominator tree
    // starting at rightAnchorVertex.
    dominatorTreePath.push_back(rightAnchorVertex);
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
    SHASTA_ASSERT(dominatorTreePath.front() == leftAnchorVertex);
    SHASTA_ASSERT(dominatorTreePath.back() == rightAnchorVertex);
}
