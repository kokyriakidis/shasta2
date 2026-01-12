// Shasta.
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "Assembler.hpp"
#include "Options.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "LocalAnchorGraph.hpp"
#include "LocalReadGraph.hpp"
#include "LocalReadAnchorGraph.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "ReadId.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>



void Assembler::exploreAnchor(const vector<string>& request, ostream& html)
{
    const uint64_t k = assemblerInfo->k;

    // Get the request parameters.
    string anchorIdString;
    const bool anchorIdStringIsPresent = HttpServer::getParameterValue(request, "anchorIdString", anchorIdString);
    boost::trim(anchorIdString);

    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    // Begin the form.
    html <<
        "<h2>Anchor information</h2>"
        "<form>"
        "<table>"

    // AnchorId
        "<tr><th class=left>Anchor id<td class=centered>"
        "<input type=text name=anchorIdString required style='text-align:center'";
    if(anchorIdStringIsPresent) {
        html << " value='" << anchorIdString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>";

    // Assembly stage for annotations.
    html <<
        "<tr title='Leave blank for no annotations'>"
        "<th class=left>Assembly stage for annotations"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center'";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Show anchor details'> "
        "</form>";



    // If the anchor id missing or invalid, stop here.
    if(not anchorIdStringIsPresent) {
        return;
    }
    const AnchorId anchorId = anchorIdFromString(anchorIdString);

    if((anchorId == invalid<AnchorId>) or (anchorId >= anchors().size())) {
        html << "<p>Invalid anchor id. Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }


    html << "<h1>Anchor " << anchorIdString << "</h1>";

    const auto markerInfos = anchors()[anchorId];
    const uint64_t coverage = markerInfos.size();
    const vector<Base> kmerSequence = anchors().anchorKmerSequence(anchorId);


    vector<AnchorId> parents;
    vector<uint64_t> parentsCoverage;
    anchors().findParents(journeys(), anchorId, parents, parentsCoverage);

    vector<AnchorId> children;
    vector<uint64_t> childrenCoverage;
    anchors().findChildren(journeys(), anchorId, children, childrenCoverage);

    // Write a summary table.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "<tr><th class=left>K-mer sequence<td class=centered style='font-family:monospace'>";
    copy(kmerSequence.begin(), kmerSequence.end(), ostream_iterator<Base>(html));

    html <<
        "<tr><th class=left>Parent anchors"
        "<td class=centered>";
    for(uint64_t i=0; i<parents.size(); i++) {
        const string parentAnchorIdString = anchorIdToString(parents[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(parentAnchorIdString) << "'>" <<
            parentAnchorIdString << "</a> coverage " << parentsCoverage[i];

    }

    html <<
        "<tr><th class=left>Children anchors"
        "<td class=centered>";
    for(uint64_t i=0; i<children.size(); i++) {
        const string childAnchorIdString = anchorIdToString(children[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(childAnchorIdString) << "'>" <<
            childAnchorIdString << "</a> coverage " << childrenCoverage[i];

    }

    html <<
        "<tr><th class=left>Forward read following"
        "<td class=centered>";
    {
        AnchorPair anchorPair;
        uint32_t offset;
        const bool found = anchors().readFollowing(
            journeys(), anchorId, 0,
            httpServerData.options->minAnchorGraphEdgeCoverage,
            httpServerData.options->aDrift,
            httpServerData.options->bDrift,
            anchorPair, offset);
        if(found) {
            const string s = anchorIdToString(anchorPair.anchorIdB);
            html <<
                "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(s) << "'>" << s << "</a>"
                ", coverage " << anchorPair.size() <<
                ", offset " << offset;
        }
    }

    html <<
        "<tr><th class=left>Backward read following"
        "<td class=centered>";
    {
        AnchorPair anchorPair;
        uint32_t offset;
        const bool found = anchors().readFollowing(
            journeys(), anchorId, 1,
            httpServerData.options->minAnchorGraphEdgeCoverage,
            httpServerData.options->aDrift,
            httpServerData.options->bDrift,
            anchorPair, offset);
        if(found) {
            const string s = anchorIdToString(anchorPair.anchorIdA);
            html <<
                "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(s) << "'>" << s << "</a>"
                ", coverage " << anchorPair.size() <<
                ", offset " << offset;
        }
    }

    html <<
        "<tr><th class=left>Forward read following<br>(multiple)"
        "<td class=centered>";
    {
        vector< pair<AnchorPair, uint32_t> > anchorPairs;
        anchors().readFollowing(
            journeys(), anchorId, 0,
            httpServerData.options->minAnchorGraphEdgeCoverage,
            httpServerData.options->minAnchorGraphContinueReadFollowingCount,
            httpServerData.options->aDrift,
            httpServerData.options->bDrift,
            anchorPairs);
        for(uint64_t i=0; i<anchorPairs.size(); i++) {
            const auto& p = anchorPairs[i];
            const AnchorPair& anchorPair = p.first;
            const uint32_t offset = p.second;
            const string s = anchorIdToString(anchorPair.anchorIdB);
            html <<
                "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(s) << "'>" << s << "</a>"
                ", coverage " << anchorPair.size() <<
                ", offset " << offset;
            if(i != anchorPairs.size() - 1) {
                html << "<br>";
            }
        }
    }

    html <<
        "<tr><th class=left>Backward read following<br>(multiple)"
        "<td class=centered>";
    {
        vector< pair<AnchorPair, uint32_t> > anchorPairs;
        anchors().readFollowing(
            journeys(), anchorId, 1,
            httpServerData.options->minAnchorGraphEdgeCoverage,
            httpServerData.options->minAnchorGraphContinueReadFollowingCount,
            httpServerData.options->aDrift,
            httpServerData.options->bDrift,
            anchorPairs);
        for(uint64_t i=0; i<anchorPairs.size(); i++) {
            const auto& p = anchorPairs[i];
            const AnchorPair& anchorPair = p.first;
            const uint32_t offset = p.second;
            const string s = anchorIdToString(anchorPair.anchorIdA);
            html <<
                "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(s) << "'>" << s << "</a>"
                ", coverage " << anchorPair.size() <<
                ", offset " << offset;
            if(i != anchorPairs.size() - 1) {
                html << "<br>";
            }
        }
    }

    html << "</table>";



    // Assembly graph annotations, if requested.
    if(not assemblyStage.empty()) {

        // Get the AssemblyGraph for this assembly stage.
        const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
            assemblyStage,
            *httpServerData.options);
        const auto annotations = assemblyGraph.getAnnotations(anchorId);

        html << "<h2>Assembly graph annotations at assembly stage " << assemblyStage << "</h2>";

        if(annotations.empty()) {
            html << "This Anchor is not referenced in assembly stage " << assemblyStage;
        } else {
            html << "<ul>";

            for(const auto& annotation: annotations) {
                html << "<li>";

                if(annotation.v == AssemblyGraph::null_vertex()) {
                    // This AnchorId is used in a step.
                    const AssemblyGraphEdge& edge = assemblyGraph[annotation.e];
                    const string segmentUrl = "exploreSegment?assemblyStage=" + assemblyStage +
                        "&segmentName=" +to_string(edge.id);
                    const string stepUrl = "exploreSegmentStep?assemblyStage=" + assemblyStage +
                        "&segmentName=" + to_string(edge.id) + "&stepId=" + to_string(annotation.step);
                    html <<
                        "Segment <a href='" << segmentUrl << "'>" << edge.id << "</a>"
                        ", step <a href='" << stepUrl << "'>" << annotation.step << "</a>"
                        " of " << edge.size() <<
                        ", " <<
                        (annotation.isAnchorIdA ? "first" : "second") <<
                        " anchor.";

                } else {

                    // This AnchorId is used in a vertex.
                    const AssemblyGraph::vertex_descriptor v = annotation.v;
                    html << "Assembly graph vertex with";

                    // Incoming segments.
                    if(in_degree(v, assemblyGraph) == 0) {
                        html << " no incoming segments";
                    } else {
                        html << " incoming segments";
                        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
                            const AssemblyGraphEdge& edge = assemblyGraph[e];
                            const string segmentUrl = "exploreSegment?assemblyStage=" + assemblyStage +
                                "&segmentName=" + to_string(edge.id);
                            html << " <a href='" << segmentUrl << "'>" << edge.id << "</a>";
                        }
                        html << ",";
                    }

                    // Outgoing segments.
                    if(out_degree(v, assemblyGraph) == 0) {
                        html << " no outgoing segments";
                    } else {
                        html << " outgoing segments";
                        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
                            const AssemblyGraphEdge& edge = assemblyGraph[e];
                            const string segmentUrl = "exploreSegment?assemblyStage=" + assemblyStage +
                                "&segmentName=" + to_string(edge.id);
                            html << " <a href='" << segmentUrl << "'>" << edge.id << "</a>";
                        }
                        html << ".";
                    }
                }
            }

            html << "</ul>";
        }
    }



    std::map<pair<AnchorId, AnchorId>, uint64_t> tangleMatrix;

    // Write the marker intervals of this Anchor.
    html <<
        "<h2>Marker intervals</h2>"
        "<table>"
        "<tr>"
        "<th>Index"
        "<th>Oriented<br>read<br>id"
        "<th>Position<br>in<br>journey"
        "<th>Ordinal"
        "<th>Position"
        "<th>Previous<br>anchor<br>in journey"
        "<th>Next<br>anchor<br>in journey";

    // Loop over the marker intervals.
    for(uint64_t i=0; i<coverage; i++) {
        const AnchorMarkerInfo& markerInfo = markerInfos[i];
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys()[orientedReadId];

        const uint32_t ordinal = markerInfo.ordinal;

        const auto orientedReadMarkers = markers()[orientedReadId.getValue()];
        const uint32_t position = orientedReadMarkers[ordinal].position;

        AnchorId previousAnchorInJourney = invalid<AnchorId>;
        if(markerInfo.positionInJourney > 0) {
            previousAnchorInJourney = journey[markerInfo.positionInJourney - 1];
        }
        AnchorId nextAnchorInJourney = invalid<AnchorId>;
        if(markerInfo.positionInJourney < journey.size() - 1) {
            nextAnchorInJourney = journey[markerInfo.positionInJourney + 1];
        }

        html <<
            "<tr>"
            "<td class=centered>" << i;

        // The OrientedReadId is written with an hyperlink that will
        // display the portion of the oriented read around this Anchor.
        const string url =
            "exploreReadSequence?"
            "readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand()) +
            "&beginPosition=" + to_string((position > 2 * k) ? (position - 2 * k) : 0) +
            "&endPosition=" + to_string(position + 3 * k - 1);
        html <<
            "<td class=centered>" <<
            "<a href='" << url << "'>" <<
            orientedReadId << "</a>";

       html <<
            "<td class=centered>" << markerInfo.positionInJourney <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position;

       // Previous anchor in journey.
       html << "<td class=centered>";
       if(previousAnchorInJourney != invalid<AnchorId>) {
           html << anchorIdToString(previousAnchorInJourney);
       }

       // Next anchor in journey.
       html << "<td class=centered>";
       if(nextAnchorInJourney != invalid<AnchorId>) {
           html << anchorIdToString(nextAnchorInJourney);
       }

       auto it = tangleMatrix.find(make_pair(previousAnchorInJourney, nextAnchorInJourney));
       if(it == tangleMatrix.end()) {
           tangleMatrix.insert(make_pair(make_pair(previousAnchorInJourney, nextAnchorInJourney), 1));
       } else {
           ++(it->second);
       }
    }
    html << "</table>";



    html << "<h2>Tangle matrix</h2>";
    html << "<table><tr><th>In<th>Out<th>Coverage";
    for(const auto& p: tangleMatrix) {
        const AnchorId previousAnchorInJourney = p.first.first;
        const AnchorId nextAnchorInJourney = p.first.second;
        const uint64_t coverage = p.second;

        html << "<tr><td class=centered>";
        if(previousAnchorInJourney != invalid<AnchorId>) {
            html << anchorIdToString(previousAnchorInJourney);
        }

        html << "<td class=centered>";
        if(nextAnchorInJourney != invalid<AnchorId>) {
            html << anchorIdToString(nextAnchorInJourney);
        }

        html << "<td class=centered>" << coverage;
    }
}



void Assembler::exploreAnchorPair2(const vector<string>& request, ostream& html)
{
    // Get the parameters for the request
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    boost::trim(anchorIdAString);
    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);
    boost::trim(anchorIdBString);

    string adjacentInJourneyString;
    const bool adjacentInJourney = HttpServer::getParameterValue(request,
        "adjacentInJourney", adjacentInJourneyString);



    // Write the form.
    html << "<form><table>";

    html <<
        "<tr><th class=left>Anchor A"
        "<td class=centered><input type=text name=anchorIdAString required style='text-align:center'";
    if(anchorIdAStringIsPresent) {
        html << " value='" << anchorIdAString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr><th class=left>Anchor B"
        "<td class=centered><input type=text name=anchorIdBString required style='text-align:center'";
    if(anchorIdBStringIsPresent) {
        html << " value='" << anchorIdBString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>"

        "<tr><th>Adjacent in journey"
        "<td class=centered><input type=checkbox name=adjacentInJourney" <<
        (adjacentInJourney ? " checked" : "") <<
        ">"

        "</table>"
        "<input type=submit value='Get anchor pair information'>"
        "</form>";



    // Check the AnchorIds
    if(not (anchorIdAStringIsPresent and anchorIdBStringIsPresent)) {
        return;
    }
    const AnchorId anchorIdA = anchorIdFromString(anchorIdAString);
    const AnchorId anchorIdB = anchorIdFromString(anchorIdBString);

    if((anchorIdA == invalid<AnchorId>) or (anchorIdA >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdAString << ". Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if((anchorIdB == invalid<AnchorId>) or (anchorIdB >= anchors().size())) {
        html << "<p>Invalid anchor id " << anchorIdBString << " .Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }

    if(anchorIdA == anchorIdB) {
        html << "Specify two distinct anchors.";
        return;
    }



    // Create a AnchorPair from these two anchors.
    const AnchorPair anchorPair(anchors(), anchorIdA, anchorIdB, adjacentInJourney);

    // Output to html.
    html << "<h2>Anchor pair</h2>";
    anchorPair.writeAllHtml(html, anchors(), journeys());
}



void Assembler::exploreJourney(const vector<string>& request, ostream& html)
{

    html <<
        "<h2>Oriented read journey</h2>"
        "<p>The journey of an oriented read is the sequence of anchors it visits.";


    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = HttpServer::getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = HttpServer::getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = HttpServer::getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = HttpServer::getParameterValue(request, "endPosition", endPosition);

    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads().readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= reads().readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(reads().getReadSequenceLength(readId));
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }

    // Access the information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const span<const Marker> orientedReadMarkers = markers()[orientedReadId.getValue()];
    SHASTA2_ASSERT(journeys().isOpen());
    SHASTA2_ASSERT(journeys().size() == 2 * reads().readCount());
    const span<const AnchorId> journey = journeys()[orientedReadId];

    // Page title.
    html << "<h2>Journey of oriented read " << orientedReadId << "</h2>";

    // Begin the main table.
    html <<
        "<table><tr>"
        "<th>Position<br>in journey"
        "<th>Anchor"
        "<th>Anchor<br>coverage"
        "<th>Marker<br>ordinal"
        "<th>Marker<br>position";

    // Loop over the anchors in the journey of this oriented read.
    for(uint64_t positionInJourney=0; positionInJourney<journey.size(); positionInJourney++) {
        const AnchorId anchorId = journey[positionInJourney];
        const uint64_t anchorCoverage = anchors()[anchorId].coverage();
        const string anchorIdString = anchorIdToString(anchorId);

        const uint64_t ordinal = anchors().getOrdinal(anchorId, orientedReadId);

        const auto orientedReadMarkers = markers()[orientedReadId.getValue()];
        const uint32_t position = orientedReadMarkers[ordinal].position;

        if(position < beginPosition) {
            continue;
        }
        if(position >= endPosition) {
            continue;
        }

        html <<
            "<tr>"
            "<td class=centered>" << positionInJourney <<
            "<td class=centered>" <<
            "<a href='exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "'>" <<
            anchorIdString << "</a>"
            "<td class=centered>" << anchorCoverage <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position;
    }

    html << "</table>";
}





void Assembler::exploreReadFollowing(const vector<string>& request, ostream& html)
{
    // Get the request parameters.
    string anchorIdString;
    const bool anchorIdStringIsPresent = HttpServer::getParameterValue(request, "anchorIdString", anchorIdString);
    boost::trim(anchorIdString);

    uint64_t direction = 0;
    HttpServer::getParameterValue(request, "direction", direction);

    uint64_t minCommon = 4;
    HttpServer::getParameterValue(request, "minCommon", minCommon);

    double minJaccard = 0.;
    HttpServer::getParameterValue(request, "minJaccard", minJaccard);

    double minCorrectedJaccard = 0.8;
    HttpServer::getParameterValue(request, "minCorrectedJaccard", minCorrectedJaccard);



    // Begin the form.
    html <<
        "<h2>Read following</h2>"
        "<form>"
        "<table>"

    // AnchorId
        "<tr><th class=left>Anchor id<td class=centered>"
        "<input type=text name=anchorIdString required style='text-align:center'";
    if(anchorIdStringIsPresent) {
        html << " value='" << anchorIdString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>";

    // Read following parameters
    html <<
        "<tr><th class=left>Direction"
        "<td class=centered><input type=text name=direction size=8 value='" << direction << "' style='text-align:center'>"
        "<tr><th class=left>minCommon"
        "<td class=centered><input type=text name=minCommon size=8 value='" << minCommon << "' style='text-align:center'>"
        "<tr><th class=left>minJaccard"
        "<td class=centered><input type=text name=minJaccard size=8 value='" << minJaccard << "' style='text-align:center'>"
        "<tr><th class=left>minCorrectedJaccard"
        "<td class=centered><input type=text name=minCorrectedJaccard size=8 value='" << minCorrectedJaccard << "' style='text-align:center'>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Show anchor details'> "
        "</form>";



    // If the anchor id missing or invalid, stop here.
    if(not anchorIdStringIsPresent) {
        return;
    }
    const AnchorId anchorId0 = anchorIdFromString(anchorIdString);

    if((anchorId0 == invalid<AnchorId>) or (anchorId0 >= anchors().size())) {
        html << "<p>Invalid anchor id. Must be a number between 0 and " <<
            anchors().size() / 2 - 1 << " followed by + or -.";
        return;
    }


    html << "<h2>Read following</h2>";

    // Do the read following.
    vector< pair<AnchorId, AnchorPairInfo> > anchorInfos;
    anchors().followOrientedReads(
        journeys(),
        anchorId0,
        direction,
        minCommon,
        minJaccard,
        minCorrectedJaccard,
        anchorInfos);



    // Write out the results.
    html << "<p>Found " << anchorInfos.size() << " anchors.";
    html <<
        "<p><table>"
        "<tr><th>AnchorId<th>Offset<th>Common<th>Jaccard<th>Corrected<br>Jaccard";


    for(const auto& p: anchorInfos) {
        const AnchorId anchorId1 = p.first;
        const AnchorPairInfo& info = p.second;

        html <<
            "<tr>"
            "<td class=centered>" << anchorIdToString(anchorId1) <<
            "<td class=centered>" << info.offsetInBases <<
            "<td class=centered>" << info.common <<
            "<td class=centered>" << info.jaccard() <<
            "<td class=centered>" << info.correctedJaccard();

    }

    html << "</table>";


}




void Assembler::exploreLocalAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the options that control graph creation.
    string anchorIdsString;
    HttpServer::getParameterValue(request, "anchorIdsString", anchorIdsString);
    boost::trim(anchorIdsString);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    uint64_t minCoverage = 1;
    HttpServer::getParameterValue(request, "minCoverage", minCoverage);

    string useCompleteAnchorGraphString;
    bool useCompleteAnchorGraph = HttpServer::getParameterValue(request,
        "useCompleteAnchorGraph", useCompleteAnchorGraphString);

    string includeEdgesNotMarkedForAssemblyString;
    bool includeEdgesNotMarkedForAssembly = HttpServer::getParameterValue(request,
        "includeEdgesNotMarkedForAssembly", includeEdgesNotMarkedForAssemblyString);


    // Get the options that control graph display.
    const LocalAnchorGraphDisplayOptions displayOptions(request);



    // Start the form.
    html << "<form><table>";

    // Form items for options that control graph creation.
    html <<
        "<tr title='Enter comma or space separated anchor ids, each a number between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>"
        "<th class=left>Starting anchor ids"
        "<td class=centered><input type=text name=anchorIdsString style='text-align:center' required";
    if(not anchorIdsString.empty()) {
        html << " value='" << anchorIdsString + "'";
    }
    html <<
        " size=8 title='Enter comma separated anchor ids, each between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'>";

    html << "<tr>"
        "<th class=left>Distance"
        "<td class=centered>"
        "<input type=text name=distance style='text-align:center' required size=8 value=" <<
        distance << ">";

    html << "<tr>"
        "<th class=left>Minimum coverage"
        "<td class=centered>"
        "<input type=text name=minCoverage style='text-align:center' required size=8 value=" <<
        minCoverage << ">";

    if(completeAnchorGraphPointer) {
        html << "<tr>"
            "<th class=left>Use the complete AnchorGraph"
            "<td class=centered>"
            "<input type=checkbox name=useCompleteAnchorGraph" <<
            (useCompleteAnchorGraph ? " checked" : "") <<
            ">";
    }

    html << "<tr>"
        "<th class=left>Include edges not marked for use in assembly"
        "<td class=centered>"
        "<input type=checkbox name=includeEdgesNotMarkedForAssembly" <<
        (includeEdgesNotMarkedForAssembly ? " checked" : "") <<
        ">";

    // Form items for options that control graph display.
    displayOptions.writeForm(html);

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Create local anchor graph'>"
        "</form>";



    // If the anchor id are missing, stop here.
    if(anchorIdsString.empty()) {
        return;
    }


    // Extract the AnchorIds.
    vector<AnchorId> anchorIds;
    boost::tokenizer< boost::char_separator<char> > tokenizer(anchorIdsString, boost::char_separator<char>(", "));

    for(const string& anchorIdString: tokenizer) {
        const AnchorId anchorId = anchorIdFromString(anchorIdString);

        if((anchorId == invalid<AnchorId>) or (anchorId >= anchors().size())) {
            html << "<p>Invalid anchor id " << anchorIdString << ". Must be a number between 0 and " <<
                anchors().size() / 2 - 1 << " followed by + or -.";
            return;
        }

        anchorIds.push_back(anchorId);
    }
    deduplicate(anchorIds);

    // Access the AnchorGraph we are going to use.
    const AnchorGraph& anchorGraph = useCompleteAnchorGraph ? *completeAnchorGraphPointer : *anchorGraphPointer;


    // If needed, get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor* assemblyGraphPointer = 0;
    if(displayOptions.vertexColoring == "byAssemblyAnnotations") {
        const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
            displayOptions.assemblyStage,
            *httpServerData.options);
        assemblyGraphPointer = &assemblyGraph;
    }



    // Create the LocalAnchorGraph starting from these AnchorIds and moving
    // away up to the specified distance.
    LocalAnchorGraph graph(
        anchors(),
        anchorGraph,
        anchorIds,
        distance,
		minCoverage,
        not includeEdgesNotMarkedForAssembly);

    html << "<h1>Local anchor graph</h1>";
    html << "The local anchor graph has " << num_vertices(graph) <<
         " vertices and " << num_edges(graph) << " edges.";

    // Write it to html.
    graph.writeHtml(html, displayOptions, assemblyGraphPointer);

}



void Assembler::exploreLocalReadAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    LocalReadAnchorGraph graph(anchors(), journeys(), request, html);
}



void Assembler::exploreLocalReadGraph(
    const vector<string>& request,
    ostream& html)
{
    if(not readGraphPointer) {
        throw runtime_error("The ReadGraph is not available.");
    }
    LocalReadGraph graph(readGraph(), request, html);
}
