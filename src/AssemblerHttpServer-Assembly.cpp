// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "Anchor.hpp"
#include "LocalAssembly.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

// Standard library.
#include <tuple.hpp>



void Assembler::exploreLocalAssembly(
    const vector<string>& request,
    ostream& html)
{

    const auto& localAssemblyOptions =
        httpServerData.assemblerOptions->localAssemblyOptions;

    LocalAssemblyDisplayOptions options(html);

    // Get the parameters for the request.
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    boost::trim(anchorIdAString);

    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);
    boost::trim(anchorIdBString);

    string useAString;
    const bool useA = HttpServer::getParameterValue(request, "useA", useAString);

    string useBString;
    const bool useB = HttpServer::getParameterValue(request, "useB", useBString);

    uint64_t minVertexCoverage = 0;
    HttpServer::getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    string showOrientedReadsString;
    options.showOrientedReads = HttpServer::getParameterValue(request, "showOrientedReads", showOrientedReadsString);

    string showMarkersString;
    options.showMarkers = HttpServer::getParameterValue(request, "showMarkers", showMarkersString);

    string showGraphString;
    options.showGraph = HttpServer::getParameterValue(request, "showGraph", showGraphString);

    string showVerticesString;
    options.showVertices = HttpServer::getParameterValue(request, "showVertices", showVerticesString);

    string showVertexLabelsString;
    options.showVertexLabels = HttpServer::getParameterValue(request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    options.showEdgeLabels = HttpServer::getParameterValue(request, "showEdgeLabels", showEdgeLabelsString);

    string showAssemblyDetailsString;
    options.showAssemblyDetails = HttpServer::getParameterValue(request, "showAssemblyDetails", showAssemblyDetailsString);

    string showDebugInformationString;
    options.showDebugInformation = HttpServer::getParameterValue(request, "showDebugInformation", showDebugInformationString);



    // Write the form.
    html <<
        "<form>"
        "<table>";

    html <<
        "<tr><th class=left>Anchor A"
        "<td class=centered><input type=text name=anchorIdAString required";
    if(anchorIdAStringIsPresent) {
        html << " value='" << anchorIdAString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr><th class=left>Anchor B"
        "<td class=centered><input type=text name=anchorIdBString required";
    if(anchorIdBStringIsPresent) {
        html << " value='" << anchorIdBString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors().size() / 2 - 1 << " followed by + or -.'><br>";

    html <<
        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on anchor A"
        "<td class=centered><input type=checkbox name=useA" <<
        (useA ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Use for assembly oriented reads that appear only on anchor B"
        "<td class=centered><input type=checkbox name=useB" <<
        (useB ? " checked" : "") << ">"

        "<tr><th class=left>Minimum vertex coverage<br>(0 = automatic)<td class=centered>"
        "<input type=text required name=minVertexCoverage size=8 style='text-align:center' "
        "value='" << minVertexCoverage << "'>"

        "<tr>"
        "<th class=left>Display the oriented reads"
        "<td class=centered><input type=checkbox name=showOrientedReads" <<
        (options.showOrientedReads ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the markers"
        "<td class=centered><input type=checkbox name=showMarkers" <<
        (options.showMarkers ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the graph"
        "<td class=centered><input type=checkbox name=showGraph" <<
        (options.showGraph ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display the vertices"
        "<td class=centered><input type=checkbox name=showVertices" <<
        (options.showVertices ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display vertex labels"
        "<td class=centered><input type=checkbox name=showVertexLabels" <<
        (options.showVertexLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display edge labels"
        "<td class=centered><input type=checkbox name=showEdgeLabels" <<
        (options.showEdgeLabels ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display assembly details"
        "<td class=centered><input type=checkbox name=showAssemblyDetails" <<
        (options.showAssemblyDetails ? " checked" : "") << ">"

        "<tr>"
        "<th class=left>Display debug information"
        "<td class=centered><input type=checkbox name=showDebugInformation" <<
        (options.showDebugInformation ? " checked" : "") << ">"

        "</table>"
        "<br><input type=submit value='Do it'>"
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

    // Local assembly for this assembly step.
    LocalAssembly localAssembly(
        assemblerInfo->k, reads(), markers(), anchors(),
        anchorIdA, anchorIdB, minVertexCoverage,
        options,
        localAssemblyOptions,
        useA, useB);
}



void Assembler::exploreSegment(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string displaySteps = "none";
    HttpServer::getParameterValue(request, "displaySteps", displaySteps);

    string beginString;
    HttpServer::getParameterValue(request, "begin", beginString);

    string endString;
    HttpServer::getParameterValue(request, "end", endString);

    string firstStepsCountString = "5";
    HttpServer::getParameterValue(request, "firstStepsCountString", firstStepsCountString);

    string lastStepsCountString = "5";
    HttpServer::getParameterValue(request, "lastStepsCount", lastStepsCountString);


    // Start the form.
    html << "<h2>Assembly graph segment</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";



    // Options to control which segment steps are shown.
    html <<
        "<tr>"
        "<th class=left>Show segment steps"
        "<td class=left>"

        "<input type=radio required name=displaySteps value='none'" <<
        (displaySteps == "none" ? " checked=on" : "") << "> None"

        "<br><input type=radio required name=displaySteps value='all'" <<
        (displaySteps == "all" ? " checked=on" : "") << "> All"

        "<br><input type=radio required name=displaySteps value='range'" <<
        (displaySteps == "range" ? " checked=on" : "") << "> Steps in position range "
        "<input type=text name=begin size=8 style='text-align:center' value='" << beginString << "'> to "
        "<input type=text name=end size=8 style='text-align:center' value='" << endString << "'>"

        "<br><input type=radio required name=displaySteps value='first'" <<
        (displaySteps == "first" ? " checked=on" : "") << "> First "
        "<input type=text name=firstStepsCount size=8 style='text-align:center' value='" << firstStepsCountString << "'>"
        " steps"

        "<br><input type=radio required name=displaySteps value='last'" <<
        (displaySteps == "last" ? " checked=on" : "") << "> Last "
        "<input type=text name=lastStepsCount size=8 style='text-align:center' value='" << lastStepsCountString << "'>"
        " steps"
        ;

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Get segment information'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    const AnchorId anchorIdA = edge.front().anchorPair.anchorIdA;
    const AnchorId anchorIdB = edge.back().anchorPair.anchorIdB;
    SHASTA_ASSERT(anchorIdA == assemblyGraph[source(e, assemblyGraph)].anchorId);
    SHASTA_ASSERT(anchorIdB == assemblyGraph[target(e, assemblyGraph)].anchorId);

    html << "<h2>Segment " << segmentId << " at assembly stage " << assemblyStage << "</h2>";

    // Summary table.
    html <<
        "<table>"
        "<tr><th class=left>First anchor<td class = centered>" << anchorIdToString(anchorIdA) <<
        "<tr><th class=left>Last anchor<td class = centered>" << anchorIdToString(anchorIdB) <<
        "<tr><th class=left>Number of steps<td class = centered>" << edge.size() <<
        "<tr><th class=left>Estimated length<td class = centered>" << edge.length() <<
        "<tr><th class=left>Lower limit on estimated length<td class = centered>" << edge.minOffset <<
        "<tr><th class=left>Upper limit on estimated length<td class = centered>" << edge.maxOffset <<
        "</table>";

    // Details table showing the requested steps.
    if(displaySteps == "none") {
        return;
    }



    // Figure out the step position range to use.
    uint64_t begin = invalid<uint64_t>;
    uint64_t end = invalid<uint64_t>;
    if(displaySteps == "all") {
        begin = 0;
        end = edge.size();
    } else if(displaySteps == "range") {
        try {
            begin = atoul(beginString);
        } catch(std::exception& e) {
            throw runtime_error("Begin " + beginString + " is not valid. Must be a number.");
        }
        try {
            end = atoul(endString);
        } catch(std::exception& e) {
            throw runtime_error("End " + endString + " is not valid. Must be a number.");
        }
        if(begin >= edge.size()) {
            begin = edge.size() - 1;
        }
        if(end > edge.size()) {
            end = edge.size();
        }
        if(end < begin) {
            end = begin + 1;
        }
    } else if(displaySteps == "first") {
        begin = 0;
        try {
            end = atoul(firstStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("First steps count " + firstStepsCountString + " is not valid. Must be a number.");
        }
        if(end > edge.size()) {
            end = edge.size();
        }
    } else if(displaySteps == "last") {
        end = edge.size();
        uint64_t count = invalid<uint64_t>;
        try {
            count = atoul(lastStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("Last steps count " + lastStepsCountString + " is not valid. Must be a number.");
        }
        if(count > edge.size()) {
            begin = 0;
        } else {
            begin = end - count;
        }
    } else {
        SHASTA_ASSERT(0);
    }
    SHASTA_ASSERT(end > begin);


    html << "<p><table>"
        "<tr><th>Position<th>AnchorIdA<th>AnchorIdB<th>Coverage"
        "<th>Average<br>offset"
        "<th>Min<br>offset"
        "<th>Max<br>offset";

    for(uint64_t i=begin; i!=end; i++) {
        const AssemblyGraphStep& step = edge[i];
        const string url =
            "exploreSegmentStep?"
            "assemblyStage=" + assemblyStage +
            "&segmentName=" + segmentName +
            "&positionInSegment=" + to_string(i);

        html <<
            "<tr>"
            "<td class=centered>" <<
            "<a href='" << url << "'>" << i << "</a>" <<
            "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
            "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<
            "<td class=centered>" << step.anchorPair.orientedReadIds.size() <<
            "<td class=centered>" << step.averageOffset <<
            "<td class=centered>" << step.minOffset <<
            "<td class=centered>" << step.maxOffset;
    }

    html << "</table>";
}


AssemblyGraphPostprocessor& Assembler::getAssemblyGraph(const string& assemblyStage)
{
    auto it = assemblyGraphTable.find(assemblyStage);
    if(it == assemblyGraphTable.end()) {
        shared_ptr<AssemblyGraphPostprocessor> p =
            make_shared<AssemblyGraphPostprocessor>(anchors(), assemblyStage);
        tie(it, ignore) = assemblyGraphTable.insert(make_pair(assemblyStage, p));
    }
    return *(it->second);
}



void Assembler::exploreSegmentStep(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    uint64_t positionInSegment;
    HttpServer::getParameterValue(request, "positionInSegment", positionInSegment);

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    if(positionInSegment >= edge.size()) {
        html << "<p>Invalid position in segment.";
        return;
    }

    const AssemblyGraphStep& step = edge[positionInSegment];

    html << "<h2>Step " << positionInSegment << " for segment " << segmentId <<
        " at assembly stage " << assemblyStage << "</h2>";

    // Summary table.
    html <<
        "<p><table>"
        "<tr><th class=left>AnchorIdA<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
        "<tr><th class=left>AnchorIdB<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<
        "<tr><th class=left>Coverage<td class=centered>" << step.anchorPair.orientedReadIds.size() <<
        "<tr><th class=left>Average offset<td class=centered>" << step.averageOffset <<
        "<tr><th class=left>Min offset<td class=centered>" << step.minOffset <<
        "<tr><th class=left>max offset<td class=centered>" << step.maxOffset <<
        "</table>";



    // Details table.
    using Positions = AnchorPair::Positions;
    vector< pair<Positions, Positions> > positions;
    step.anchorPair.get(anchors(), positions);
    html <<
        "<p><table><tr>"
        "<th class=centered>Oriented<br>read id"
        "<th class=centered>Position<br>in<br>journey<br>A"
        "<th class=centered>Position<br>in<br>journey<br>B"
        "<th class=centered>Journey<br>offset"
        "<th class=centered>OrdinalA"
        "<th class=centered>OrdinalB"
        "<th class=centered>Ordinal<br>offset"
        "<th class=centered>A<br>middle<br>position"
        "<th class=centered>B<br>middle<br>position"
        "<th class=centered>Sequence<br>length"
        ;
    for(uint64_t i=0; i<step.anchorPair.orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = step.anchorPair.orientedReadIds[i];
        const auto& p = positions[i];
        const Positions& positionsA = p.first;
        const Positions& positionsB = p.second;
        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << positionsA.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney - positionsA.positionInJourney <<
            "<td class=centered>" << positionsA.ordinal <<
            "<td class=centered>" << positionsB.ordinal <<
            "<td class=centered>" << positionsB.ordinal - positionsA.ordinal <<
            "<td class=centered>" << positionsA.basePosition <<
            "<td class=centered>" << positionsB.basePosition <<
            "<td class=centered>" << positionsB.basePosition - positionsA.basePosition;

    }
    html << "</table>";
}



void Assembler::exploreTangleMatrix(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);
    boost::trim(assemblyStage);

    string entrancesString;
    HttpServer::getParameterValue(request, "entrances", entrancesString);
    boost::trim(entrancesString);

    string exitsString;
    HttpServer::getParameterValue(request, "exits", exitsString);
    boost::trim(exitsString);



    // Start the form.
    html << "<h2>Assembly graph tangle matrix</h2><form>";

    html <<
        "<table>"
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr title='Enter assembly segment names separated by commas or spaces.'>"
        "<th class=left>Entrances"
        "<td class=centered>"
        "<input type=text name=entrances style='text-align:center'";
    if(not entrancesString.empty()) {
        html << " value='" << entrancesString << "'";
    }
    html << " size=30>";

    html <<
        "<tr title='Enter assembly segment names separated by commas or spaces.'>"
        "<th class=left>Exits"
        "<td class=centered>"
        "<input type=text name=exits style='text-align:center'";
    if(not exitsString.empty()) {
        html << " value='" << exitsString << "'";
    }
    html << " size=30>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Compute tangle matrix'>"
        "</form>";

    if(entrancesString.empty() or exitsString.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage);



    // Find AssemblyGraph edges corresponding to the entrances.
    vector<AssemblyGraph::edge_descriptor> entrances;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(entrancesString, boost::char_separator<char>(", "));
        for(const string& edgeIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(edgeIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << edgeIdString << ". Must be a number.";
                return;
            }

            // Find the AssemblyGraphEdge corresponding to the requested segment.
            auto it = assemblyGraph.edgeMap.find(segmentId);
            if(it == assemblyGraph.edgeMap.end()) {
                html << "<p>Assembly graph at stage " << assemblyStage <<
                    " does not have segment " << segmentId;
                return;
            }

            entrances.push_back(it->second);
        }

    }



    // Find AssemblyGraph edges corresponding to the exits.
    vector<AssemblyGraph::edge_descriptor> exits;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(exitsString, boost::char_separator<char>(", "));
        for(const string& edgeIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(edgeIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << edgeIdString << ". Must be a number.";
                return;
            }

            // Find the AssemblyGraphEdge corresponding to the requested segment.
            auto it = assemblyGraph.edgeMap.find(segmentId);
            if(it == assemblyGraph.edgeMap.end()) {
                html << "<p>Assembly graph at stage " << assemblyStage <<
                    " does not have segment " << segmentId;
                return;
            }

            exits.push_back(it->second);
        }
    }


    // Create the TangleMatrix.
    const TangleMatrix tangleMatrix(assemblyGraph, entrances, exits);
    tangleMatrix.writeHtml(assemblyGraph, html);
}



// The vertex is specified using one of the segments that have the vertex as their target.
void Assembler::exploreVertexTangle(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    // Start the form.
    html << "<h2>Assembly graph vertex tangle</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr>"
        "<th class=left>Name of a segment with target<br>at the desired vertex"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";



   // End the form.
    html <<
        "</table>"
        "<input type=submit value='Get tangle information'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;



    // Create a Tangle at the target of this edge.
    const AssemblyGraph::vertex_descriptor v = target(e, assemblyGraph);
    const Tangle tangle(assemblyGraph, v);

    // Write out the TangleMatrix.
    tangle.tangleMatrix.writeHtml(assemblyGraph, html);
}



void Assembler::exploreEdgeTangle(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    // Start the form.
    html << "<h2>Assembly graph edge tangle</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";



   // End the form.
    html <<
        "</table>"
        "<input type=submit value='Get tangle information'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(assemblyStage);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;



    // Create a Tangle at this edge.
    const Tangle tangle(assemblyGraph, e);

    // Write out the TangleMatrix.
    tangle.tangleMatrix.writeHtml(assemblyGraph, html);
}
