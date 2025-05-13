// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "AssemblyGraph2Postprocessor.hpp"
#include "LocalAssembly.hpp"
#include "LocalAssembly1.hpp"
#include "LocalAssembly2.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>



void Assembler::exploreLocalAssembly(
    const vector<string>& request,
    ostream& html)
{
    html << "<h2>LocalAssembly</h2>";

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



void Assembler::exploreSegments(
    const vector<string>& request,
    ostream& html)
{
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    html << "<h2>Assembly graph segments</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "</table>"
        "<input type=submit value='Get information'>"
        "</form>";

    if(assemblyStage.empty()) {
        return;
    }

    // Get the AssemblyGraph2 for this assembly stage.
    const AssemblyGraph2Postprocessor& assemblyGraph2 = getAssemblyGraph2(
        assemblyStage,
        *httpServerData.assemblerOptions);

    html <<
        "<h2>Assembly graph at stage " << assemblyStage << " </h2>"
        "<p>The assembly graph at stage " << assemblyStage <<
        " has " << num_vertices(assemblyGraph2) << " vertices (segments) and " <<
        num_edges(assemblyGraph2) << " edges (links)." << endl;

    html << "<table><tr><th>Vertex<br>(segment)<br>id<th>Number<br>of<br>steps<th>Estimated<br>length<th>Actual<br>length";

    BGL_FORALL_VERTICES(v, assemblyGraph2, AssemblyGraph2) {
        const AssemblyGraph2Vertex& vertex = assemblyGraph2[v];
        html <<
            "<tr>"
            "<td class=centered>" << vertex.id <<
            "<td class=centered>" << vertex.size() <<
            "<td class=centered>" << vertex.offset() <<
            "<td class=centered>";
        if(vertex.wasAssembled) {
            html << vertex.sequenceLength();
        }
    }

    html << "</table>";
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

    string stepBeginString;
    HttpServer::getParameterValue(request, "stepBegin", stepBeginString);

    string stepEndString;
    HttpServer::getParameterValue(request, "stepEnd", stepEndString);

    string firstStepsCountString = "5";
    HttpServer::getParameterValue(request, "firstStepsCount", firstStepsCountString);

    string lastStepsCountString = "5";
    HttpServer::getParameterValue(request, "lastStepsCount", lastStepsCountString);

    string showSequenceString;
    const bool showSequence = getParameterValue(request, "showSequence", showSequenceString);

    string showSequenceDetailsString;
    const bool showSequenceDetails = getParameterValue(request, "showSequenceDetails", showSequenceDetailsString);

    uint64_t sequenceBegin = 0;
    const bool sequenceBeginIsPresent = getParameterValue(request, "sequenceBegin", sequenceBegin);

    uint64_t sequenceEnd = 0;
    const bool sequenceEndIsPresent = getParameterValue(request, "sequenceEnd", sequenceEnd);


    // Start the form.
    html << "<h2>Assembly graph segment</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

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
        "<input type=text name=stepBegin size=8 style='text-align:center' value='" << stepBeginString << "'> to "
        "<input type=text name=stepEnd size=8 style='text-align:center' value='" << stepEndString << "'>"

        "<br><input type=radio required name=displaySteps value='first'" <<
        (displaySteps == "first" ? " checked=on" : "") << "> First "
        "<input type=text name=firstStepsCount size=8 style='text-align:center' value='" << firstStepsCountString << "'>"
        " steps"

        "<br><input type=radio required name=displaySteps value='last'" <<
        (displaySteps == "last" ? " checked=on" : "") << "> Last "
        "<input type=text name=lastStepsCount size=8 style='text-align:center' value='" << lastStepsCountString << "'>"
        " steps"

        "<br><input type=checkbox name=showSequenceDetails" << (showSequenceDetails ? " checked" : "") <<
        "> Also show sequence details for these steps";
        ;

    // Options to control the sequence display.
    html <<
        "<tr>"
        "<th class=left>Show sequence"
        "<td class=left>"
        "<input type=checkbox name=showSequence" << (showSequence ? " checked" : "") << "> Show sequence";

    html << "<br><input type=text name=sequenceBegin size=8";
    if(sequenceBeginIsPresent) {
        html << " value=" << sequenceBegin;
    }
    html << "> Begin position";

    html << "<br><input type=text name=sequenceEnd size=8";
    if(sequenceEndIsPresent) {
        html << " value=" << sequenceEnd;
    }
    html << "> End position";



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
    const AssemblyGraph2Postprocessor& assemblyGraph2 = getAssemblyGraph2(
        assemblyStage,
        *httpServerData.assemblerOptions);

    // Find the AssemblyGraph2Vertex corresponding to the requested segment.
    auto it = assemblyGraph2.vertexMap.find(segmentId);
    if(it == assemblyGraph2.vertexMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph2::vertex_descriptor v = it->second;
    const AssemblyGraph2Vertex& vertex = assemblyGraph2[v];



    html << "<h2>Segment " << segmentId << " at assembly stage " << assemblyStage << "</h2>";

    // Summary table.
    html <<
        "<table>"
        "<tr><th class=left>First anchor<td class = centered>" << anchorIdToString(vertex.front().anchorPair.anchorIdA) <<
        "<tr><th class=left>Last anchor<td class = centered>" << anchorIdToString(vertex.back().anchorPair.anchorIdB) <<
        "<tr><th class=left>Number of steps<td class = centered>" << vertex.size() <<
        "<tr><th class=left>Estimated length<td class = centered>" << vertex.offset();
    if(vertex.wasAssembled) {
        html <<
            "<tr><th class=left>Assembled length<td class = centered>" << vertex.sequenceLength();

    }
    html << "</table>";


    // Figure out the step position range to use.
    uint64_t stepBegin = invalid<uint64_t>;
    uint64_t stepEnd = invalid<uint64_t>;
    if(displaySteps == "all") {
        stepBegin = 0;
        stepEnd = vertex.size();
    } else if(displaySteps == "range") {
        try {
            stepBegin = atoul(stepBeginString);
        } catch(std::exception& e) {
            throw runtime_error("Begin " + stepBeginString + " is not valid. Must be a number.");
        }
        try {
            stepEnd = atoul(stepEndString);
        } catch(std::exception& e) {
            throw runtime_error("End " + stepEndString + " is not valid. Must be a number.");
        }
        if(stepBegin >= vertex.size()) {
            stepBegin = vertex.size() - 1;
        }
        if(stepBegin > vertex.size()) {
            stepBegin = vertex.size();
        }
        if(stepEnd < stepBegin) {
            stepEnd = stepBegin + 1;
        }
        if(stepEnd > vertex.size()) {
            stepEnd = vertex.size();
        }
    } else if(displaySteps == "first") {
        stepBegin = 0;
        try {
            stepEnd = atoul(firstStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("First anchors count " + firstStepsCountString + " is not valid. Must be a number.");
        }
        if(stepEnd > vertex.size()) {
            stepEnd = vertex.size();
        }
    } else if(displaySteps == "last") {
        stepEnd = vertex.size();
        uint64_t count = invalid<uint64_t>;
        try {
            count = atoul(lastStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("Last anchors count " + lastStepsCountString + " is not valid. Must be a number.");
        }
        if(count > vertex.size()) {
            stepBegin = 0;
        } else {
            stepBegin = stepEnd - count;
        }
    }



    // Details table showing the requested anchors.
    if(displaySteps != "none") {

        html <<
            "<p>"
            "<table>"
            "<tr><th>Step<th>AnchorIdA<th>AnchorIdB<th>Coverage<th>Estimated<br>Length";
        if(vertex.wasAssembled) {
            html << "<th>Actual<br>Length";
            if(showSequenceDetails) {
                html <<
                    "<th>Sequence<br>begin"
                    "<th>Sequence<br>end"
                    "<th>Sequence";
            }
        }

        uint64_t sequencePosition = 0;
        if(showSequenceDetails) {
            for(uint64_t i=0; i<stepBegin; i++) {
                sequencePosition += vertex[i].sequence.size();
            }
        }

        for(uint64_t stepId=stepBegin; stepId!=stepEnd; ++stepId) {
            const AssemblyGraph2VertexStep& step = vertex[stepId];

            html <<
                "<tr>"
                "<td class=centered>" << stepId <<
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<
                "<td class=centered>" << step.anchorPair.orientedReadIds.size() <<
                "<td class=centered>" << step.offset;
            if(vertex.wasAssembled) {
                html << "<td class=centered>" << step.sequence.size();
                if(showSequenceDetails) {
                    html <<
                        "<td class=centered>" << sequencePosition <<
                        "<td class=centered>" << sequencePosition + step.sequence.size() <<
                        "<td style='font-family:monospace'>";
                    copy(step.sequence.begin(), step.sequence.end(), ostream_iterator<Base>(html));
                    sequencePosition += step.sequence.size();
                }
            }


        }

        html << "</table>";
    }



    // Sequence, if requested.
    if(showSequence and vertex.wasAssembled) {
        html << "<h2>Assembled sequence</h2>";

        vector<Base> sequence;
        vertex.getSequence(sequence);

        if(not sequenceBeginIsPresent) {
            sequenceBegin = 0;
        }
        if(not sequenceEndIsPresent) {
            sequenceEnd = sequence.size();
        }

        if(sequenceBegin >= sequence.size()) {
            sequenceBegin = sequence.size() - 1;
        }
        if(sequenceEnd > sequence.size()) {
            sequenceEnd = sequence.size();
        }
        if(sequenceEnd < sequenceBegin) {
            sequenceEnd = sequenceBegin;
        }

        html << "<div style='font-family:monospace'>";
        html << ">" << segmentName << "-" << sequenceBegin << "-" << sequenceEnd <<
            ", length " << sequenceEnd - sequenceBegin << "<br>";
        copy(sequence.begin() + sequenceBegin, sequence.begin() + sequenceEnd,
            ostream_iterator<Base>(html));
        html << "</div>";

        // Also write the sequence to LocalAssembly.fasta.
        ofstream fasta("LocalAssembly.fasta");
        fasta << ">" << segmentName << "-" << sequenceBegin << "-" << sequenceEnd <<
            " length " << sequenceEnd - sequenceBegin << endl;
        copy(sequence.begin() + sequenceBegin, sequence.begin() + sequenceEnd,
            ostream_iterator<Base>(fasta));
    }

#if 0
    // Sequence details, if requested.
    if(showSequenceDetails and (displayAnchors != "none")) {
        html <<
            "<h2>Sequence assembly details</h2>"
            "<table>"
            "<tr>"
            "<th>Step<th>AnchorIdA<th>AnchorIdB"
            "<th>CoverageA<th>CoverageB<th>Common<br>coverage<th>"
            "Sequence<br>Length<th>Position<br>begin<th>Position<br>end<th>Sequence";

        uint64_t positionBegin = 0;
        for(uint64_t step=0; step<edge.sequences.size(); step++) {
            const vector<Base>& sequence = edge.sequences[step];
            const uint64_t sequenceLength = sequence.size();
            const uint64_t positionEnd = positionBegin + sequenceLength;

            if((step >= begin) and (step + 1 < end)) {
                const AnchorId anchorIdA = edge[step];
                const AnchorId anchorIdB = edge[step + 1];
                const uint64_t coverageA = anchors()[anchorIdA].coverage();
                const uint64_t coverageB = anchors()[anchorIdB].coverage();
                const uint64_t commonCount = anchors().countCommon(anchorIdA, anchorIdB);

                html <<
                    "<tr>"

                    "<td class=centered>" <<

                    "<a href='"
                    "exploreLocalAssembly2?anchorIdAString=" << HttpServer::urlEncode(anchorIdToString(anchorIdA)) <<
                    "&anchorIdBString=" << HttpServer::urlEncode(anchorIdToString(anchorIdB)) <<
                    "'>" << step << "</a>"

                    "<td class=centered>" << anchorIdToString(anchorIdA) <<
                    "<td class=centered>" << anchorIdToString(anchorIdB) <<
                    "<td class=centered>" << coverageA <<
                    "<td class=centered>" << coverageB <<
                    "<td class=centered>" << commonCount <<
                    "<td class=centered>" << sequenceLength <<
                    "<td class=centered>" << positionBegin <<
                    "<td class=centered>" << positionEnd <<
                    "<td style='font-family:monospace'>";
                copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(html));
            }

            positionBegin = positionEnd;
        }

        html << "</table>";

    }
#endif
}



void Assembler::exploreSegmentStep(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string stepIdString;
    const bool stepIdStringIsPresent = HttpServer::getParameterValue(request, "stepId", stepIdString);
    boost::trim(stepIdString);

    string showAlignmentString;
    const bool showAlignment = getParameterValue(request, "showAlignment", showAlignmentString);

    string debugString;
    const bool debug = getParameterValue(request, "debug", debugString);


    // Start the form.
    html << "<h2>Assembly graph segment step</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>Segment step id"
        "<td class=centered><input type=text name=stepId style='text-align:center' required";
    if(stepIdStringIsPresent) {
        html << " value='" << stepIdString + "'";
    }
    html <<
        ">"

        "<tr><th>Show the alignment<td class=centered><input type=checkbox name=showAlignment" <<
        (showAlignment ? " checked" : "") << ">"

        "<tr><th>Show debug information<td class=centered><input type=checkbox name=debug" <<
        (debug ? " checked" : "") << ">";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Get segment step information'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }
    if(not stepIdStringIsPresent) {
        return;
    }

    uint64_t stepId;
    try {
        stepId = atoul(stepIdString);
    } catch(std::exception& e) {
        throw runtime_error("Step id " + stepIdString + " is not valid. Must be a number.");
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
    const AssemblyGraph2Postprocessor& assemblyGraph2 = getAssemblyGraph2(
        assemblyStage,
        *httpServerData.assemblerOptions);

    // Find the AssemblyGraph2Vertex corresponding to the requested segment.
    auto it = assemblyGraph2.vertexMap.find(segmentId);
    if(it == assemblyGraph2.vertexMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }

    const AssemblyGraph2::vertex_descriptor v = it->second;
    const AssemblyGraph2Vertex& vertex = assemblyGraph2[v];

    if(stepId >= vertex.size()) {
        html << "<p>Step " << stepId << " is not valid for this segment, which has " <<
            vertex.size() << " steps.";
        return;
    }

    html << "<h2>Step " << stepId << " of segment " << segmentId << " at assembly stage " <<
        assemblyStage << "</h2>";



    // Do the local assembly for this step.
    LocalAssembly2 localAssembly(
        anchors(), html, debug,
        httpServerData.assemblerOptions->aDrift,
        httpServerData.assemblerOptions->bDrift,
        vertex[stepId].anchorPair);
    localAssembly.run(showAlignment, httpServerData.assemblerOptions->localAssemblyOptions.maxAbpoaLength);



    // Also output the sequence to fasta.
    vector<Base> sequence;
    localAssembly.getSequence(sequence);

    ofstream fasta("LocalAssembly.fasta");
    fasta << ">LocalAssembly " << sequence.size() << endl;
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
}




AssemblyGraphPostprocessor& Assembler::getAssemblyGraph(
    const string& assemblyStage,
    const AssemblerOptions& assemblerOptions)
{
    auto it = assemblyGraphTable.find(assemblyStage);
    if(it == assemblyGraphTable.end()) {
        shared_ptr<AssemblyGraphPostprocessor> p =
            make_shared<AssemblyGraphPostprocessor>(assemblerOptions, anchors(), assemblyStage);
        tie(it, ignore) = assemblyGraphTable.insert(make_pair(assemblyStage, p));
    }
    return *(it->second);
}



AssemblyGraph2Postprocessor& Assembler::getAssemblyGraph2(
    const string& assemblyStage,
    const AssemblerOptions& assemblerOptions)
{
    auto it = assemblyGraph2Table.find(assemblyStage);
    if(it == assemblyGraph2Table.end()) {
        shared_ptr<AssemblyGraph2Postprocessor> p =
            make_shared<AssemblyGraph2Postprocessor>(anchors(), assemblerOptions, assemblyStage);
        tie(it, ignore) = assemblyGraph2Table.insert(make_pair(assemblyStage, p));
    }
    return *(it->second);
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
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.assemblerOptions);



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
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.assemblerOptions);

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
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.assemblerOptions);

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



void Assembler::exploreLocalAssembly1(
    const vector<string>& request,
    ostream& html)
{
    html << "<h2>LocalAssembly1</h2>";

    // Get the parameters for the request.
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    boost::trim(anchorIdAString);

    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);
    boost::trim(anchorIdBString);

    string showAlignmentString;
    const bool showAlignment = getParameterValue(request, "showAlignment", showAlignmentString);

    uint64_t maxAbpoaLength = httpServerData.assemblerOptions->localAssemblyOptions.maxAbpoaLength;
    getParameterValue(request, "maxAbpoaLength", maxAbpoaLength);


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
        anchors().size() / 2 - 1 << " followed by + or -.'><br>"

        "<tr><th>Show the alignment<td class=centered><input type=checkbox name=showAlignment" << (showAlignment ? " checked" : "") <<
        ">"

        "<tr><th>Maximum length for abpoa<br>(switch to poasta above that)"
        "<td class=centered><input type=text name=maxAbpoaLength size=8 value=" << maxAbpoaLength << ">"

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


    LocalAssembly1 localAssembly(
        anchors(), anchorIdA, anchorIdB,
        showAlignment,
        maxAbpoaLength,
        html);
}



void Assembler::exploreLocalAssembly2(
    const vector<string>& request,
    ostream& html)
{
    html << "<h2>LocalAssembly2</h2>";

    // Get the parameters for the request.
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    boost::trim(anchorIdAString);

    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);
    boost::trim(anchorIdBString);

    string showAlignmentString;
    const bool showAlignment = getParameterValue(request, "showAlignment", showAlignmentString);

    string debugString;
    const bool debug = getParameterValue(request, "debug", debugString);

    uint64_t maxAbpoaLength = httpServerData.assemblerOptions->localAssemblyOptions.maxAbpoaLength;
    getParameterValue(request, "maxAbpoaLength", maxAbpoaLength);


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
        anchors().size() / 2 - 1 << " followed by + or -.'><br>"

        "<tr><th>Show the alignment<td class=centered><input type=checkbox name=showAlignment" << (showAlignment ? " checked" : "") <<
        ">"

        "<tr><th>Show debug information<td class=centered><input type=checkbox name=debug" << (debug ? " checked" : "") <<
        ">"

        "<tr><th>Maximum length for abpoa<br>(switch to poasta above that)"
        "<td class=centered><input type=text name=maxAbpoaLength size=8 value=" << maxAbpoaLength << ">"

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



    // Write to fasta the oriented read sequences to be used in this local Assembly.
    {
        const AnchorPair anchorPair(anchors(), anchorIdA, anchorIdB, false);
        vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
        vector< vector<Base> > sequences;
        anchorPair.get(anchors(), positions, sequences);
        ofstream fasta("LocalAssemblyInput.fasta");
        for(uint64_t i=0; i<anchorPair.orientedReadIds.size(); i++) {
            fasta << ">" << anchorPair.orientedReadIds[i] << " " << i << " " <<
                positions[i].first.basePosition << "-" <<
                positions[i].second.basePosition << " " <<
                positions[i].second.basePosition - positions[i].first.basePosition<< "\n";
           copy(sequences[i].begin(), sequences[i].end(), ostream_iterator<Base>(fasta));
           fasta << "\n";
        }

    }



    LocalAssembly2 localAssembly(
        anchors(),
        html,
        debug,
        httpServerData.assemblerOptions->aDrift,
        httpServerData.assemblerOptions->bDrift,
        anchorIdA, anchorIdB);
    localAssembly.run(showAlignment, httpServerData.assemblerOptions->localAssemblyOptions.maxAbpoaLength);



    vector<Base> sequence;
    localAssembly.getSequence(sequence);

    ofstream fasta("LocalAssembly.fasta");
    fasta << ">LocalAssembly " << sequence.size() << endl;
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
}
