// Shasta.
#include "Assembler.hpp"
#include "areSimilarSequences.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "graphvizToHtml.hpp"
#include "GTest.hpp"
#include "LocalAssembly.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
#include "TangleMatrix1.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>



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

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    html <<
        "<h2>Assembly graph at stage " << assemblyStage << " </h2>"
        "<p>The assembly graph at stage " << assemblyStage <<
        " has " << num_edges(assemblyGraph) << " edges (segments)." << endl;

    html << "<table><tr><th>Vertex<br>(segment)<br>id<th>Number<br>of<br>steps"
        "<th>Average<br>coverage<th>Estimated<br>length<th>Actual<br>length";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        const string url = "exploreSegment?assemblyStage=" + assemblyStage + "&segmentName=" + to_string(edge.id);
        html <<
            "<tr>"
            "<td class=centered><a href='" << url << "'>" << edge.id << "</a>"
            "<td class=centered>" << edge.size() <<
            "<td class=centered>" << coverage <<
            "<td class=centered>" << edge.offset() <<
            "<td class=centered>";
        if(edge.wasAssembled) {
            html << edge.sequenceLength();
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
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];
    const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));



    html << "<h2>Segment " << segmentId << " at assembly stage " << assemblyStage << "</h2>";

    // Summary table.
    html <<
        "<table>"
        "<tr><th class=left>First anchor<td class = centered>" << anchorIdToString(edge.front().anchorPair.anchorIdA) <<
        "<tr><th class=left>Last anchor<td class = centered>" << anchorIdToString(edge.back().anchorPair.anchorIdB) <<
        "<tr><th class=left>Number of steps<td class = centered>" << edge.size() <<
        "<tr><th class=left>Average coverage<td class = centered>" << coverage <<
        "<tr><th class=left>Estimated length<td class = centered>" << edge.offset() <<
        "<tr><th class=left>Assembled<td class = centered>" << (edge.wasAssembled ? "Yes" : "No");
    if(edge.wasAssembled) {
        html <<
            "<tr><th class=left>Assembled length<td class = centered>" << edge.sequenceLength();

    }
    html << "</table>";

    html <<
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(anchorIdToString(edge.front().anchorPair.anchorIdA)) <<
        "'>See the first anchor in the local anchor graph</a>"
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(anchorIdToString(edge.back().anchorPair.anchorIdB)) <<
        "'>See the last anchor in the local anchor graph</a>"
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(
            anchorIdToString(edge.front().anchorPair.anchorIdA) + " " +
            anchorIdToString(edge.back().anchorPair.anchorIdB)) <<
        "'>See the first and last anchor in the local anchor graph</a>";


    // Figure out the step position range to use.
    uint64_t stepBegin = invalid<uint64_t>;
    uint64_t stepEnd = invalid<uint64_t>;
    if(displaySteps == "all") {
        stepBegin = 0;
        stepEnd = edge.size();
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
        if(stepBegin >= edge.size()) {
            stepBegin = edge.size() - 1;
        }
        if(stepBegin > edge.size()) {
            stepBegin = edge.size();
        }
        if(stepEnd < stepBegin) {
            stepEnd = stepBegin + 1;
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "first") {
        stepBegin = 0;
        try {
            stepEnd = atoul(firstStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("First anchors count " + firstStepsCountString + " is not valid. Must be a number.");
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "last") {
        stepEnd = edge.size();
        uint64_t count = invalid<uint64_t>;
        try {
            count = atoul(lastStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("Last anchors count " + lastStepsCountString + " is not valid. Must be a number.");
        }
        if(count > edge.size()) {
            stepBegin = 0;
        } else {
            stepBegin = stepEnd - count;
        }
    }



    // Details table showing the requested steps.
    if(displaySteps != "none") {

        html <<
            "<p>"
            "<table>"
            "<tr><th>Step<th>AnchorIdA<th>AnchorIdB<th>Coverage<th>Estimated<br>Length";
        if(edge.wasAssembled) {
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
                sequencePosition += edge[i].sequence.size();
            }
        }

        for(uint64_t stepId=stepBegin; stepId!=stepEnd; ++stepId) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const string url = "exploreSegmentStep?assemblyStage=" +
                assemblyStage + "&segmentName=" + to_string(segmentId) + "&stepId=" + to_string(stepId);

            html <<
                "<tr>"
                "<td class=centered><a href='" << url << "'>" << stepId << "</a>"
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<
                "<td class=centered>" << step.anchorPair.orientedReadIds.size() <<
                "<td class=centered>" << step.offset;
            if(edge.wasAssembled) {
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

        // Link to the local anchor graph showing these anchors.
        {
            const AssemblyGraphEdgeStep& firstStep = edge[stepBegin];
            const AnchorId anchorIdA = firstStep.anchorPair.anchorIdA;
            string urlAnchors = anchorIdToString(anchorIdA);
            for(uint64_t stepId=stepBegin; stepId!=stepEnd; ++stepId) {
                const AssemblyGraphEdgeStep& step = edge[stepId];
                const AnchorId anchorIdB = step.anchorPair.anchorIdB;
                urlAnchors += " ";
                urlAnchors += anchorIdToString(anchorIdB);
            }

            html << "<br><a href='exploreLocalAnchorGraph?anchorIdsString="<< HttpServer::urlEncode(urlAnchors) <<
            "'>See these anchors in the local anchor graph</a>";
        }
    }




    // Sequence, if requested.
    if(showSequence and edge.wasAssembled) {
        html << "<h2>Assembled sequence</h2>";

        vector<Base> sequence;
        edge.getSequence(sequence);

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
        ofstream fasta("Segment.fasta");
        fasta << ">" << segmentName << "-" << sequenceBegin << "-" << sequenceEnd <<
            " length " << sequenceEnd - sequenceBegin << endl;
        copy(sequence.begin() + sequenceBegin, sequence.begin() + sequenceEnd,
            ostream_iterator<Base>(fasta));
    }

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
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }

    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    if(stepId >= edge.size()) {
        html << "<p>Step " << stepId << " is not valid for this segment, which has " <<
            edge.size() << " steps.";
        return;
    }



    // Write the AnchorPair to html.
    html << "<h2>AnchorPair for step " << stepId << " of segment " << segmentId << " at assembly stage " <<
        assemblyStage << "</h2>";
    edge[stepId].anchorPair.writeAllHtml(html, anchors(), journeys());



    // Write the local assembly to html.
    html << "<h2>Local assembly for step " << stepId << " of segment " << segmentId << " at assembly stage " <<
        assemblyStage << "</h2>";

    // Do the local assembly for this step.
    LocalAssembly localAssembly(
        anchors(), html, debug,
        httpServerData.options->aDrift,
        httpServerData.options->bDrift,
        edge[stepId].anchorPair);
    localAssembly.run(showAlignment, httpServerData.options->maxAbpoaLength);



    // Also output the sequence to fasta.
    vector<Base> sequence;
    localAssembly.getSequence(sequence);

    ofstream fasta("LocalAssembly.fasta");
    fasta << ">LocalAssembly " << sequence.size() << endl;
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
}



void Assembler::exploreBridgeSegmentSteps(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);
    boost::trim(assemblyStage);

    string segmentNameA;
    HttpServer::getParameterValue(request, "segmentNameA", segmentNameA);
    boost::trim(segmentNameA);

    uint64_t trimA = 0;
    HttpServer::getParameterValue(request, "trimA", trimA);

    string segmentNameB;
    HttpServer::getParameterValue(request, "segmentNameB", segmentNameB);
    boost::trim(segmentNameB);

    uint64_t trimB = 0;
    HttpServer::getParameterValue(request, "trimB", trimB);



    // Start the form.
    html << "<h2>Bridge between segment steps</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered colspan=2><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "<tr>"
        "<th class=left>Step A"
        "<td class=left>Segment <input type=text size=8 name=segmentNameA style='text-align:center' required";
    if(not segmentNameA.empty()) {
        html << " value='" << segmentNameA + "'";
    }
    html <<
        ">"
        "<td class=left>Trim <input type=text size=8 name=trimA style='text-align:center' required"
        " value='" << trimA << "'"
        "> steps at end";

    html <<
        "<tr>"
        "<th class=left>Step B"
        "<td class=left>Segment <input type=text size=8 name=segmentNameB style='text-align:center' required";
    if(not segmentNameB.empty()) {
        html << " value='" << segmentNameB + "'";
    }
    html <<
        ">"
        "<td class=left>Trim <input type=text size=8 name=trimB style='text-align:center' required"
        " value='" << trimB << "'"
        "> steps at beginning";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Compute'>"
        "</form>";

    if(segmentNameA.empty() or segmentNameB.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Get the segment ids.
    uint64_t segmentIdA = invalid<uint64_t>;
    try {
        segmentIdA = std::stol(segmentNameA);
    } catch(exception&) {
    }
    if(segmentIdA == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }
    uint64_t segmentIdB = invalid<uint64_t>;
    try {
        segmentIdB = std::stol(segmentNameB);
    } catch(exception&) {
    }
    if(segmentIdB == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }



    // Find the AssemblyGraphEdges corresponding to the requested segments.
    auto itA = assemblyGraph.edgeMap.find(segmentIdA);
    if(itA == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentIdA;
        return;
    }
    const AssemblyGraph::edge_descriptor eA = itA->second;
    const AssemblyGraphEdge& edgeA = assemblyGraph[eA];

    auto itB = assemblyGraph.edgeMap.find(segmentIdB);
    if(itB == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentIdB;
        return;
    }
    const AssemblyGraph::edge_descriptor eB = itB->second;
    const AssemblyGraphEdge& edgeB = assemblyGraph[eB];


    // Get the steps and the AnchorPairs.
    if(trimA >= edgeA.size()) {
        html << "Trim for segment " << segmentIdA <<
            " is too large. Segment has " << edgeA.size() << " steps.";
        return;
    }
    if(trimB >= edgeB.size()) {
        html << "Trim for segment " << segmentIdB <<
            " is too large. Segment has " << edgeB.size() << " steps.";
        return;
    }
    const AssemblyGraphEdgeStep& stepA = edgeA[edgeA.size() - 1 - trimA];
    const AssemblyGraphEdgeStep& stepB = edgeB[trimB];
    const AnchorPair& anchorPairA = stepA.anchorPair;
    const AnchorPair& anchorPairB = stepB.anchorPair;

    // Compute the bridge AnchorPair.
    const AnchorPair bridgeAnchorPair = anchors().bridge(
        anchorPairA, anchorPairB,
        httpServerData.options->aDrift,
        httpServerData.options->bDrift);
    html << "<p>" << anchorPairA.size() << " " << anchorPairB.size() << " " << bridgeAnchorPair.size();

}



AssemblyGraphPostprocessor& Assembler::getAssemblyGraph(
    const string& assemblyStage,
    const Options& options)
{
    auto it = assemblyGraphTable.find(assemblyStage);
    if(it == assemblyGraphTable.end()) {
        shared_ptr<AssemblyGraphPostprocessor> p =
            make_shared<AssemblyGraphPostprocessor>(anchors(), journeys(), options, assemblyStage);
        tie(it, ignore) = assemblyGraphTable.insert(make_pair(assemblyStage, p));
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

    uint64_t maxTrim = 0;
    HttpServer::getParameterValue(request, "maxTrim", maxTrim);

    double epsilon = httpServerData.options->detangleEpsilon;
    HttpServer::getParameterValue(request, "epsilon", epsilon);



    // Start the form.
    html << "<h2>Assembly graph tangle matrix (old)</h2><form>";

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

    html << "<tr><th class=left>Maximum number of trimmed steps"
        "<td class=centered>"
        "<input type=text name=maxTrim style='text-align:center' value='" << maxTrim << "'>";

    html << "<tr><th class=left>Epsilon for G-test evaluation"
        "<td class=centered>"
        "<input type=text name=epsilon style='text-align:center' value='" << epsilon << "'>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Compute tangle matrix'>"
        "</form>";

    if(entrancesString.empty() or exitsString.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);



    // Find AssemblyGraph edges corresponding to the entrances.
    vector<AssemblyGraph::edge_descriptor> entrances;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(entrancesString, boost::char_separator<char>(", "));
        for(const string& vertexIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(vertexIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << vertexIdString << ". Must be a number.";
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
    std::ranges::sort(entrances, assemblyGraph.orderById);



    // Find AssemblyGraph edge corresponding to the exits.
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
    std::ranges::sort(exits, assemblyGraph.orderById);



    // Create the TangleMatrix.
    TangleMatrix tangleMatrix(assemblyGraph, entrances, exits,
        maxTrim,
        httpServerData.options->aDrift,
        httpServerData.options->bDrift);
    tangleMatrix.writeHtml(assemblyGraph, html);


    // Likelihood ratio test (G test).
    vector< vector<uint64_t> > tangleMatrixCoverage;
    tangleMatrix.getTangleMatrixCoverage(tangleMatrixCoverage);
    const GTest gTest(tangleMatrixCoverage, epsilon);
    gTest.writeHtml(html);
}



void Assembler::exploreTangleMatrix1(const vector<string>& request, ostream& html)
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

    double epsilon = httpServerData.options->detangleEpsilon;
    HttpServer::getParameterValue(request, "epsilon", epsilon);



    // Start the form.
    html << "<h2>Assembly graph tangle matrix (new)</h2><form>";

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

    html << "<tr><th class=left>Epsilon for G-test evaluation"
        "<td class=centered>"
        "<input type=text name=epsilon style='text-align:center' value='" << epsilon << "'>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Compute tangle matrix'>"
        "</form>";

    if(entrancesString.empty() or exitsString.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);



    // Find AssemblyGraph edges corresponding to the entrances.
    vector<AssemblyGraph::edge_descriptor> entrances;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(entrancesString, boost::char_separator<char>(", "));
        for(const string& vertexIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(vertexIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << vertexIdString << ". Must be a number.";
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
    std::ranges::sort(entrances, assemblyGraph.orderById);



    // Find AssemblyGraph edge corresponding to the exits.
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
    std::ranges::sort(exits, assemblyGraph.orderById);


    // Make sure we have the information we need.
    if(assemblyGraph.orientedReadEdgeInformation.empty()) {
        assemblyGraph.findOrientedReadEdgeInformation();
    }

    // Compute the tangle matrix.
    const TangleMatrix1 tangleMatrix(assemblyGraph, entrances, exits, html);
    GTest gTest(tangleMatrix.tangleMatrix, epsilon);
    gTest.writeHtml(html);



    // Create a RestrictedAnchorGraph for each element of the top hypothesis
    // that is set to 1.
    if(gTest.hypotheses.empty()) {
        return;
    }
    const auto& bestConnectivityMatrix = gTest.hypotheses.front().connectivityMatrix;
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {

                const AssemblyGraph::edge_descriptor eEntrance = entrances[iEntrance];
                const AssemblyGraph::edge_descriptor eExit = exits[iExit];

                const AssemblyGraphEdge& entranceEdge = assemblyGraph[eEntrance];
                const AssemblyGraphEdge& exitEdge = assemblyGraph[eExit];

                const AnchorId entranceAnchorId = entranceEdge.back().anchorPair.anchorIdB;
                const AnchorId exitAnchorId = exitEdge.front().anchorPair.anchorIdA;

                if(entranceAnchorId == exitAnchorId) {
                    html << "<br>The two anchors are coincident.";
                } else {

                    html << "<h4>RestrictedAnchorGraph to connect entrance " <<
                        entranceEdge.id <<
                        " with exit " << exitEdge.id << "</h4>"
                        "Last AnchorId on entrance is " << anchorIdToString(entranceAnchorId) <<
                        "<br>First AnchorId on exit is " << anchorIdToString(exitAnchorId);


                    RestrictedAnchorGraph restrictedAnchorGraph(
                        anchors(), journeys(), tangleMatrix, iEntrance, iExit, html);
                    restrictedAnchorGraph.keepBetween(entranceAnchorId, exitAnchorId);
                    restrictedAnchorGraph.removeCycles();
                    restrictedAnchorGraph.keepBetween(entranceAnchorId, exitAnchorId);

                    html << "<br>The RestrictedAnchorGraph has " << num_vertices(restrictedAnchorGraph) <<
                        " vertices and " << num_edges(restrictedAnchorGraph) << " edges ";

                    // Find the longest path in the RestrictedAnchorGraph.
                    vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                    // restrictedAnchorGraph.findLongestPath(longestPath);
                    restrictedAnchorGraph.findOptimalPath(entranceAnchorId, exitAnchorId, longestPath);

                    // Write it out in Graphviz format.
                    const string uuid = to_string(boost::uuids::random_generator()());
                    const string dotFileName = tmpDirectory() + uuid + ".dot";
                    restrictedAnchorGraph.writeGraphviz(dotFileName, {entranceAnchorId, exitAnchorId});


                    // Display it in html in svg format.
                    const double timeout = 30.;
                    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
                    html << "<p>";
                    graphvizToHtml(dotFileName, "dot", timeout, options, html);
                }

            }
        }

    }
}



void Assembler::exploreLocalAssembly(
    const vector<string>& request,
    ostream& html)
{
    html << "<h2>Local assembly</h2>";

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

    uint64_t maxAbpoaLength = httpServerData.options->maxAbpoaLength;
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



    LocalAssembly localAssembly(
        anchors(),
        html,
        debug,
        httpServerData.options->aDrift,
        httpServerData.options->bDrift,
        anchorIdA, anchorIdB);
    localAssembly.run(showAlignment, httpServerData.options->maxAbpoaLength);



    vector<Base> sequence;
    localAssembly.getSequence(sequence);

    ofstream fasta("LocalAssembly.fasta");
    fasta << ">LocalAssembly " << sequence.size() << endl;
    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
}



void Assembler::exploreSimilarSequences(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string sequence0String;
    HttpServer::getParameterValue(request, "sequence0", sequence0String);
    boost::trim(sequence0String);

    string sequence1String;
    HttpServer::getParameterValue(request, "sequence1", sequence1String);
    boost::trim(sequence1String);

    // Write the form.
    html <<
        "<h2>Similar sequences</h2><form>"
        "<table>"
        "<tr>"
        "<th class=left>First sequence"
        "<td class=centered><input type=text name=sequence0 style='text-align:center;font-family:monospace' required size=100"
        " value='" << sequence0String << "'>"
        "<tr>"
        "<th class=left>Second sequence"
        "<td class=centered><input type=text name=sequence1 style='text-align:center;font-family:monospace' required size=100"
        " value='" << sequence1String << "'>"
        "</table>"
        "<input type=submit value='Analyze'>"
        "</form>";

    if(sequence0String.empty() or sequence1String.empty()) {
        return;
    }

    // Fill in the Base sequences.
    vector<Base> sequence0;
    for(const char c: sequence0String) {
        sequence0.push_back(Base::fromCharacter(c));
    }
    vector<Base> sequence1;
    for(const char c: sequence1String) {
        sequence1.push_back(Base::fromCharacter(c));
    }

    // EXPOSE WHEN CODE STABILIZES.
    const vector<uint64_t> minRepeatCount = {0, 2, 2, 2, 2, 2, 2};
    const bool areSimilar = areSimilarSequences(sequence0, sequence1, minRepeatCount, html);
    if(areSimilar) {
        html << "<p>These sequence are similar and their differences are likely caused by sequencing errors.";
    }
}
