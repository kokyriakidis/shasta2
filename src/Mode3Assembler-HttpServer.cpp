// Shasta.
#include "Mode3Assembler.hpp"
#include "deduplicate.hpp"
#include "HttpServer.hpp"
#include "Markers.hpp"
#include "mode3-LocalAssembly.hpp"
#include "mode3-LocalAnchorGraph.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>



namespace shasta {
    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.
}

// Boost libraries.
#include <boost/tokenizer.hpp>

// Standard library.
#include "fstream.hpp"



void Mode3Assembler::exploreLocalAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the options that control graph creation.
    string anchorIdsString;
    HttpServer::getParameterValue(request, "anchorIdsString", anchorIdsString);
    boost::trim(anchorIdsString);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    uint64_t minCoverage = 0;
    HttpServer::getParameterValue(request, "minCoverage", minCoverage);

    string filterEdgesByCoverageLossString;
    const bool filterEdgesByCoverageLoss = HttpServer::getParameterValue(request,
        "filterEdgesByCoverageLoss", filterEdgesByCoverageLossString);

    double maxCoverageLoss =  options.primaryGraphOptions.maxLoss;
    HttpServer::getParameterValue(request, "maxCoverageLoss", maxCoverageLoss);

    // Get the options that control graph display.
    const LocalAnchorGraphDisplayOptions displayOptions(request);



    // Start the form.
    html << "<form><table>";

    // Form items for options that control graph creation.
    html <<
        "<tr title='Enter comma separated anchor ids, each a number between 0 and " <<
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

    html <<
        "<tr><th>Edge filtering"
        "<td>"
        "Minimum coverage "
        "<input type=text name=minCoverage style='text-align:center' required size=8 value=" <<
        minCoverage << ">"
        "<br>"
        "<input type=checkbox name=filterEdgesByCoverageLoss" <<
        (filterEdgesByCoverageLoss ? " checked" : "") <<
        ">Filter edges by coverage loss"
        "<br><input type=text name=maxCoverageLoss style='text-align:center' required size=6 value=" <<
        maxCoverageLoss << "> Maximum coverage loss"
        "<hr>";

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
    boost::tokenizer< boost::char_separator<char> > tokenizer(anchorIdsString, boost::char_separator<char>(","));

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



    // Create the LocalAnchorGraph starting from these AnchorIds and moving
    // away up to the specified distance.
    LocalAnchorGraph graph(
        anchors(),
        anchorIds,
        distance,
        filterEdgesByCoverageLoss,
        maxCoverageLoss,
        minCoverage
        );

    html << "<h1>Local anchor graph</h1>";
    html << "The local anchor graph has " << num_vertices(graph) <<
         " vertices and " << num_edges(graph) << " edges.";

    // Write it to html.
    graph.writeHtml(html, displayOptions);

}
