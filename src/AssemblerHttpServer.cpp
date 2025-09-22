
// Shasta.
#include "Assembler.hpp"
#include "filesystem.hpp"
#include "MarkerKmers.hpp"
#include "tmpDirectory.hpp"
#include "Reads.hpp"
#include "ReadSummary.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <filesystem>
#include "fstream.hpp"



#define SHASTA_ADD_TO_FUNCTION_TABLE(name) httpServerData.functionTable[string("/") + #name ] = &Assembler::name



// Associate http keywords with member functions.
void Assembler::fillServerFunctionTable()
{
    httpServerData.functionTable[""]        = &Assembler::exploreSummary;
    httpServerData.functionTable["/"]       = &Assembler::exploreSummary;
    httpServerData.functionTable["/index"]  = &Assembler::exploreSummary;

    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadRaw);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreLookupRead);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadSequence);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadMarkers);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerKmers);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerKmerPair);

    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAnchor);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAnchorPair);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAnchorPair1);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAnchorPair2);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreJourney);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadFollowing);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreLocalAssembly);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSimilarSequences);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreLocalAnchorGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreRestrictedAnchorGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreLocalAssemblyGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSuperbubble);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSegments);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSegmentSequence);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSegmentSteps);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSegmentOrientedReads);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSegmentStep);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreBridgeSegmentSteps);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreTangleMatrix);
}
#undef SHASTA_ADD_TO_FUNCTION_TABLE



void Assembler::processRequest(
    const vector<string>& request,
    ostream& html,
    const BrowserInformation&)
{
    const string& keyword = request.front();

    // Look up the keyword to find the function that will process this request.
    // Note that the keyword includes the initial "/".
    const auto it = httpServerData.functionTable.find(keyword);
    if(it == httpServerData.functionTable.end()) {
        writeHtmlBegin(html);
        writeNavigation(html);
        html << "Unsupported keyword " << keyword;
        cout << "Valid keywords are:";
        for(const auto& p: httpServerData.functionTable) {
            cout << " " << p.first;
        }
        cout << endl;
        writeHtmlEnd(html);
        return;
    }


    // We found the keyword. Call the function that processes this keyword.
    // The processing function is only responsible for writing the html body.
    writeHtmlBegin(html);
    writeNavigation(html);
    try {
        const auto function = it->second;
        (this->*function)(request, html);
    } catch(const std::exception& e) {
        html << "<br><br><span style='color:purple'>" << e.what() << "</span>";
    }
    writeHtmlEnd(html);
}



void Assembler::writeMakeAllTablesCopyable(ostream& html) const
{
    html << R"###(
    <script>

    // Copy to the clipboard the table that generated the event.
    function copyToClipboard(event)
    {
        // If the CTRL key is not pressed, don't do anything.
        if(!event.ctrlKey) {
            return;
        }

        // Prevent default behavior.
        // event.preventDefault();
        // event.stopPropagation();
        // event.returnValue = false;
        
        // Get the table element.
        var element = event.currentTarget;
         
        // Remove any previous selection.
        var selection = window.getSelection();
        selection.removeAllRanges();
        
        // Select the table.
        var range = document.createRange();
        range.selectNodeContents(element);
        selection.addRange(range);
        
        // Copy it to the clipboard.
        document.execCommand("copy");

        // Unselect it.
        selection.removeAllRanges();

        window.alert("The table was copied to the clipboard");
    }

    // Make a table copyable by Ctrl-click.
    function makeCopyable(element)
    {
        element.addEventListener('click', copyToClipboard);
        element.title = 'Ctrl-click anywhere on the table to copy the entire table to the clipboard';
    }

    // Make all tables copyable by Ctrl-click.
    function makeAllTablesCopyable()
    {
        var tables = document.getElementsByTagName('table');
        var i;
        for(i=0; i<tables.length; i++) {
            makeCopyable(tables[i]);
        }
    }
    </script>
    )###";


#if 0
    html << R"###(
<script>

// Make all tables selectable by double click.
// This must be called after all tables have
// already been created, so it can be called during onload.

// This function is called when the user double clicks on a table.
function selectElement(table)
{
    var selection = window.getSelection();
    selection.removeAllRanges();
    var range = document.createRange();
    range.selectNode(table);
    selection.addRange(range);
}

// Attach the above function to the double click event
// for all tables in the document.
// Also add to each table a title that displays a tooltip 
// explaining that the table can be selected via double click.
function makeAllTablesSelectableByDoubleClick()
{
    var allTables = document.getElementsByTagName("table");
    for (var i=0; i<allTables.length; i++) {
        var table = allTables[i];
        table.ondblclick = function() {selectElement(this);};
        table.setAttribute("title", 
        "Double click to select the entire table. You can then paste it into a spreadsheet.");
    }
}
</script>
    )###";
#endif
}



void Assembler::writeNavigation(ostream& html) const
{
    html << "<ul class=navigationMenu>";

    // Reads menu.
    writeNavigation(html, "Reads", {
        {"Sequence", "exploreReadSequence"},
        {"Look up a read by name", "exploreLookupRead"},
        });



    // Markers menu.
    writeNavigation(html, "Markers", {
        {"Markers", "exploreReadMarkers"},
        {"Marker k-mers", "exploreMarkerKmers"},
        {"Marker k-mer pair", "exploreMarkerKmerPair"}
        });



    // Anchors menu.
    writeNavigation(html, "Anchors", {
        {"Anchor", "exploreAnchor"},
        {"Anchor pair", "exploreAnchorPair2"},
        {"Anchor pair (old)", "exploreAnchorPair"},
        {"Journey", "exploreJourney"},
        {"Local journey analysis", "exploreRestrictedAnchorGraph"},
        {"Read following on anchors", "exploreReadFollowing"},
        {"Local anchor graph", "exploreLocalAnchorGraph"},
        });

    // Assembly menu.
    writeNavigation(html, "Assembly", {
        {"Local assembly graph", "exploreLocalAssemblyGraph"},
        {"Superbubble", "exploreSuperbubble"},
        {"Segments", "exploreSegments"},
        {"Segment sequence", "exploreSegmentSequence"},
        {"Segment steps", "exploreSegmentSteps"},
        {"Segment oriented reads", "exploreSegmentOrientedReads"},
        {"Segment step", "exploreSegmentStep"},
        {"Bridge segment steps", "exploreBridgeSegmentSteps"},
        {"Tangle matrix", "exploreTangleMatrix"},
        {"Local assembly", "exploreLocalAssembly"},
        {"Sequence similarity", "exploreSimilarSequences"},
        });

    html << "</ul>";
}



void Assembler::writeNavigation(
    ostream& html,
    const string& title,
    const vector<pair <string, string> >& items) const
{
    html <<
        "<li class=navigationMenuEntry>"
        "<div class=navigationButton>" << title << "</div>"
        "<div class=navigationItems>";

    for(const auto& item: items) {
        html << "<a class=navigationItem href=" << item.second << ">" << item.first << "</a>";
    }

    html << "</div></li>";

}



// Write to html an img tag displaying a png file.
void Assembler::writePngToHtml(
    ostream& html,
    const string& pngFileName,
    const string useMap)
{
    // Convert the png file to base64.
    const string base64FileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    const string base64Command = "base64 " + pngFileName + " > " +
        base64FileName;
    const int errorCode = ::system(base64Command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " +
            to_string(errorCode) + " " + strerror(errorCode) +
            "\nrunning command: " + base64Command);
    }

    // Write the base64 file to html in an img tag.
    html << "<p><img ";
    if(not useMap.empty()) {
        html << "usemap='" << useMap << "'";
    }
    html << " src=\"data:image/png;base64,";
    ifstream png(base64FileName);
    SHASTA_ASSERT(png);
    html << png.rdbuf();
    html << "\"/>";

    // Remove the base64 file.
    std::filesystem::remove(base64FileName);

}



void Assembler::exploreSummary(const vector<string>&, ostream&)
{
    return;
}



void Assembler::writeGnuPlotPngToHtml(
    ostream& html,
    int width,
    int height,
    const string& gnuplotCommands)
{

    // Create a file to contain gnuplot commands.
    const string gnuplotFileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    const string pngFileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    {
        ofstream gnuplotFile(gnuplotFileName);
        gnuplotFile <<
            "set terminal pngcairo size " << width << "," << height <<
            " font 'Noto Serif'\n"
            "set output '" << pngFileName << "'\n" <<
            gnuplotCommands;
    }

    // Invoke gnuplot.
    const string command = "gnuplot " + gnuplotFileName;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " +
            to_string(errorCode) + " " + strerror(errorCode) +
            "\nrunning command: " + command);
    }

    // Write the png file to html.
    writePngToHtml(html, pngFileName);

    // Remove the files we created.
    std::filesystem::remove(gnuplotFileName);
    std::filesystem::remove(pngFileName);
}



// Access all available assembly data, without throwing exceptions
void Assembler::accessAllSoft()
{

    bool allDataAreAvailable = true;

    try {
        accessReadSummaries();
    } catch(const exception& e) {
        cout << "The read summaries not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessKmerChecker();
    } catch(const exception& e) {
        cout << "The k-mer checker is not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkers();
    } catch(const exception& e) {
        cout << "Markers are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerKmers();
    } catch(const exception& e) {
        cout << "Marker k-mers are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAnchors(false);
    } catch(const exception& e) {
        cout << "The anchors are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessJourneys();
    } catch(const exception& e) {
        cout << "The journeys are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAnchorGraph();
    } catch(const exception& e) {
        cout << "The AnchorGraph is not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessSimpleAnchorGraph();
    } catch(const exception& e) {
    }


    if(!allDataAreAvailable) {
        cout << "Not all assembly data are accessible." << endl;
        cout << "Some functionality is not available." << endl;
    }
}




void Assembler::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: monospace;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 1px solid #b8b5c7d9;
        padding: 2px;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.left {
        text-align: left;
    }
    td.right {
        text-align: right;
    }
    td.smaller {
        font-size: smaller;
    }
    a {
        color: DarkSlateBlue;
    }

    /* This can be used to get vertical text in table cells. */
    span.rotated 
    {
      writing-mode: vertical-rl;
      transform: rotate(180deg);
    }

    ul.navigationMenu {
        list-style-type: none;
        margin: 0px 0px 12px 0px;
        padding: 0;
        overflow: hidden;
        background-color: #404040;
    }

    div.navigationButton {
        display: inline-block;
        color: white;
        text-align: center;
        padding: 14px 16px;
        text-decoration: none;
        // min-width: 120px;
    }

    .navigationMenuEntry:hover .navigationButton {
        background-color: black;
    }
    
    li.navigationMenuEntry {
        display: inline-block;
    }
    
    .navigationItems {
        display: none;
        position: absolute;
        background-color: DodgerBlue;
        // min-width: 120px;
        box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
        z-index: 1;
    }
    
    a.navigationItem {
        color: black;
        padding: 12px 16px;
        text-decoration: none;
        display: block;
        text-align: left;
    }
    
    .navigationItems a:hover {background-color: SteelBlue}
    
    .navigationMenuEntry:hover .navigationItems {
        display: block;
    }

    input[type=submit] {
        background-color: #89bef2;
        padding: 4px;
        margin: 2px;
        border-radius: 8px;
    }

    input[type=button] {
        padding: 4px;
    }

    input[type=text], input[type=radio] {
        background-color: #ecf1f0;
        border-width: thin;
    }

    button {
        background-color: #89bef2;
        padding: 4px;
        margin: 2px;
        border-radius: 8px;
    }

</style>
    )%";
}


void Assembler::writeHtmlBegin(ostream& html) const
{
    html <<
        "\r\n"
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<meta charset='UTF-8'>"
        "<title>Shasta assembler</title>";
    writeStyle(html);
    writeMakeAllTablesCopyable(html);
    html <<
        "</head>"
        "<body onload='makeAllTablesCopyable()'>";
}



void Assembler::writeHtmlEnd(ostream& html) const
{
    html << "</body>";
    html << "</html>";
}



void shasta::writeStrandSelection(
    ostream& html,          // The html stream to write the form to.
    const string& name,     // The selection name.
    bool select0,           // Whether strand 0 is selected.
    bool select1)           // Whether strand 1 is selected.
{
    html <<
        "<select name=" << name <<
        " title='Choose 0 (+) for the input read or 1 (-) for its reverse complement'>"

        "<option value=0"
        << (select0 ? " selected" : "") <<
        ">0 (+)</option>"

        "<option value=1"
        << (select1 ? " selected" : "") <<
        ">1 (-)</option>"

        "</select>";

}
