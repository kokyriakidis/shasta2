#include "graphvizToHtml.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;

#include <filesystem>
#include "fstream.hpp"
#include <stdexcept.hpp>
#include "string.hpp"



// This  can be used to display a graph in Graphviz format to html.
// It does the following:
// - It constructs a Graphviz command to render the graph in svg format.
// - It invokes the Graphviz command with the given timeout, in seconds.
// - It writes the svg output to the given html ostream.
// - Finally, it removes the Graphviz file and the svg file.
// The options argument should not include the options to specify the
// input and output file name and the svg format. They should only include
// rendering options.
// If successful, thi removes the input Graphviz file and the svg file.
// Otherwise, it throws an exception and does not remove the input graphviz file.

void shasta::graphvizToHtml(

    // The name of the Graphviz file (dot format) containing the graph to be displayed.
    const string& dotFileName,

    // The Graphviz layout to be used (dot, sfdp, etc.).
    const string& layoutName,

    // The timeout in seconds.
    double timeout,

    // Graphviz options. This should not include the options to specify the
    // input and output file name and the svg format. They should only include
    // rendering options.
    const string& options,

    // Svg output goes here.
    shasta::ostream& html
    )
{
    // Construct the Graphviz command.
    const string svgFileName = dotFileName + ".svg";
    const string command = layoutName + " -T svg " + dotFileName + " -o " + svgFileName + " " + options;

    // Run the command with the given timeout.
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    if(signalOccurred) {
        throw runtime_error("Error during graph layout. Command was " + command);
    }
    if(timeoutTriggered) {
        throw runtime_error("Timeout during graph layout. Command was " + command);
    }
    if(returnCode!=0 ) {
        throw runtime_error("Error during graph layout. Command was " + command);
    }

    // Success, we can remove the Graphviz file.
    std::filesystem::remove(dotFileName);

    // Write the svg to html.
    html << "<div style='display:inline-block'>";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();
    html << "</div>";

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}
