// Shasta.
#include "LocalReadGraph.hpp"
#include "computeLayout.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "ReadGraph.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>

// Standard library.
#include <queue>


LocalReadGraph::LocalReadGraph(
    const ReadGraph& readGraph,
    const vector<string>& request,
    ostream& html) :
    readGraph(readGraph),
    html(html)

{
    LocalReadGraph& graph = *this;

    findOutIfCustomLayoutIsAvailable();
    getRequestOptions(request);

    writeHeader();
    writeForm();

    parseStartOrientedReadIds();
    if(startOrientedReadIds.empty()) {
        return;
    }

    createVertices();
    createEdges();
    html << "<br>The local read graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

    computeLayout();
    computeBoundingBox();
    write();
}



void LocalReadGraph::findOutIfCustomLayoutIsAvailable()
{
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA2_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);

    customLayoutIsAvailable = (returnCode == 0);
}



void LocalReadGraph::getRequestOptions(const vector<string>& request)
{
    HttpServer::getParameterValue(request, "startOrientedReadIdsString", startOrientedReadIdsString);
    HttpServer::getParameterValue(request, "maxDistance", maxDistance);
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);

    string suppressEdgesBetweenVerticesAtMaxDistanceString;
    suppressEdgesBetweenVerticesAtMaxDistance = HttpServer::getParameterValue(request,
        "suppressEdgesBetweenVerticesAtMaxDistance", suppressEdgesBetweenVerticesAtMaxDistanceString);

    string suppressCrossStrandEdgesString;
    suppressCrossStrandEdges = HttpServer::getParameterValue(request,
        "suppressCrossStrandEdges", suppressCrossStrandEdgesString);

    // The default layout is custom, if available, instead of sfdp.
    if(customLayoutIsAvailable) {
        layoutMethod = "custom";
    }
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);
}



void LocalReadGraph::writeHeader()
{
    html << "<h2>Local read graph</h2>";
}




void LocalReadGraph::writeForm()
{
    html <<
        "<form><table>"
        "<tr>"
        "<th class=left>Starting oriented read ids"
        "<td class=centered><input type=text name=startOrientedReadIdsString style='text-align:center' required";
    if(not startOrientedReadIdsString.empty()) {
        html << " value='" << startOrientedReadIdsString + "'";
    }
    html <<
        " size=20>";

    html << "<tr>"
        "<th class=left>Maximum distance"
        "<td class=centered>"
        "<input type=text name=maxDistance style='text-align:center' required size=8 value=" <<
        maxDistance << ">";

    html <<
        "<tr>"
        "<th class=left title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>";

    html <<
        "<tr>"
        "<th class=left>Layout method"
        "<td class=left>"
        "<input type=radio required name=layoutMethod value='sfdp'" <<
        (layoutMethod == "sfdp" ? " checked=on" : "") <<
        ">sfdp"
        "<br><input type=radio required name=layoutMethod value='fdp'" <<
        (layoutMethod == "fdp" ? " checked=on" : "") <<
        ">fdp"
        "<br><input type=radio required name=layoutMethod value='neato'" <<
        (layoutMethod == "neato" ? " checked=on" : "") <<
        ">neato"
        "<br><input type=radio required name=layoutMethod value='dot'" <<
        (layoutMethod == "dot" ? " checked=on" : "") <<
        ">dot";

    // If command "customLayout" is available, add an option for that.
    if(customLayoutIsAvailable) {
        html <<
            "<br><input type=radio required name=layoutMethod value='custom'" <<
            (layoutMethod == "custom" ? " checked=on" : "") <<
            ">custom";
    }

    html <<
        "<tr><th class=left>Suppress edges between vertices at maximum distance"
        "<td class=centered><input type=checkbox name=suppressEdgesBetweenVerticesAtMaxDistance" <<
        (suppressEdgesBetweenVerticesAtMaxDistance ? " checked" : "") <<
        ">";

    html <<
        "<tr><th class=left>Suppress cross-strand edges"
        "<td class=centered><input type=checkbox name=suppressCrossStrandEdges" <<
        (suppressCrossStrandEdges ? " checked" : "") <<
        ">";

    html <<
        "</table>"
        "<br>"
        "<input type=submit value='Create the local read graph'>"
        "</form>";

}



void LocalReadGraph::parseStartOrientedReadIds()
{
    // Loop over tokens in the startVerticesString.
    boost::trim(startOrientedReadIdsString);
    boost::tokenizer< boost::char_separator<char> > tokenizer(startOrientedReadIdsString, boost::char_separator<char>(", "));
    for(const string& startOrientedReadIdString: tokenizer) {
        const OrientedReadId orientedReadId = OrientedReadId(startOrientedReadIdString);

        // Check that it is a valid OrientedReadId for this assembly.
        if(orientedReadId.getReadId() >= readGraph.edgePairs.size()) {
            throw runtime_error(startOrientedReadIdString + " is not a valid oriented read id for this assembly.");
        }

        startOrientedReadIds.push_back(orientedReadId);
    }

    // Sort and remove duplicates.
    deduplicate(startOrientedReadIds);
}



void LocalReadGraph::createVertices()
{
    LocalReadGraph& graph = *this;

    // BFS queue.
    std::queue<vertex_descriptor> q;

    // Create the vertices corresponding to the initial OrientedReadIds and add them to the queue.
    for(const OrientedReadId orientedReadId: startOrientedReadIds) {
        const vertex_descriptor v = add_vertex(LocalReadGraphVertex(orientedReadId, 0), graph);
        vertexMap.insert(make_pair(orientedReadId, v));
        q.push(v);
    }



    // BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalReadGraphVertex& vertex0 = graph[v0];
        const OrientedReadId orientedReadId0 = vertex0.orientedReadId;
        const ReadId readId0 = orientedReadId0.getReadId();
        const ReadId strand0 = orientedReadId0.getStrand();
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        const span<const uint64_t> edgePairIndexes = readGraph.connectivityTable[readId0];
        for(const uint64_t edgePairIndex: edgePairIndexes) {
            const ReadGraph::EdgePair& edgePair = readGraph.edgePairs[edgePairIndex];
            if(suppressCrossStrandEdges and edgePair.isCrossStrand) {
                continue;
            }

            ReadId readId1;
            Strand strand1;
            if(readId0 == edgePair.readId0) {
                readId1 = edgePair.readId1;
                strand1 = (edgePair.isSameStrand ? strand0 : 1 - strand0);
            } else if(readId0 == edgePair.readId1) {
                readId1 = edgePair.readId0;
                strand1 = (edgePair.isSameStrand ? strand0 : 1 - strand0);
            } else {
                cout << "Assertion failing with:" << endl;
                cout << "orientedReadId0 " << orientedReadId0 << endl;
                cout << "EdgePair " << edgePair.readId0 << " " << edgePair.readId1 << " " << int(edgePair.isSameStrand) << endl;
                SHASTA2_ASSERT(0);
            }
            const OrientedReadId orientedReadId1(readId1, strand1);

            const auto it = vertexMap.find(orientedReadId1);
            if(it == vertexMap.end()) {
                const vertex_descriptor v1 =
                    add_vertex(LocalReadGraphVertex(orientedReadId1, distance1), graph);
                vertexMap.insert(make_pair(orientedReadId1, v1));
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }

    }

}



void LocalReadGraph::createEdges()
{

    LocalReadGraph& graph = *this;

    BGL_FORALL_VERTICES(v0, graph, LocalReadGraph) {

        const OrientedReadId orientedReadId0 = graph[v0].orientedReadId;
        const ReadId readId0 = orientedReadId0.getReadId();
        const ReadId strand0 = orientedReadId0.getStrand();

        const span<const uint64_t> edgePairIndexes = readGraph.connectivityTable[readId0];
        for(const uint64_t edgePairIndex: edgePairIndexes) {
            const ReadGraph::EdgePair& edgePair = readGraph.edgePairs[edgePairIndex];
            if(suppressCrossStrandEdges and edgePair.isCrossStrand) {
                continue;
            }

            ReadId readId1;
            Strand strand1;
            if(readId0 == edgePair.readId0) {
                readId1 = edgePair.readId1;
                strand1 = (edgePair.isSameStrand ? strand0 : 1 - strand0);
            } else if(readId0 == edgePair.readId1) {
                readId1 = edgePair.readId0;
                strand1 = (edgePair.isSameStrand ? strand0 : 1 - strand0);
            } else {
                cout << "Assertion failing with:" << endl;
                cout << "orientedReadId0 " << orientedReadId0 << endl;
                cout << "EdgePair " << edgePair.readId0 << " " << edgePair.readId1 << " " << int(edgePair.isSameStrand) << endl;
                SHASTA2_ASSERT(0);
            }
            const OrientedReadId orientedReadId1(readId1, strand1);

            if(orientedReadId0 < orientedReadId1) {

                const auto it1 = vertexMap.find(orientedReadId1);
                if(it1 != vertexMap.end()) {
                    const vertex_descriptor v1 = it1->second;

                    // If requested, don't include edges between vertices at maximum distance.
                    if(suppressEdgesBetweenVerticesAtMaxDistance and
                        (graph[v0].distance == maxDistance) and
                        (graph[v1].distance == maxDistance)) {
                        continue;
                    }

                    add_edge(v0, v1, LocalReadGraphEdge(edgePairIndex), graph);
                }
            }
        }
    }
}



void LocalReadGraph::computeLayout()
{
    LocalReadGraph& graph = *this;

    // Create a map containing the desired length for each edge.
    std::map<edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
        edgeLengthMap.insert(make_pair(e, 1));
    }



    // Compute the graph layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    const double timeout = 60.;
    if(layoutMethod == "custom") {
        const int quality = 2;
        computeLayoutCustom(
            graph,
            edgeLengthMap,
            layout,
            quality,
            timeout);
    } else {
        const string additionalOptions = "";
        computeLayoutGraphviz(
            graph,
            layoutMethod,
            timeout,
            layout,
            additionalOptions,
            &edgeLengthMap);
    }



    // Store the layout in the vertices.
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        LocalReadGraphVertex& vertex = graph[v];
        const array<double, 2>& position = layout[v];
        vertex.x = position[0];
        vertex.y = position[1];
    }
}



void LocalReadGraph::Box::extend(double factor)
{
    const double extend = factor * max(xSize(), ySize());
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
}



void LocalReadGraph::Box::makeSquare()
{
    if(xSize() > ySize()) {
        const double delta = (xSize() - ySize()) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = (ySize() - xSize()) / 2.;
        xMin -= delta;
        xMax += delta;
    }
}



void LocalReadGraph::computeBoundingBox()
{
    LocalReadGraph& graph = *this;

    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = std::numeric_limits<double>::min();
    boundingBox.yMin = std::numeric_limits<double>::max();
    boundingBox.yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        LocalReadGraphVertex& vertex = graph[v];
        const double x = vertex.x;
        const double y = vertex.y;
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }
}






void LocalReadGraph::write() const
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    // Define the viewport box.
    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const string svgId = "LocalReadAnchorGraph";
    html <<
        "\n<br><br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first to avoid obscuring the vertices.
    writeEdges();
    writeVertices();

    // Finish the svg.
    html << "</svg></div>";
    writeSvgControls();
}



void LocalReadGraph::writeSvgControls() const
{
    // Add drag and zoom.
    addSvgDragAndZoom(html);

    // Side panel.
    html << "<div style='display:inline-block;margin-left:20px'>"
        "<p><table>";

    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Vertex size<td>
    <button type='button' onClick='changeVertexSize(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize(10.)' style='width:3em'>+++</button>

    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>

    <tr><th class=left>Zoom<td>
    <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
    <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>

    <script>
        function changeVertexSize(factor)
        {
            var vertexGroup = document.getElementById('vertices');
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            numericWidth = parseFloat(edgeGroup.style.strokeWidth);
            edgeGroup.style.strokeWidth = numericWidth * factor + 'px';
        }
    </script>
    )stringDelimiter";


    // End of side panel.
    html << "</table>";
    html << "</div>";

    // Scroll down to the svg.
    const string svgId = "LocalReadAnchorGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<br>Use Ctrl+Click to pan."
        "<br>Use Ctrl-Wheel or the above buttons to zoom.";
}





void LocalReadGraph::writeVertices() const
{
    const LocalReadGraph& graph = *this;


    html << "\n<g id='vertices' style='stroke:none'>";
    const double vertexRadius = 0.05;
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        const LocalReadGraphVertex& vertex = graph[v];
        const string vertexColor = (vertex.distance == maxDistance ? "Cyan" : "Green");

        html << "\n<circle cx='" << vertex.x << "' cy='" << vertex.y <<
            "' fill='" << vertexColor <<
            "' r='" << vertexRadius << "' >"
            "<title>" <<
            vertex.orientedReadId <<
            "</title>"
            "</circle>";
    }
    html << "\n</g>";

}



void LocalReadGraph::writeEdges() const
{
    const LocalReadGraph& graph = *this;

    html << "\n<g id='edges' style='stroke:Black;stroke-width:0.005'>";

    BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const LocalReadGraphVertex& vertex0 = graph[v0];
        const LocalReadGraphVertex& vertex1 = graph[v1];

        html <<
            "\n<line "
            "x1='" << vertex0.x << "' y1='" << vertex0.y <<
            "' x2='" << vertex1.x << "' y2='" << vertex1.y << "'";

        const uint64_t edgePairIndex = graph[e].edgePairIndex;
        if(readGraph.edgePairs[edgePairIndex].isCrossStrand) {
            html << " stroke=\"red\"";
        }

        html << " />";
    }

    html << "\n</g>";
}
