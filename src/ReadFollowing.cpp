// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "Options.hpp"
using namespace shasta;
using namespace ReadFollowing;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>



Graph::Graph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    // Initial creation with all possible vertices and edges.
    createVertices();
    createEdges();
    write("A");

    // Remove edges with negative offsets.
    removeNegativeOffsetEdges();
    write("B");

    // Remove edges with low commonCount.
    removeLowCommonCountEdges(minCommonCount);
    write("C");

    // Remove edges with low correctedJaccard.
    removeLowCommonCorrectedJaccardEdges(minCorrectedJaccard);
    write("D");

    // Prune short leaves.
    prune();
    write("E");

    // Remove edges that have both isLowestOffset0 and isLowestOffset1 set to false.
    removeNonLowestOffsetEdges();
    write("F");
}



// Create vertices of the ReadFollowing graph.
// Each vertex corresponds to a Segment of the AssemblyGraph.
void Graph::createVertices()
{
    Graph& graph = *this;

    // Each Segment generates a Vertex.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
        vertexMap.insert(make_pair(segment, v));
    }

}



Vertex::Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment) :
    segment(segment)
{
    const AssemblyGraphEdge& edge = assemblyGraph[segment];
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }

    // Compute initial/final support.
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, initialSupport);
    SegmentStepSupport::getFinalLast   (assemblyGraph, segment, representativeRegionStepCount, finalSupport  );
}



// Create edges of the ReadFollowing graph.
// An edge v0->v1 is created if the final support of v0
// shares at least one OrientedReadId with the initial support of v1.
void Graph::createEdges()
{
    Graph& graph = *this;
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();

    // For each OrientedReadId, gather the vertices that the OrientedReadId
    // appears in, in the initial/final support.
    vector< vector<vertex_descriptor> > initialSupportVertices(orientedReadCount);
    vector< vector<vertex_descriptor> > finalSupportVertices(orientedReadCount);
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];

        for(const SegmentStepSupport& s: vertex.initialSupport) {
            initialSupportVertices[s.orientedReadId.getValue()].push_back(v);
        }

        for(const SegmentStepSupport& s: vertex.finalSupport) {
            finalSupportVertices[s.orientedReadId.getValue()].push_back(v);
        }
    }



    // An edge v0->v1 will be created if the final support of v0
    // shares at least one OrientedReadId with the initial support of v1.
    std::set< pair<vertex_descriptor, vertex_descriptor> > vertexPairs;
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const vector<vertex_descriptor>& initialVertices = initialSupportVertices[i];
        const vector<vertex_descriptor>& finalVertices = finalSupportVertices[i];
        for(const vertex_descriptor v0: finalVertices) {
            for(const vertex_descriptor v1: initialVertices) {
                if(v1 == v0) {
                    continue;
                }

                // This OrientedReadId appears in the final support of v0 and in the
                // initial support of v1, so we will create an edge v0->v1.
                vertexPairs.insert({v0, v1});
            }
        }
    }


    // Generate an edge for each of these pairs.
    for(const auto& p: vertexPairs) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;

        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        add_edge(v0, v1, Edge(assemblyGraph, segment0, segment1), graph);
    }
}



Edge::Edge(
    const AssemblyGraph& assemblyGraph,
    Segment segment0,
    Segment segment1)
{
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    ostream html(0);
    segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
        html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

    SHASTA_ASSERT(segmentPairInformation.commonCount > 0);
}



void Graph::removeLowCommonCountEdges(uint64_t minCommonCount)
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].segmentPairInformation.commonCount < minCommonCount) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



void Graph::removeNegativeOffsetEdges()
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].segmentPairInformation.segmentOffset < 0) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



// Remove edges that have both isLowestOffset0 and isLowestOffset1 set to false.
void Graph::removeNonLowestOffsetEdges()
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];
        if(not (edge.isLowestOffset0 or edge.isLowestOffset1)) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



void Graph::removeLowCommonCorrectedJaccardEdges(double minCorrectedJaccard)
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].segmentPairInformation.correctedJaccard < minCorrectedJaccard) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



void Graph::write(const string& name)
{
    cout << "ReadFollowing-" << name << ": " << num_vertices(*this) <<
        " vertices, " << num_edges(*this) << " edges." << endl;
    writeGraphviz(name);
    writeCsv(name);
}



void Graph::writeGraphviz(const string& name)
{
    const Graph& graph = *this;

    // Make sure the lowest offset flags are valid.
    setLowestOffsetFlags();

    ofstream dot("ReadFollowing-" + name + ".dot");
    dot << "digraph ReadFollowing1 {\n";
    dot << std::fixed << std::setprecision(2);

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id <<
            " [label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            vertex.initialSupport.size() << "/" <<
            vertex.finalSupport.size() <<
            "\"";
        dot <<
            "]"
            ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];
        const int32_t offset = edge.segmentPairInformation.segmentOffset;

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

        // Begin attributes.
        dot << "[";

        // Label.
        dot << "label=\"" <<
            edge.segmentPairInformation.commonCount << "/" <<
            std::fixed << std::setprecision(2) <<
            edge.segmentPairInformation.correctedJaccard << "\\n" <<
            offset << "\"";

        // Thickness.
        dot << " penwidth=" << 0.2 * double(edge.segmentPairInformation.commonCount);

        // Color.
        string color;
        if(edge.isLowestOffset0) {
            if (edge.isLowestOffset1) {
                color = "Green";
            } else {
                color = "Cyan";
            }
        } else {
            if (edge.isLowestOffset1) {
                color = "Magenta";
            } else {
                color = "Black";
            }
        }
        dot << " color=" << color;

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";

    }

    dot << "}\n";
}



void Graph::writeCsv(const string& name) const
{
    writeVerticesCsv(name);
    writeEdgesCsv(name);
}



void Graph::writeVerticesCsv(const string& name) const
{
    const Graph& graph = *this;

    ofstream csv("ReadFollowing-Vertices-" + name + ".csv");
    csv << "Segment,Length,InitialSupport,FinalSupport,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];

        csv << assemblyGraphEdge.id << ",";
        csv << vertex.length << ",";
        csv << vertex.initialSupport.size() << ",";
        csv << vertex.finalSupport.size() << ",";
        csv << "\n";
    }
}



void Graph::writeEdgesCsv(const string& name) const
{
    const Graph& graph = *this;

    ofstream csv("ReadFollowing-Edges-" + name + ".csv");
    csv << "Segment0,Segment1,Length0,Length1,FinalSupport0,InitialSupport1,"
        "Common,Missing0,Missing1,MissingTotal,CorrectedJaccard,Offset,\n";
    csv << std::fixed << std::setprecision(2);

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];

        const Segment segment0 = vertex0.segment;
        const Segment segment1 = vertex1.segment;

        csv << assemblyGraph[segment0].id << ",";
        csv << assemblyGraph[segment1].id << ",";
        csv << vertex0.length << ",";
        csv << vertex1.length << ",";
        csv << vertex0.finalSupport.size() << ",";
        csv << vertex1.initialSupport.size() << ",";
        csv << edge.segmentPairInformation.commonCount << ",";
        csv << edge.segmentPairInformation.missing0 << ",";
        csv << edge.segmentPairInformation.missing1 << ",";
        csv << edge.segmentPairInformation.missing0 + edge.segmentPairInformation.missing1 << ",";
        csv << edge.segmentPairInformation.correctedJaccard << ",";
        csv << edge.segmentPairInformation.segmentOffset << ",";
        csv << "\n";
    }
}



// Find a minimum offset path starting at the given vertex and
// ending if one of the terminalVertices is encountered.
// Direction is 0 for forward and 1 backward.
void Graph::findPath(
    Segment segment, uint64_t direction,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& terminalVertices) const
{
    if(direction == 0) {
        findForwardPath(segment, path, terminalVertices);
    } else {
        findBackwardPath(segment, path, terminalVertices);
    }
}



void Graph::findForwardPath(
    Segment segment,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& terminalVertices) const
{
    const Graph& graph = *this;

    // Locate the vertex  corresponding to this segment.
    const auto it = vertexMap.find(segment);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);



    // At each iteration, add one vertex to the path.
    // Use the edge with minimum offset.
    while((not terminalVertices.contains(v)) and (out_degree(v, assemblyGraph) > 0)) {

        // Find the edge with lowest offset.
        int32_t lowestOffset = std::numeric_limits<int32_t>::max();
        edge_descriptor eLowestOffset;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const int32_t offset = graph[e].segmentPairInformation.segmentOffset;
            if(offset < lowestOffset) {
                lowestOffset = offset;
                eLowestOffset = e;
            }
        }
        SHASTA_ASSERT(lowestOffset != std::numeric_limits<int32_t>::max());

        // Add to the path the target of this vertex and continue from here.
        v = target(eLowestOffset, graph);
        path.push_back(v);
    }
}



void Graph::findBackwardPath(
    Segment segment,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& terminalVertices) const
{
    const Graph& graph = *this;

    // Locate the vertex  corresponding to this segment.
    const auto it = vertexMap.find(segment);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);



    // At each iteration, add one vertex to the path.
    // Use the edge with minimum offset.
    while((not terminalVertices.contains(v)) and (in_degree(v, assemblyGraph) > 0)) {

        // Find the edge with lowest offset.
        int32_t lowestOffset = std::numeric_limits<int32_t>::max();
        edge_descriptor eLowestOffset;
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            const int32_t offset = graph[e].segmentPairInformation.segmentOffset;
            if(offset < lowestOffset) {
                lowestOffset = offset;
                eLowestOffset = e;
            }
        }
        SHASTA_ASSERT(lowestOffset != std::numeric_limits<int32_t>::max());

        // Add to the path the source of this vertex and continue from here.
        v = source(eLowestOffset, graph);
        path.push_back(v);
    }

    // Reverse the path so it goes forward.
    std::ranges::reverse(path);
}



void Graph::writePath(Segment segment, uint64_t direction) const
{
    const Graph& graph = *this;

    vector<vertex_descriptor> path;
    std::set<vertex_descriptor> terminalVertices;
    findPath(segment, direction, path, terminalVertices);

    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;
}




void Graph::prune()
{
    while(pruneIteration());
}



bool Graph::pruneIteration()
{
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(graph[v].length < pruneLength) {
            const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
            if(isLeaf) {
                verticesToBeRemoved.push_back(v);
            }
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

    return not verticesToBeRemoved.empty();
}



// For each edge v0->v1:
// - isLowestOffset0 is set if this edge has the lowest offset out of all out-edges of v0.
// - isLowestOffset1 is set if this edge has the lowest offset out of all in-edges of v1.
void Graph::setLowestOffsetFlags()
{
    Graph& graph = *this;

    // First set all the flags to false;
    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        edge.isLowestOffset0 = false;
        edge.isLowestOffset1 = false;
    }



    // Then loop over all vertices to set the flags.
    BGL_FORALL_VERTICES(v, graph, Graph) {

        // Set the isLowestOffset0 flag for out-edge with the lowest offset.
        edge_descriptor eLowest;
        int32_t lowestOffset = std::numeric_limits<int32_t>::max();
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const int32_t offset = graph[e].segmentPairInformation.segmentOffset;
            if(offset < lowestOffset) {
                lowestOffset = offset;
                eLowest = e;
            }
        }
        if(lowestOffset != std::numeric_limits<int32_t>::max()) {
            graph[eLowest].isLowestOffset0 = true;
        }

        // Set the isLowestOffset1 flag for in-edge with the lowest offset.
        lowestOffset = std::numeric_limits<int32_t>::max();
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            const int32_t offset = graph[e].segmentPairInformation.segmentOffset;
            if(offset < lowestOffset) {
                lowestOffset = offset;
                eLowest = e;
            }
        }
        if(lowestOffset != std::numeric_limits<int32_t>::max()) {
            graph[eLowest].isLowestOffset1 = true;
        }
    }
}
