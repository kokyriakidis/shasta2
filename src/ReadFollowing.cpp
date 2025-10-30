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
    const uint64_t minCommonCount = 6;
    removeLowCommonCountEdges(minCommonCount);
    write("C");

    // Remove edges with low correctedJaccard.
    const double minCorrectedJaccard = 0.7;
    removeLowCommonCorrectedJaccardEdges(minCorrectedJaccard);
    write("D");

    // Prune short leaves.
    const uint64_t minimumLength = 100000;
    prune(minimumLength);
    write("E");
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



void Graph::write(const string& name) const
{
    cout << "ReadFollowing-" << name << ": " << num_vertices(*this) <<
        " vertices, " << num_edges(*this) << " edges." << endl;
    writeGraphviz(name);
    writeCsv(name);
}



void Graph::writeGraphviz(const string& name) const
{
    const Graph& graph = *this;

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
            edge.segmentPairInformation.segmentOffset << "\"";

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



#if 0
void Graph::findPath(Segment segment, uint64_t direction, vector<vertex_descriptor>& path) const
{
    if(direction == 0) {
        findForwardPath(segment, path);
    } else if(direction == 1) {
        findBackwardPath(segment, path);
    } else {
        SHASTA_ASSERT(0);
    }
}



void Graph::findForwardPath(Segment segment, vector<vertex_descriptor>& path) const
{
    const Graph& graph = *this;

    // Find the start vertex.
    const auto it = vertexMap.find(segment);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Each iteration adds one vertex to the path.
    path.clear();
    path.push_back(v);
    std::set<vertex_descriptor> pathVertices;
    pathVertices.insert(v);
    while(true) {

         // Find the best next vertex.
         vertex_descriptor vNext = null_vertex();
         double bestJaccardSum = 0.;
         BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
             vertex_descriptor v1 = target(e, graph);
             if(pathVertices.contains(v1)) {
                 continue;
             }

             // Compute the sum of Jaccard similarities between all
             // previous vertices in the path and v1.
             double jaccardSum = 0.;
             for(const vertex_descriptor vPrevious: path) {
                 edge_descriptor ePrevious;
                 bool edgeWasFound;
                 tie(ePrevious, edgeWasFound) = edge(vPrevious, v1, graph);
                 if(edgeWasFound) {
                    jaccardSum += graph[ePrevious].jaccard;
                 }
             }

             if(jaccardSum > bestJaccardSum) {
                 vNext = v1;
                 bestJaccardSum = jaccardSum;
             }

         }

         if(vNext == null_vertex()) {
             break;
         }

         v = vNext;

         path.push_back(v);
         pathVertices.insert(v);
     }

}



void Graph::findBackwardPath(Segment, vector<vertex_descriptor>& /* path */) const
{
    SHASTA_ASSERT(0);
}
#endif



void Graph::writePath(Segment /* segment */, uint64_t /* direction */) const
{
    cout << "Not implemented." << endl;

#if 0
    const Graph& graph = *this;

    vector<vertex_descriptor> path;
    findPath(segment, direction, path);

    if(direction == 1) {
        std::ranges::reverse(path);
    }

    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << " ";
    }
    cout << endl;
#endif
}




void Graph::prune(uint64_t minimumLength)
{
    while(pruneIteration(minimumLength));
}



bool Graph::pruneIteration(uint64_t minimumLength)
{
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(graph[v].length < minimumLength) {
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
