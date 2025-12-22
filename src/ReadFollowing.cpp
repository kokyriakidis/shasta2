// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing.hpp"
#include "findLinearChains.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "Options.hpp"
using namespace shasta2;
using namespace ReadFollowing;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>



Graph::Graph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
	const bool debug = true;

    // Initial creation with all possible vertices and edges.
    createVertices();
    createEdges();
    computeEdgeScores();
    if(debug) {
    	write("A");
    }

    // Remove edges with negative offsets.
    removeNegativeOffsetEdges();
    if(debug) {
    	write("B");
    }

    // Remove edges with low commonCount.
    removeLowCommonCountEdges(assemblyGraph.options.readFollowingMinCommonCount);
    if(debug) {
    	write("C");
    }

    // Remove edges with low correctedJaccard.
    removeLowCommonCorrectedJaccardEdges(assemblyGraph.options.readFollowingMinCorrectedJaccard);
    if(debug) {
    	write("D");
    }

    // Prune short leaves.
    prune();
    if(debug) {
    	write("E");
    }

    // Remove edges that don't have the best score at least one direction.
    removeNonBestScoreEdges();
    if(debug) {
    	write("F");
    }
}



// Create vertices of the ReadFollowing graph.
// Each vertex corresponds to a Segment of the AssemblyGraph.
void Graph::createVertices()
{
    // EXPOSE WHEN CODE STABILIZES.
    const double maxCoverage = 18.;

    Graph& graph = *this;

    // Each Segment generates a Vertex.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const Vertex vertex(assemblyGraph, segment);
        if(vertex.coverage < maxCoverage) {
        const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
        vertexMap.insert(make_pair(segment, v));
        }
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

    coverage = edge.lengthWeightedAverageCoverage();

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

    SHASTA2_ASSERT(segmentPairInformation.commonCount > 0);
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



// Remove edges that are not flagged as "best edge" in at least one direction.
void Graph::removeNonBestScoreEdges()
{
    Graph& graph = *this;

    // Make sure the best edge flags are valid.
    computeEdgeScores();
    setBestEdgeFlags();

    // Gather the edges to be removed.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];
        if(not (edge.isBest0 or edge.isBest1)) {
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
    dot << "digraph ReadFollowing {\n";
    dot << std::fixed << std::setprecision(1);

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id;

        // Begin attributes.
        dot << " [";

        // Label.
        dot <<
            "label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            vertex.coverage << "x" <<
            "\"";

        // Color.
        if(vertex.length >= assemblyGraph.options.readFollowingSegmentLengthThreshold) {
            dot << " style=filled fillcolor=cyan";
        }

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    dot << std::fixed << std::setprecision(2);
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

        // Thickness is proportional to commonCount.
        dot << " penwidth=" << 0.2 * double(edge.segmentPairInformation.commonCount);

        // Color is determined by correctedJaccard for the edge.
        // Green = 1
        // Red = assemblyGraph.options.readFollowingMinCorrectedJaccard.
        double hue;
        if(edge.segmentPairInformation.correctedJaccard >= 1.) {
            hue = 1.;
        } else if(edge.segmentPairInformation.correctedJaccard <= assemblyGraph.options.readFollowingMinCorrectedJaccard) {
            hue = 0.;
        } else {
            hue =
                (edge.segmentPairInformation.correctedJaccard - assemblyGraph.options.readFollowingMinCorrectedJaccard) /
                (1. - assemblyGraph.options.readFollowingMinCorrectedJaccard);
        }
        hue /= 3.;
        dot << std::fixed << std::setprecision(3) << " color=\""  << hue << " 1. 1.\"";

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



// Find a best path starting at the given vertex and
// ending if one of the terminalVertices is encountered.
// Direction is 0 for forward and 1 backward.
void Graph::findPath(
    vertex_descriptor v, uint64_t direction,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices) const
{
    if(direction == 0) {
        findForwardPath(v, path, stopVertices);
    } else {
        findBackwardPath(v, path, stopVertices);
    }
}



void Graph::findForwardPath(
    vertex_descriptor v,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices) const
{
    const Graph& graph = *this;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);



    // At each iteration, add one vertex to the path.
    // Use the edge with best score
    while(out_degree(v, assemblyGraph) > 0) {

        // Find the edge with best score.
        double bestScore = 0.;
        edge_descriptor eBest = edge_descriptor({0, 0, 0});
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const double score = graph[e].score;
            if(score > bestScore) {
                bestScore = score;
                eBest = e;
            }
        }

        // Add to the path the target of this vertex and continue from here.
        v = target(eBest, graph);
        path.push_back(v);

        if(stopVertices.contains(v)) {
            break;
        }
    }
}



void Graph::findBackwardPath(
    vertex_descriptor v,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices) const
{
    const Graph& graph = *this;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);



    // At each iteration, add one vertex to the path.
    // Use the edge with minimum offset.
    while(in_degree(v, assemblyGraph) > 0) {

        // Find the edge with best score.
        double bestScore = 0.;
        edge_descriptor eBest  = edge_descriptor({0, 0, 0});
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            const double score = graph[e].score;
            if(score > bestScore) {
                bestScore = score;
                eBest = e;
            }
        }

        // Add to the path the source of this vertex and continue from here.
        v = source(eBest, graph);
        path.push_back(v);

        if(stopVertices.contains(v)) {
            break;
        }
    }

    // Reverse the path so it goes forward.
    std::ranges::reverse(path);
}



void Graph::writePath(Segment segment, uint64_t direction) const
{
    const Graph& graph = *this;

    const auto it = vertexMap.find(segment);
    SHASTA2_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;

    vector<vertex_descriptor> path;
    std::set<vertex_descriptor> stopVertices;
    findPath(v, direction, path, stopVertices);

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
        if(graph[v].length < assemblyGraph.options.readFollowingPruneLength) {
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
// - isBest0 is set if this edge has the best score among all out-edges of v0.
// - isBest1 is set if this edge has the best score among all in-edges of v1.
void Graph::setBestEdgeFlags()
{
    Graph& graph = *this;

    // First set all the flags to false;
    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        edge.isBest0 = false;
        edge.isBest1 = false;
    }



    // Then loop over all vertices to set the flags.
    BGL_FORALL_VERTICES(v, graph, Graph) {

        // Set the isBest0 flag for the out-edge with the best score.
        if(out_degree(v, graph) > 0) {
            edge_descriptor eBest;
            double bestScore = 0.;
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                const double score = graph[e].score;
                if(score > bestScore) {
                    bestScore = score;
                    eBest = e;
                }
            }
            graph[eBest].isBest0 = true;
        }

        // Set the isBest1 flag for the in-edge with the best score.
        if(in_degree(v, graph) > 0) {
            edge_descriptor eBest;
            double bestScore = 0.;
            BGL_FORALL_INEDGES(v, e, graph, Graph) {
                const double score = graph[e].score;
                if(score > bestScore) {
                    bestScore = score;
                    eBest = e;
                }
            }
            graph[eBest].isBest1 = true;
        }

    }
}



// Find assembly paths.
// These are paths between vertices corresponding to long segments.
// Note these are paths in the ReadFollowing::Graph but not in the AssemblyGraph.
void Graph::findPaths(vector< vector<Segment> >& assemblyPaths) const
{
    const Graph& graph = *this;
    const bool debug = true;


    // A graph to store the paths we find.
    // Each vertex corresponds to a long segment.
    // An edge u0->u1 contains a path that starts at segment(u0)
    // and ends at segment(u1). We only keep one path u0->u1,
    // even though in many/most cases we will find two,
    // one in each direction.
    class PathGraphVertex {
    public:
        Segment segment;
    };
    class PathGraphEdge {
    public:
        vector<Segment> path;
    };
    using PathGraph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        PathGraphVertex,
        PathGraphEdge>;
    PathGraph pathGraph;
    std::map<Segment, PathGraph::vertex_descriptor> pathGraphVertexMap;

    // Each long Segment generates a PathGraphVertex.
    std::set<vertex_descriptor> longSegments;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        if(vertex.length >= assemblyGraph.options.readFollowingSegmentLengthThreshold) {
            longSegments.insert(v);
            const PathGraph::vertex_descriptor u = boost::add_vertex({vertex.segment}, pathGraph);
            pathGraphVertexMap.insert({vertex.segment, u});
        }
    }



    // For each PathGraphVertex, compute a path in each direction,
    // always stopping when another long segment is encountered.
    // Each path generates a PathGraphEdge, as long as an edge between
    // the same two PathGraph vertices does not already exist.
    vector<vertex_descriptor> path;
    BGL_FORALL_VERTICES(u, pathGraph, PathGraph) {
        const Segment segment = pathGraph[u].segment;
        const auto it = vertexMap.find(segment);
        SHASTA2_ASSERT(it != vertexMap.end());
        const vertex_descriptor v = it->second;
        for(uint64_t direction=0; direction<2; direction++) {
            findPath(v, direction, path, longSegments);

            if(path.size() > 1) {

                // Find the first/last segment of the path.
                const vertex_descriptor v0 = path.front();
                const vertex_descriptor v1 = path.back();
                const Segment segment0 = graph[v0].segment;
                const Segment segment1 = graph[v1].segment;

                // Locate the corresponding PathGraph vertices.
                const auto it0 = pathGraphVertexMap.find(segment0);
                const auto it1 = pathGraphVertexMap.find(segment1);
                if((it0 != pathGraphVertexMap.end()) and (it1 != pathGraphVertexMap.end())) {
                    const PathGraph::vertex_descriptor u0 = it0->second;
                    const PathGraph::vertex_descriptor u1 = it1->second;

                    // Create an edge, as long as an edge between
                    // the same two PathGraph vertices does not already exist.
                    bool edgeExists = false;
                    tie(ignore, edgeExists) = boost::edge(u0, u1, pathGraph);
                    if(not edgeExists) {

                        // Ok, we are going to add a PathGraph edge.

                        edge_descriptor e;
                        tie(e, ignore) = boost::add_edge(u0, u1, pathGraph);

                        // Fill in the path of this PathGraphEdge.
                        PathGraphEdge& pathGraphEdge = pathGraph[e];
                        for(const vertex_descriptor v: path) {
                            pathGraphEdge.path.push_back(graph[v].segment);
                        }
                    }
                }
            }
        }
    }



    if(debug) {
    	ofstream dot("PathGraph.dot");
        dot << "digraph PathGraph {\n";

        BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        	const Segment segment = pathGraph[v].segment;
        	dot << assemblyGraph[segment].id << ";\n";
        }

        BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        	const PathGraph::vertex_descriptor v0 = source(e, pathGraph);
        	const PathGraph::vertex_descriptor v1 = target(e, pathGraph);
        	const Segment segment0 = pathGraph[v0].segment;
        	const Segment segment1 = pathGraph[v1].segment;
        	dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id << ";\n";
        }

    	dot << "}\n";
    }



    // Find linear chains of vertices in the PathGraph.
    vector< vector<PathGraph::vertex_descriptor> > chains;
    findLinearVertexChains(pathGraph, chains);



    // Each linear vertex chain generates a Segment sequence,
    // that is, a sequence of Segments that
    // should be assembled into a single Segment.
    assemblyPaths.clear();
    for(const vector<PathGraph::vertex_descriptor>& chain: chains) {
    	if(chain.size() < 2) {
    		continue;
    	}
        const Segment segment0 = pathGraph[chain.front()].segment;
        const Segment segment1 = pathGraph[chain.back()].segment;

        if(debug) {
			cout << "Found an assembly path that begins at " <<
				assemblyGraph[segment0].id << " and ends at " <<
				assemblyGraph[segment1].id << endl;
        }

        // Create a new assembly path.
        assemblyPaths.emplace_back();
        vector<Segment>& assemblyPath = assemblyPaths.back();
        for(uint64_t i1=1; i1<chain.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const PathGraph::vertex_descriptor u0 = chain[i0];
            const PathGraph::vertex_descriptor u1 = chain[i1];

            // Locate the PathGraph edge.
            PathGraph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(u0, u1, pathGraph);
            SHASTA2_ASSERT(edgeExists);
            const PathGraphEdge& pathGraphEdge = pathGraph[e];
            const vector<Segment>& path = pathGraphEdge.path;

            SHASTA2_ASSERT(pathGraph[u0].segment == path.front());
            SHASTA2_ASSERT(pathGraph[u1].segment == path.back());

            // The path for this vertex already contains the segments
            // corresponding to u0 and u1. So, to avoid duplications,
            // for each edge except the last we copy the path
            // without the last segment.
            // For the last edge we copy the entire path.
            auto end = path.end();
            if(i1 != chain.size() - 1) {
                --end;
            }
            copy(path.begin(), end, back_inserter(assemblyPath));
        }
    }
}



void Graph::writePaths() const
{
    vector< vector<Segment> > assemblyPaths;
    findPaths(assemblyPaths);

    ofstream csv("AssemblyPaths.csv");
    cout << "Found " << assemblyPaths.size() << " assembly paths." << endl;
    for(const vector<Segment>& assemblyPath: assemblyPaths) {
        cout << "Assembly path with " << assemblyPath.size() <<
        " segments beginning at " << assemblyGraph[assemblyPath.front()].id <<
        " and ending at " << assemblyGraph[assemblyPath.back()].id << endl;

        for(const Segment& segment: assemblyPath) {
            csv << assemblyGraph[segment].id << ",";
        }
        csv << "\n";
    }
}



// This will require some experimentation.
void Graph::computeEdgeScores()
{
    Graph& graph = *this;

    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        edge.score = double(edge.segmentPairInformation.commonCount) * edge.segmentPairInformation.correctedJaccard;
    }
}
