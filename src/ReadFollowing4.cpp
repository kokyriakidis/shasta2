// Shasta.
#include "ReadFollowing4.hpp"
#include "Journeys.hpp"
#include "Options.hpp"
using namespace shasta2;
using namespace ReadFollowing4;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "fstream.hpp"
#include <queue>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<ReadFollower>;



ReadFollower::ReadFollower(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<ReadFollower>(*this),
    assemblyGraph(assemblyGraph)
{
    fillSupportMaps();
    findSegmentPairs();
    createVertices();
    createEdges();

    for(uint64_t direction=0; direction<2; direction++) {
        cout << "The initial read following graph for direction " << direction <<
            " has " << num_vertices(graphs[direction]) <<
            " vertices and " << num_edges(graphs[direction]) << " edges." << endl;
    }

    // Prune.
    for(uint64_t direction=0; direction<2; direction++) {
        graphs[direction].prune();
        graphs[direction].writeGraphviz(assemblyGraph, "Pruned-" + to_string(direction));
    }

    for(uint64_t direction=0; direction<2; direction++) {
        cout << "After pruning, the read following graph for direction " << direction <<
            " has " << num_vertices(graphs[direction]) <<
            " vertices and " << num_edges(graphs[direction]) << " edges." << endl;
    }

    // Before we can compute shortest paths we have to create the vertex index map
    // for each of the two graphs.
    for(uint64_t direction=0; direction<2; direction++) {
        graphs[direction].createVertexIndexMap();
    }
}



void ReadFollower::fillSupportMaps()
{
    const uint32_t representativeRegionStepCount = uint32_t(assemblyGraph.options.representativeRegionStepCount);

    vector<SegmentStepSupport> support;
    vector<SegmentStepSupport> finalSupport;

    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, support);
        initialSupportMap.insert(make_pair(segment, support));

        SegmentStepSupport::getFinalLast   (assemblyGraph, segment, representativeRegionStepCount, support);
        finalSupportMap.insert(make_pair(segment, support));
    }
}



void ReadFollower::findSegmentPairs()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;

    // For each OrientedReadId, gather the Segments that the OrientedReadId
    // appears in, in the initial support.
    vector< vector<Segment> > initialSupportSegments(orientedReadCount);
    for(const auto&[segment, support]: initialSupportMap) {
        for(const SegmentStepSupport& s: support) {
            initialSupportSegments[s.orientedReadId.getValue()].push_back(segment);
        }

    }

    // For each OrientedReadId, gather the Segments that the OrientedReadId
    // appears in, in the final support.
    vector< vector<Segment> > finalSupportSegments(orientedReadCount);
    for(const auto&[segment, support]: finalSupportMap) {
        for(const SegmentStepSupport& s: support) {
            finalSupportSegments[s.orientedReadId.getValue()].push_back(segment);
        }

    }



    // Store all segment pairs (segment0, segment1) such that the final support
    // of segment0 shares one or more OrientedReadIds with the initial support of segment1.
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const vector<Segment>& initialSegments = initialSupportSegments[i];
        const vector<Segment>& finalSegments = finalSupportSegments[i];
        for(const Segment segment0: finalSegments) {
            for(const Segment segment1: initialSegments) {
                if(segment0 != segment1) {
                    segmentPairs.push_back(make_pair(segment0, segment1));
                }
            }
        }
    }



    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(segmentPairs, count, minCommonCount);
    SHASTA2_ASSERT(segmentPairs.size() == count.size());

    cout << "Found " << segmentPairs.size() << " segment pairs." << endl;
}



// Each edge of the AssemblyGraph (Segment) generates a vertex in
// each of our two Graphs.
void ReadFollower::createVertices()
{
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        for(uint64_t direction=0; direction<2; direction++) {
            graphs[direction].createVertex(assemblyGraph, segment);
        }
    }
}



void Graph::createVertex(const AssemblyGraph& assemblyGraph, Segment segment)
{
    SHASTA2_ASSERT(not vertexMap.contains(segment));
    Graph& graph = *this;
    const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
    vertexMap.insert(make_pair(segment, v));
}



Vertex::Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment) :
    segment(segment)
{
    // Get the Segment corresponding to this Vertex.
    const AssemblyGraphEdge& edge = assemblyGraph[segment];

    // Get the length.
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }

    // Store the isLong flag
    isLong = (length >= assemblyGraph.options.readFollowingSegmentLengthThreshold);
}



void ReadFollower::createEdges()
{
    uint64_t threadCount = assemblyGraph.options.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    setupLoadBalancing(segmentPairs.size(), 100);
    runThreads(&ReadFollower::createEdgesThreadFunction, threadCount);

}



void ReadFollower::createEdgesThreadFunction([[maybe_unused]] uint64_t threadId)
{
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);

    ostream html(0);

    // Prepare vectors of edges to be added to the two graphs.
    // We will add them all at the end so we only have to acquire the mutex once.
    class EdgeToBeAdded {
    public:
        Graph::vertex_descriptor v0;
        Graph::vertex_descriptor v1;
        Edge edge;
    };
    array<vector<EdgeToBeAdded>, 2> edgesToBeAdded;

    // Loop over batches of segment pairs assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segment pairs in this batch.
        for(uint64_t i=begin; i<end; i++) {
            const auto& [segment0, segment1] = segmentPairs[i];

            const SegmentPairInformation segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
                html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

            // If it does not satisfy our requirements, skip it.
            if(segmentPairInformation.commonCount < minCommonCount) {
                continue;
            }
            if(segmentPairInformation.segmentOffset < 0) {
                continue;
            }
            if(not assemblyGraph.canConnect(segment0, segment1)) {
                continue;
            }

            // See if we can generate an edge for the forward Graph.
            {
                const double logP = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing0);
                if(logP >= logPThreshold) {
                    Graph& graph = graphs[0];
                    EdgeToBeAdded& edgeToBeAdded = edgesToBeAdded[0].emplace_back();
                    edgeToBeAdded.v0 = graph.vertexMap.at(segment0);
                    edgeToBeAdded.v1 = graph.vertexMap.at(segment1);
                    Edge& edge = edgeToBeAdded.edge;
                    edge.commonCount = segmentPairInformation.commonCount;
                    edge.missingCount = segmentPairInformation.missing0;
                    edge.logP = logP;
                    edge.weight = std::pow(10., -0.1 * logP);
                }
            }

            // See if we can generate an edge for the backward Graph.
            {
                const double logP = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing1);
                if(logP >= logPThreshold) {
                    Graph& graph = graphs[1];
                    EdgeToBeAdded& edgeToBeAdded = edgesToBeAdded[1].emplace_back();
                    // The edge is in the opposie direction.
                    edgeToBeAdded.v0 = graph.vertexMap.at(segment1);
                    edgeToBeAdded.v1 = graph.vertexMap.at(segment0);
                    Edge& edge = edgeToBeAdded.edge;
                    edge.commonCount = segmentPairInformation.commonCount;
                    edge.missingCount = segmentPairInformation.missing0;
                    edge.logP = logP;
                    edge.weight = std::pow(10., -0.1 * logP);
                }
            }

        }
    }

    // Now grab the mutex and add the edges we found.
    std::lock_guard<std::mutex> lock(mutex);
    for(uint64_t direction=0; direction<2; direction++) {
        Graph& graph = graphs[direction];
        for(const EdgeToBeAdded& edgeToBeAdded: edgesToBeAdded[direction]) {
            add_edge(edgeToBeAdded.v0, edgeToBeAdded.v1, edgeToBeAdded.edge, graph);
        }
    }
}



// Prune removes all vertices that are not accessible from long
// vertices in both directions.
void Graph::prune()
{
    Graph& graph = *this;

    // Loop over both directions.
    array<std::set<vertex_descriptor>, 2> reachedVertices;
    vector<vertex_descriptor> neighbors;
    for(uint64_t direction=0; direction<2; direction++) {

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        BGL_FORALL_VERTICES(v, graph, Graph) {
            if(graph[v].isLong) {
                q.push(v);
                reachedVertices[direction].insert(v);
            }
        }

        // BFS loop in this direction.
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();

            // Find the neighbors in this direction.
            neighbors.clear();
            if(direction == 0) {
                // Forward.
                BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
                    const vertex_descriptor v1 = target(e, graph);
                    neighbors.push_back(v1);
                }
            } else {
                // Backward.
                BGL_FORALL_INEDGES(v0, e, graph, Graph) {
                    const vertex_descriptor v1 = source(e, graph);
                    neighbors.push_back(v1);
                }
            }

            // Loop over the neighbors.
            for(const vertex_descriptor v1: neighbors) {
                if(not reachedVertices[direction].contains(v1)) {
                    reachedVertices[direction].insert(v1);
                    q.push(v1);
                }
            }
        }
    }


    // Remove vertices that are not reachable in both directions.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not (reachedVertices[0].contains(v) and reachedVertices[1].contains(v))) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexMap.erase(graph[v].segment);
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

    // Sanity check: all leafs must be long vertices.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
        if(isLeaf) {
            SHASTA2_ASSERT(graph[v].isLong);
        }
    }


}



void Graph::writeGraphviz(const AssemblyGraph& assemblyGraph, const string& name) const
{
    const Graph& graph = *this;

    ofstream dot("ReadFollowing-" + name + ".dot");
    dot << "digraph ReadFollowing {\n";

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
            "\"";

        // Color.
        string color;
        if(vertex.isLong) {
            color = "cyan";
        }
        if(not color.empty()) {
            dot << " style=filled fillcolor=" << color;
        }

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
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


        // Tooltip.
        dot << " tooltip=\"" <<
            edge.commonCount << "/" <<
            edge.missingCount << "/" <<
            std::fixed << std::setprecision(1) <<
            edge.logP << "\"";

        // Thickness is determined to logP.
        double logPClipped = max(1., edge.logP);
        logPClipped = min(100., logPClipped);
        const double thickness = 0.05 * logPClipped;
        dot << std::fixed << std::setprecision(2) << " penwidth=" << thickness;

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";

    }

    dot << "}\n";
}



// This finds a shortest path starting at segment0 and ending at a long Segment,
// with path length defined by Edge::weight.
void ReadFollower::findShortestPath(
    Segment segment0,
    uint64_t direction,     // 0 = forward, 1 = backward
    vector<Segment>& path
    ) const
{
    if(direction == 0) {
        findShortestPathForward(segment0, path);
    } else {
        findShortestPathBackward(segment0, path);
    }
}



void ReadFollower::findShortestPathForward(
    Segment segment0,
    vector<Segment>& path
    ) const
{
    graphs[0].findShortestPath(segment0, path);
}



void ReadFollower::findShortestPathBackward(
    Segment segment0,
    vector<Segment>& path
    ) const
{
    graphs[1].findShortestPath(segment0, path);
    std::ranges::reverse(path);
}



void ReadFollower::findAndWriteShortestPath(Segment segment0, uint64_t direction) const
{

    vector<Segment> path;
    findShortestPath(segment0, direction, path);

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const Segment segment: path) {
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;

}



void Graph::findShortestPath(Segment segment0, vector<Segment>& path) const
{
    using namespace boost;
    const Graph& graph = *this;

    // Find the vertex corresponding to this segment.
    const vertex_descriptor v0 = vertexMap.at(segment0);

    // If no outgoing edges, return a path consisting of just segment0.
    if(out_degree(v0, graph) == 0) {
        path.clear();
        path.push_back(segment0);
        return;
    }

    // An exception class used to stop the shortest path computation
    // when a long vertex is encountered.
    class LongVertexReached {
    public:
        vertex_descriptor v;
        LongVertexReached(vertex_descriptor v) : v(v) {}
    };

    // The DijkstraVisitor class throws LongVertexReached when a long vertex is encountered.
    class DijkstraVisitor : public boost::dijkstra_visitor<> {
    public:
        vertex_descriptor v0;
        DijkstraVisitor(vertex_descriptor v0) : v0(v0) {}
        void examine_vertex(vertex_descriptor v, const Graph& graph)
        {
            if((v != v0) and (graph[v].isLong)) {
                throw LongVertexReached(v);
            }
        }
    };
    DijkstraVisitor dijkstraVisitor(v0);

    // The predecessorMap is filled in by the call to dag_shortest_paths
    // and can be used to reconstruct the path from v0 to
    // the first long edge encountered.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

    // Compute the shortest path using Edge::weight.
    vertex_descriptor v1 = null_vertex();
    try {
        dijkstra_shortest_paths(graph, v0,
           weight_map(boost::get(&Edge::weight, graph)).
           vertex_index_map(make_assoc_property_map(vertexIndexMap)).
           predecessor_map(make_assoc_property_map(predecessorMap)).
           visitor(dijkstraVisitor)
           );
    } catch(LongVertexReached& longVertexReached) {
        v1 = longVertexReached.v;
    } catch(std::exception& e) {
        SHASTA2_ASSERT(0);
    }

    // Use the predecessor map to construct the path.
    path.clear();
    vertex_descriptor v = v1;
    while(true) {
        path.push_back(graph[v].segment);
        if(v == v0) {
            break;
        }
        v = predecessorMap[v];
    }
    std::ranges::reverse(path);
}



void Graph::createVertexIndexMap()
{
    Graph& graph = *this;

    vertexIndexMap.clear();
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }
}

