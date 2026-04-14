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
            " has " << num_vertices(searchGraphs[direction]) <<
            " vertices and " << num_edges(searchGraphs[direction]) << " edges." << endl;
    }

    // Prune.
    for(uint64_t direction=0; direction<2; direction++) {
        searchGraphs[direction].prune();
    }
    searchGraphs[0].writeGraphviz(assemblyGraph, "Pruned");

    for(uint64_t direction=0; direction<2; direction++) {
        cout << "After pruning, the read following graph for direction " << direction <<
            " has " << num_vertices(searchGraphs[direction]) <<
            " vertices and " << num_edges(searchGraphs[direction]) << " edges." << endl;
    }

    cout << "The read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    graph.writeGraphviz(assemblyGraph, "A");

    // Before we can compute shortest paths we have to create the vertex index map
    // for each of the two graphs.
    for(uint64_t direction=0; direction<2; direction++) {
        searchGraphs[direction].createVertexIndexMap();
    }

    // Use the SearchGraphs to find shortest paths between long segments
    // and store them in the Graph.
    findShortestPaths();

    cout << "After finding shortest paths, the read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    graph.writeGraphviz(assemblyGraph, "B");

    graph.removeWeakEdges();
    graph.writeGraphviz(assemblyGraph, "C");
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



void ReadFollower::createVertices()
{
    const uint64_t lengthThreshold = assemblyGraph.options.readFollowingSegmentLengthThreshold;

    // Loop over all AssemblyGraph Segments.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {

        // Get the length and check if it qualifies as long.
        const uint64_t length = assemblyGraph[segment].length();
        const bool isLong = (length >= lengthThreshold);

        // Add vertices to the SearchGraphs.
        for(uint64_t direction=0; direction<2; direction++) {
            searchGraphs[direction].createVertex(segment, length, isLong);
        }

        // If isLong, also add a vertex to the Graph.
        if(isLong) {
            graph.createVertex(segment, length);
        }
    }
}



void SearchGraph::createVertex(Segment segment, uint64_t length, bool isLong)
{
    SHASTA2_ASSERT(not vertexMap.contains(segment));
    SearchGraph& graph = *this;
    const vertex_descriptor v = add_vertex(SearchGraphVertex(segment, length, isLong), graph);
    vertexMap.insert(make_pair(segment, v));
}



SearchGraphVertex::SearchGraphVertex(
    Segment segment,
    uint64_t length,
    bool isLong) :
    segment(segment),
    length(length),
    isLong(isLong)
{
}



void Graph::createVertex(Segment segment, uint64_t length)
{
    SHASTA2_ASSERT(not vertexMap.contains(segment));
    Graph& graph = *this;
    const vertex_descriptor v = add_vertex(GraphVertex(segment, length), graph);
    vertexMap.insert(make_pair(segment, v));
}



GraphVertex::GraphVertex(
    Segment segment,
    uint64_t length) :
    segment(segment),
    length(length)
{
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

    // No html output from analyzeSegmentPair.
    ostream html(0);



    class SegmentPair {
    public:
        Segment segment0;
        Segment segment1;
        SegmentPairInformation segmentPairInformation;
        double logP;
        double logPForward;
        double logPBackward;
        SegmentPair(
            Segment segment0,
            Segment segment1,
            const SegmentPairInformation& segmentPairInformation) :
                segment0(segment0),
                segment1(segment1),
                segmentPairInformation(segmentPairInformation)
        {
            logP         = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing());
            logPForward  = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing0);
            logPBackward = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing1);
        }
    };

    // SegmentPairs that will generate edes in the SearchGraphs.
    vector<SegmentPair> edgesToBeAdded;

    // SegmentPairs that will generate edges in the Graph.
    vector<SegmentPair> edgesToBeAddedLong;



    // Loop over batches of segment pairs assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segment pairs in this batch.
        for(uint64_t i=begin; i<end; i++) {
            const auto& [segment0, segment1] = segmentPairs[i];

            const SegmentPairInformation segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
                html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

            // If the offset is negative, discard it.
            if(segmentPairInformation.segmentOffset < 0) {
                continue;
            }

            // Check for long segments.
            const bool isLong0 = isLong(segment0);
            const bool isLong1 = isLong(segment1);

            // Tentatively create the SegmentPair without storing it.
            const SegmentPair segmentPair(segment0, segment1, segmentPairInformation);

            // If it does not satisfy our requirements, discard it.
            if(segmentPair.segmentPairInformation.commonCount < minCommonCount) {
                continue;
            }
            if(
                (segmentPair.logP < logPThreshold) and
                (segmentPair.logPForward < logPThreshold) and
                (segmentPair.logPBackward < logPThreshold)) {
                continue;
            }
            if(not assemblyGraph.canConnect(segment0, segment1)) {
                continue;
            }

            // See if we can add it to the SearchGraphs
            if(segmentPair.logP > logPThreshold) {
                edgesToBeAdded.emplace_back(segmentPair);
            }

            // See if we can add it to the long graph.
            if(isLong0 and isLong1) {
                if( (segmentPair.logP         > logPThreshold) or
                    (segmentPair.logPForward  > logPThreshold) or
                    (segmentPair.logPBackward > logPThreshold)) {
                    edgesToBeAddedLong.emplace_back(segmentPair);
                }
            }

        }
    }



    // Now grab the mutex and add the edges we found.
    std::lock_guard<std::mutex> lock(mutex);
    for(const SegmentPair& segmentPair: edgesToBeAdded) {

        const SearchGraphEdge edge(
            segmentPair.segmentPairInformation.commonCount,
            segmentPair.segmentPairInformation.missing0,
            segmentPair.segmentPairInformation.missing1);

        // Add it to the forward graph.
        SearchGraph& forwardGraph = searchGraphs[0];
        add_edge(
            forwardGraph.vertexMap.at(segmentPair.segment0),
            forwardGraph.vertexMap.at(segmentPair.segment1),
            edge,
            forwardGraph);

        // Add it to the backward graph, reversing the direction.
        SearchGraph& backwardGraph = searchGraphs[1];
        add_edge(
            backwardGraph.vertexMap.at(segmentPair.segment1),
            backwardGraph.vertexMap.at(segmentPair.segment0),
            edge,
            backwardGraph);
    }



    for(const SegmentPair& segmentPair: edgesToBeAddedLong) {

        const GraphEdge edge(
            segmentPair.segmentPairInformation.commonCount,
            segmentPair.segmentPairInformation.missing0,
            segmentPair.segmentPairInformation.missing1);

        add_edge(
            graph.vertexMap.at(segmentPair.segment0),
            graph.vertexMap.at(segmentPair.segment1),
            edge,
            graph);
    }
}



// Prune removes all vertices that are not accessible from long
// vertices in both directions.
void SearchGraph::prune()
{
    SearchGraph& graph = *this;

    // Loop over both directions.
    array<std::set<vertex_descriptor>, 2> reachedVertices;
    vector<vertex_descriptor> neighbors;
    for(uint64_t direction=0; direction<2; direction++) {

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        BGL_FORALL_VERTICES(v, graph, SearchGraph) {
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
                BGL_FORALL_OUTEDGES(v0, e, graph, SearchGraph) {
                    const vertex_descriptor v1 = target(e, graph);
                    neighbors.push_back(v1);
                }
            } else {
                // Backward.
                BGL_FORALL_INEDGES(v0, e, graph, SearchGraph) {
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
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
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
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
        if(isLeaf) {
            SHASTA2_ASSERT(graph[v].isLong);
        }
    }


}



void SearchGraph::writeGraphviz(
    const AssemblyGraph& assemblyGraph,
    const string& name) const
{
    const SearchGraph& graph = *this;

    ofstream dot("ReadFollowing-SearchGraph-" + name + ".dot");
    dot << "digraph ReadFollowingSearchGraph {\n";

    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        const SearchGraphVertex& vertex = graph[v];
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



    BGL_FORALL_EDGES(e, graph, SearchGraph) {
        const SearchGraphEdge& edge = graph[e];

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
            edge.missingCount0 << "/" <<
            edge.missingCount1 << "/" <<
            std::fixed << std::setprecision(1) <<
            edge.logP << "\"";

        // Thickness is determined to logP.
        // Color is always black.
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



void Graph::writeGraphviz(
    const AssemblyGraph& assemblyGraph,
    const string& name) const
{
    const Graph& graph = *this;

    ofstream dot("ReadFollowing-Graph-" + name + ".dot");
    dot <<
        "digraph ReadFollowingGraph {\n"
        "tooltip=\" \";\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const GraphVertex& vertex = graph[v];
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

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    // Each edge is displayed as one or two display edges:
    // - If edge.hasDirectConnection(), an edge representing the
    //   DirectConnectInformation is drawn. This edge uses a solid line.
    //   The edge color depends on edge.directConnectionType():
    //   * If Bidirectional, the edge is black.
    //   * If forward, the edge is blue.
    //   * If backward, the edge is green.
    //   * If ambiguous, the edge is red.
    //   The arrows can be empty or filled:
    //   * The source arrow is filled if directConnectInformation.logPForward  >= logPThreshold.
    //   * The target arrow is filled if directConnectInformation.logPBackward >= logPThreshold.
    // - If one or both of the edge assemblyPaths are non-empty, an edge representing
    //   these assembly paths is drawn. This edge uses a dashed line.
    //   The edge color depends on the two assembly paths:
    //   * If both assembly paths are non-trivial, the edge is black.
    //   * If only the forward  assembly path is non-trivial, the edge is blue.
    //   * If only the backward assembly path is non-trivial, the edge is green.
    //   The arrows can be empty or filled:
    //   * The source arrow is filled if the forward  assembly path is non-trivial.
    //   * The target arrow is filled if the backward assembly path is non-trivial.

    BGL_FORALL_EDGES(e, graph, Graph) {
        const GraphEdge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;


        // Draw the edge representing the DirectConnectInformation.
        if(edge.hasDirectConnection()) {

            dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

            // Begin attributes.
            dot << "[";

            // Tooltip.
            dot << " tooltip=\"" <<
                edge.directConnectInformation.commonCount << "/" <<
                edge.directConnectInformation.missingCount0 << "/" <<
                edge.directConnectInformation.missingCount1 << "/" <<
                std::fixed << std::setprecision(1) <<
                edge.directConnectInformation.logP << "/" <<
                edge.directConnectInformation.logPForward << "/" <<
                edge.directConnectInformation.logPBackward << "\"";

            // Thickness is determined to maxLogP.
            double logPClipped = max(1., edge.directConnectInformation.maxLogP());
            logPClipped = min(100., logPClipped);
            const double thickness = 0.05 * logPClipped;
            dot << std::fixed << std::setprecision(2) << " penwidth=" << thickness;

            // Color depends on the edge type.
            string color;
            switch(edge.directConnectionType()) {
            case GraphEdge::DirectConnectionType::None:
                color = "Red";
                break;
            case GraphEdge::DirectConnectionType::Bidirectional:
                color = "Black";
                break;
            case GraphEdge::DirectConnectionType::Forward:
                color = "Blue";
                break;
            case GraphEdge::DirectConnectionType::Backward:
                color = "Green";
                break;
            case GraphEdge::DirectConnectionType::Ambiguous:
                color = "Orange";
                break;
            default:
                SHASTA2_ASSERT(0);
            }
            dot << " color=" << color;

            // The source arrow is filled if directConnectInformation.logPForward  >= logPThreshold.
            // The target arrow is filled if directConnectInformation.logPBackward >= logPThreshold.
            dot << " dir=both";
            if(edge.directConnectInformation.logPForward  >= logPThreshold) {
                dot << " arrowtail=inv";
            } else {
                dot << " arrowtail=oinv";
            }
            if(edge.directConnectInformation.logPBackward  >= logPThreshold) {
                dot << " arrowhead=normal";
            } else {
                dot << " arrowhead=onormal";
            }

            // End attributes.
            dot << "]";

            // End the line for this edge.
            dot << ";\n";
        }



        // Draw the edge representing the assembly paths.
        if((edge.assemblyPaths[0].size() > 2) or (edge.assemblyPaths[1].size()  > 2)) {

            dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

            // Begin attributes.
            dot << "[";

            dot << "style=dashed";

            // Tooltip.
            dot << " tooltip=\"" <<
                edge.assemblyPaths[0].size() << "/" <<
                edge.assemblyPaths[1].size() << "\"";

            // The arrowtail is filled if the forward path exists.
            // The arrowhead is filled if the backward path exists.
            dot << " dir=both";
            if(edge.assemblyPaths[0].size() <= 2) {
                dot << " arrowtail=oinv";
            } else {
                dot << " arrowtail=inv";
            }
            if(edge.assemblyPaths[1].size() <= 2) {
                dot << " arrowhead=onormal";
            } else {
                dot << " arrowhead=normal";
            }

            // Color.
            string color = "Black";
            if(edge.assemblyPaths[0].size() <= 2) {
                color = "Green";
            }
            if(edge.assemblyPaths[1].size() <= 2) {
                color = "Blue";
            }
            dot << " color=" << color;


            // End attributes.
            dot << "]";

            // End the line for this edge.
            dot << ";\n";
        }

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
    searchGraphs[0].findShortestPath(segment0, path);
}



void ReadFollower::findShortestPathBackward(
    Segment segment0,
    vector<Segment>& path
    ) const
{
    searchGraphs[1].findShortestPath(segment0, path);
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



void SearchGraph::findShortestPath(Segment segment0, vector<Segment>& path) const
{
    using namespace boost;
    const SearchGraph& graph = *this;

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
        void examine_vertex(vertex_descriptor v, const SearchGraph& graph)
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
           weight_map(boost::get(&SearchGraphEdge::weight, graph)).
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



void SearchGraph::createVertexIndexMap()
{
    SearchGraph& graph = *this;

    vertexIndexMap.clear();
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }
}




void ReadFollower::findAndWriteAssemblyPaths() const
{
    vector< vector<Segment> > assemblyPaths;
    findAssemblyPaths(assemblyPaths);
    writeAssemblyPaths(assemblyPaths);
}



void ReadFollower::writeAssemblyPaths(const vector< vector<Segment> >&) const
{

}



void ReadFollower::findAssemblyPaths(vector< vector<Segment> >&) const
{
}



SearchGraphEdge::SearchGraphEdge(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    commonCount(commonCount),
    missingCount0(missingCount0),
    missingCount1(missingCount1)
{
    logP         = a * double(commonCount) - b * double(missingCount0 + missingCount1);
    weight = pow(10., 0.1 * logP);
}


GraphEdge::GraphEdge(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    directConnectInformation(commonCount, missingCount0, missingCount1)
{}



GraphEdge::DirectConnectInformation::DirectConnectInformation(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    commonCount(commonCount),
    missingCount0(missingCount0),
    missingCount1(missingCount1)
{
    logP         = a * double(commonCount) - b * double(missingCount0 + missingCount1);
    logPForward  = a * double(commonCount) - b * double(missingCount0);
    logPBackward = a * double(commonCount) - b * double(missingCount1);
}



bool ReadFollower::isLong(Segment segment) const
{
    return graph.vertexMap.contains(segment);
}



double GraphEdge::DirectConnectInformation::maxLogP() const
{
    return max(logP, max(logPForward, logPBackward));
}



GraphEdge::DirectConnectionType GraphEdge::directConnectionType() const
{
    if(not hasDirectConnection()) {
        return DirectConnectionType::None;
    }

    if(directConnectInformation.logP >= logPThreshold) {
        return DirectConnectionType::Bidirectional;
    }

    if(directConnectInformation.logPForward >= logPThreshold) {
        if(directConnectInformation.logPBackward < logPThreshold) {
            return DirectConnectionType::Forward;
        } else {
            return DirectConnectionType::Ambiguous;
        }
    }

    if(directConnectInformation.logPBackward >= logPThreshold) {
        if(directConnectInformation.logPForward < logPThreshold) {
            return DirectConnectionType::Backward;
        } else {
            return DirectConnectionType::Ambiguous;
        }
    }

    return DirectConnectionType::Ambiguous;
}



// Use the SearchGraphs to find shortest paths between long segments
// and store them in the Graph.
void ReadFollower::findShortestPaths()
{
    // Add edges.
    vector<Segment> path;
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        const Segment segment0 = graph[v0].segment;

        // Loop for shortest paths in both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            findShortestPath(segment0, direction, path);

            // Discard a trivial path.
            if(path.size() < 2) {
                continue;
            }

            const Segment s0 = path.front();
            const Segment s1 = path.back();
            if(direction == 0) {
                SHASTA2_ASSERT(s0 == segment0);
            } else {
                SHASTA2_ASSERT(s1 == segment0);
            }

            const Graph::vertex_descriptor u0 = graph.vertexMap.at(s0);
            const Graph::vertex_descriptor u1 = graph.vertexMap.at(s1);


            // Look for a Graph edge between the u0 and u1.
            Graph::edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(u0, u1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = boost::add_edge(u0, u1, graph);
            }
            SHASTA2_ASSERT(edgeExists);

            // Store the path.
            graph[e].assemblyPaths[direction] = path;
        }
    }

}



// This removes edges without a direct connection
// and that don't have paths in both directions.
void Graph::removeWeakEdges()
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const GraphEdge& edge = graph[e];

        // If it has a direct connection, keep it.
        if(edge.hasDirectConnection()) {
            continue;
        }

        // Otherwise, only keep it if both assembly paths are present.
        if(edge.assemblyPaths[0].empty() or edge.assemblyPaths[1].empty()) {
            edgesToBeRemoved.push_back(e);
        }
    }


    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}
