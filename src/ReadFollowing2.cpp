// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing2.hpp"
#include "Journeys.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
using namespace shasta2;
using namespace ReadFollowing2;

// Standard library.
#include <queue>



// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Graph>;



Graph::Graph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<Graph>(*this),
    assemblyGraph(assemblyGraph)
{
    Graph& graph = *this;

    createVertices();
    createEdgeCandidates();
    createEdgesMultithreaded();
    cout << "The initial read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    write("Initial");

    prune();

    cout << "The final read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    write("Final");

    fillConnectivity();

}



void Graph::fillConnectivity()
{
    Graph& graph = *this;

    // For better reproducibility, store children/parent vertices sorted by segment id.
    class VertexInfo {
    public:
        vertex_descriptor v;
        uint64_t commonCount;
        uint64_t segmentId;
        bool operator<(const VertexInfo& that) const
        {
            return segmentId < that.segmentId;
        }
    };
    vector<VertexInfo> vertexInfos;

    // Loop over all vertices in the Graph.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        Vertex& vertex0 = graph[v0];

        // Gather the children.
        vertexInfos.clear();
        BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);
            const uint64_t commonCount = graph[e].segmentPairInformation.commonCount;
            VertexInfo& vertexInfo = vertexInfos.emplace_back();
            vertexInfo.v = v1;
            vertexInfo.commonCount = commonCount;
            vertexInfo.segmentId = segmentId(v1);
        }

        // Sort them by segment id.
        sort(vertexInfos.begin(), vertexInfos.end());

        // Store them in the vertex.
        // Each child is stored commonCount times.
        vertex0.children.clear();
        for(const VertexInfo& vertexInfo: vertexInfos) {
            const uint64_t commonCount = vertexInfo.commonCount;
            for(uint64_t i=0; i<commonCount; i++) {
                vertex0.children.push_back(vertexInfo.v);
            }
        }



        // Gather the parents.
        vertexInfos.clear();
        BGL_FORALL_INEDGES(v0, e, graph, Graph) {
            const vertex_descriptor v1 = source(e, graph);
            const uint64_t commonCount = graph[e].segmentPairInformation.commonCount;
            VertexInfo& vertexInfo = vertexInfos.emplace_back();
            vertexInfo.v = v1;
            vertexInfo.commonCount = commonCount;
            vertexInfo.segmentId = segmentId(v1);
        }

        // Sort them by segment id.
        sort(vertexInfos.begin(), vertexInfos.end());

        // Store them in the vertex.
        // Each child is stored commonCount times.
        vertex0.parents.clear();
        for(const VertexInfo& vertexInfo: vertexInfos) {
            const uint64_t commonCount = vertexInfo.commonCount;
            for(uint64_t i=0; i<commonCount; i++) {
                vertex0.parents.push_back(vertexInfo.v);
            }
        }

    }
}



// Prune removes all vertices that are not accessible from long
// vertices in both directions.
void Graph::prune()
{
    Graph& graph = *this;
    const uint64_t lengthThreshold = assemblyGraph.options.readFollowingSegmentLengthThreshold;

    // Loop over both directions.
    array<std::set<vertex_descriptor>, 2> reachedVertices;
    vector<vertex_descriptor> neighbors;
    for(uint64_t direction=0; direction<2; direction++) {

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        for(const vertex_descriptor v: longVertices) {
            q.push(v);
            reachedVertices[direction].insert(v);
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

        cout << "Found " << reachedVertices[direction].size() <<
            " reachable vertices in direction " << direction << endl;
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
    cout << "Pruned " << verticesToBeRemoved.size() <<
        " vertices from the read following graph." << endl;

    // Sanity check: all leafs must be long vertices.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
        if(isLeaf) {
            SHASTA2_ASSERT(graph[v].length >= lengthThreshold);
        }
    }

}




// Create vertices of the ReadFollowing1 graph.
// Each vertex corresponds to a Segment of the AssemblyGraph.
void Graph::createVertices()
{
    Graph& graph = *this;
    const uint64_t lengthThreshold = assemblyGraph.options.readFollowingSegmentLengthThreshold;

    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
        vertexMap.insert(make_pair(segment, v));

        if(graph[v].length >= lengthThreshold) {
            graph[v].isLongVertex = true;
            longVertices.insert(v);
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



// Create edge candidates of the ReadFollowing1 graph.
// An edge candidate v0->v1 is created if the final support of v0
// shares at least one OrientedReadId with the initial support of v1.
void Graph::createEdgeCandidates()
{
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdgeCandidates begins." << endl;

    Graph& graph = *this;
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;

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
    // shares at least minCommonCount OrientedReadId with the initial support of v1.
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
                edgeCandidates.push_back(EdgeCandidate({v0, v1}));
            }
        }
    }

    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(edgeCandidates, count, minCommonCount);
    SHASTA2_ASSERT(edgeCandidates.size() == count.size());

    cout << "Found " << edgeCandidates.size() << " candidate edges for the read following graph." << endl;
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdgeCandidates ends." << endl;
}



// Generate an edge for each of edge candidate that satisfies our requirements.
void Graph::createEdgesMultithreaded()
{
    uint64_t threadCount = assemblyGraph.options.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    setupLoadBalancing(edgeCandidates.size(), 100);
    runThreads(&Graph::createEdgesThreadFunction, threadCount);
}



void Graph::createEdgesThreadFunction([[maybe_unused]] uint64_t threadId)
{
    Graph& graph = *this;
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const double minCorrectedJaccard = assemblyGraph.options.readFollowingMinCorrectedJaccard;

    // Prepare a vector of edges to be added.
    // We will add them all at the end so we only have to acquire the mutex once.
    class EdgeToBeAdded {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        Edge edge;
        EdgeToBeAdded(
            const EdgeCandidate& edgeCandidate,
            const Graph& graph
            ) :
            v0(edgeCandidate.v0),
            v1(edgeCandidate.v1),
            edge(graph.assemblyGraph, graph[v0].segment, graph[v1].segment)
        {}
    };
    vector<EdgeToBeAdded> edgesToBeAdded;

    // Loop over batches of candidate edges assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate edges assigned to this batch.
        for(uint64_t i=begin; i<end; i++) {
            const EdgeCandidate& edgeCandidate = edgeCandidates[i];

            // Tentatively add it to the edges to be added.
            EdgeToBeAdded& edgeToBeAdded = edgesToBeAdded.emplace_back(edgeCandidate, graph);
            const Edge& edge = edgeToBeAdded.edge;

            // This must be true given the way we constructed the vertex pairs.
            SHASTA2_ASSERT(edge.segmentPairInformation.commonCount >= minCommonCount);

            // If it does not satisfy our requirements, get rid of it.
            if(edge.segmentPairInformation.segmentOffset < 0) {
                edgesToBeAdded.pop_back();
                continue;
            }
            if(edge.segmentPairInformation.correctedJaccard < minCorrectedJaccard) {
                edgesToBeAdded.pop_back();
                continue;
            }
            if(not assemblyGraph.canConnect(graph[edgeCandidate.v0].segment, graph[edgeCandidate.v1].segment)) {
                edgesToBeAdded.pop_back();
                continue;
            }
        }
    }


    // Now grab the mutex and add the edges we found.
    std::lock_guard<std::mutex> lock(mutex);
    for(const EdgeToBeAdded& edgeToBeAdded: edgesToBeAdded) {
        add_edge(edgeToBeAdded.v0, edgeToBeAdded.v1, edgeToBeAdded.edge, graph);
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



uint64_t Graph::segmentId(vertex_descriptor v) const
{
    const Graph& graph = *this;
    const Segment segment = graph[v].segment;
    return assemblyGraph[segment].id;
}



void Graph::write(const string& name) const
{
    cout << "ReadFollowing1-" << name << ": " << num_vertices(*this) <<
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
        // const int32_t offset = edge.segmentPairInformation.segmentOffset;

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
            edge.segmentPairInformation.correctedJaccard << "\"";

        // Tooltip.
        dot << " tooltip=\"" <<
            edge.segmentPairInformation.commonCount << "/" <<
            std::fixed << std::setprecision(2) <<
            edge.segmentPairInformation.correctedJaccard << "\"";

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



void Graph::writeRandomPath(Segment segment, uint64_t direction)
{
    const Graph& graph = *this;

    const auto it = vertexMap.find(segment);
    SHASTA2_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;

    Path path;
    std::random_device randomGenerator;
    findRandomPath(v, direction, randomGenerator, path);

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;
}



void Graph::writePathStatistics(Segment segment, uint64_t direction)
{

    // Find the vertex corresponding to this segment.
    const auto it = vertexMap.find(segment);
    SHASTA2_ASSERT(it != vertexMap.end());
    const vertex_descriptor v0 = it->second;

    // Create random paths and record the terminal vertices.
    vector<vertex_descriptor> terminalVertices;
    Path path;
    std::random_device randomGenerator;
    for(uint64_t i=0; i<pathCount; i++) {
        findRandomPath(v0, direction, randomGenerator, path);
        const vertex_descriptor v1 = (direction == 0) ? path.back() : path.front();
        terminalVertices.push_back(v1);
    }

    // Count the occurrences of each terminal vertex.
    vector<uint64_t> count;
    deduplicateAndCount(terminalVertices, count);
    SHASTA2_ASSERT(terminalVertices.size() == count.size());

    // Write out.
    for(uint64_t i=0; i<terminalVertices.size(); i++) {
        const vertex_descriptor v1 = terminalVertices[i];
        cout << segmentId(v1) << " " << count[i] << endl;
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomPath(
    vertex_descriptor v, uint64_t direction,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path)
{
    if(direction == 0) {
        findRandomForwardPath(v, randomGenerator, path);
    } else {
        findRandomBackwardPath(v, randomGenerator, path);
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomForwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path)
{
    const Graph& graph = *this;
    // static ofstream debugOut("findRandomForwardPath.txt");
    // debugOut << "Starting a forward random path at " << segmentId(v) << "\n";

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);

    // At each iteration, add one vertex to the path.
    while(true) {

        // Pick a random child.
        const vector<vertex_descriptor>& children = graph[v].children;
        const uint64_t n = children.size();
        if(n ==0) {
            break;
        }
        std::uniform_int_distribution<uint64_t> distribution(0, n - 1);
        const uint64_t i = distribution(randomGenerator);
        v = children[i];

        // Add this vertex to the path and continue from here.
        path.push_back(v);

        // debugOut << n << " " << i << " " << segmentId(v) << "\n";

        if(graph[v].isLongVertex) {
            break;
        }
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomBackwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path)
{
    const Graph& graph = *this;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);

    // At each iteration, add one vertex to the path.
    while(true) {

        // Pick a random parent.
        const vector<vertex_descriptor>& children = graph[v].parents;
        const uint64_t n = children.size();
        if(n ==0) {
            break;
        }
        std::uniform_int_distribution<uint64_t> distribution(0, n - 1);
        const uint64_t i = distribution(randomGenerator);
        v = children[i];

        // Add this vertex to the path and continue from here.
        path.push_back(v);

        if(graph[v].isLongVertex) {
            break;
        }
    }

    // Reverse the path so it goes forward.
    std::ranges::reverse(path);
}



// Create pathCount paths in each starting direction and for each starting vertex.
// Store in each vertex the terminal vertices of the random paths.
void Graph::findRandomPaths()
{
    Graph& graph = *this;
    std::mt19937 randomGenerator;

    Path path;
    vector<uint64_t> count;
    vector<vertex_descriptor> terminalVertices;

    // Loop over all vertices.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        Vertex& vertex0 = graph[v0];

        // Loop over the two directions for this vertex.
        for(uint64_t direction=0; direction<2; direction++) {
            terminalVertices.clear();

            // Generate pathCount random paths starting at v0, in this direction.
            for(uint64_t i=0; i<pathCount; i++) {
                findRandomPath(v0, direction, randomGenerator, path);
                const vertex_descriptor v1 = (direction == 0) ? path.back() : path.front();
                terminalVertices.push_back(v1);
            }

            // Count how many times each terminal vertex occurred.
            deduplicateAndCount(terminalVertices, count);
            vertex0.randomPathInfos[direction].clear();
            for(uint64_t j=0; j<terminalVertices.size(); j++) {
                auto& randomPathInfo = vertex0.randomPathInfos[direction].emplace_back();
                randomPathInfo.v = terminalVertices[j];
                randomPathInfo.count = count[j];
                randomPathInfo.segmentId = segmentId(randomPathInfo.v);
            }
            sort(vertex0.randomPathInfos[direction].begin(), vertex0.randomPathInfos[direction].end());
        }
    }



    ofstream csv("RandomPaths.csv");
    csv << "v0,Length0,direction,v1,Length1,count,\n";
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        Vertex& vertex0 = graph[v0];

        for(uint64_t direction=0; direction<2; direction++) {
            for(const auto& randomPathInfo: vertex0.randomPathInfos[direction]) {
                csv << segmentId(v0) << ",";
                csv << vertex0.length << ",";
                csv << direction << ",";
                csv << randomPathInfo.segmentId << ",";
                csv << graph[randomPathInfo.v].length << ",";
                csv << randomPathInfo.count << ",";
                csv << "\n";
            }
        }
    }



}
