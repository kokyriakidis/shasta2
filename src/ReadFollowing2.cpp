// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing2.hpp"
#include "approximateTopologicalSort.hpp"
#include "color.hpp"
#include "DisjointSets.hpp"
#include "longestPath.hpp"
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
    assemblyGraph(assemblyGraph),
    orderById(*this)
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



void Graph::findAndWriteAssemblyPaths()
{
    vector< vector<Segment> > assemblyPaths;
    findAssemblyPaths(assemblyPaths);
    writeAssemblyPaths(assemblyPaths);
}



void Graph::createRandomPathsMap()
{
    Graph& graph = *this;
    const bool debug = true;

    // Find short vertices that are reliably preceded/followed by a single vertex.
    // These will be used to fill in assembly paths.
    randomPathsMap.clear();
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];

        if(vertex.isLongVertex) {
            continue;
        }

        if(vertex.randomPathInfos[0].empty()) {
            continue;
        }

        if(vertex.randomPathInfos[1].empty()) {
            continue;
        }

        const Vertex::RandomPathInfo& nextInfo = vertex.randomPathInfos[0].front();
        const Vertex::RandomPathInfo& previousInfo = vertex.randomPathInfos[1].front();

        if(nextInfo.count < pathCountThreshold2) {
            continue;
        }
        if(previousInfo.count < pathCountThreshold2) {
            continue;
        }

        randomPathsMap[{previousInfo.v, nextInfo.v}].push_back(v);
    }
    if(debug) {
        ofstream csv("RandomPathsMap.csv");
        for(const auto& [p, shortVertices]: randomPathsMap) {
            const vertex_descriptor v0 = p.first;
            const vertex_descriptor v1 = p.second;
            csv << segmentId(v0) << ",";
            csv << segmentId(v1) << ",";
            for(const vertex_descriptor v: shortVertices) {
                csv << segmentId(v) << ",";
            }
            csv << endl;
        }
    }

}



void Graph::findAssemblyPaths(vector< vector<Segment> >& assemblyPaths)
{
    Graph& graph = *this;

    // Create pathCount paths in each starting direction and for each starting vertex.
    // Information on the random paths found is stored in Vertex::randomPathInfos.
    findRandomPaths();

    // Find short vertices that, based on the stored random paths,
    // are reliably preceded/followed by a single vertex.
    // These will be used to fill in assembly paths.
    createRandomPathsMap();

    // Create the PathGraph. It has a vertex for each long segment.
    PathGraph pathGraph(graph);
    pathGraph.create();
    pathGraph.writeGraphviz("Initial");
    cout << "The initial PathGraph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;

    // Compute assembly paths on PathGraph edges.
    // Remove the edges for which this fails.
    pathGraph.findAssemblyPathsOnEdges();
    pathGraph.writeGraphviz("Final");
    cout << "The final PathGraph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;

    // Find the non-trivial connected components of the PathGraph.
    const vector< shared_ptr<PathGraph> > componentPointers = pathGraph.findConnectedComponents();
    cout << "The PathGraph has " << componentPointers.size() << " non-trivial connected components." << endl;

    // Process one connected component at a time.
    assemblyPaths.clear();
    for(const shared_ptr<PathGraph>& componentPointer: componentPointers) {
        vector<Segment> assemblyPath;
        componentPointer->findAssemblyPath(assemblyPath);
        if(assemblyPath.size() > 1) {
            assemblyPaths.push_back(assemblyPath);
        }
    }
}



void Graph::writeAssemblyPaths(const vector< vector<Segment> >& assemblyPaths) const
{
    ofstream csv("AssemblyPaths.csv");
    cout << "Found " << assemblyPaths.size() << " assembly paths. See AssemblyPaths.csv for details." << endl;
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



PathGraph::PathGraph(const Graph& graph) :
    graph(graph)
{
}



void PathGraph::create()
{
    createVertices();
    createEdges();
}



void PathGraph::createVertices()
{
    PathGraph& pathGraph = *this;

    // Each long Segment generates a PathGraphVertex.
    // For better reproducibility, generate the vertices in this order,
    // instead of looping over the longVertices, so they are ordered by segment id.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        if(vertex.isLongVertex) {
            const PathGraph::vertex_descriptor u = boost::add_vertex({vertex.segment}, pathGraph);
            vertexMap.insert({vertex.segment, u});
        }
    }
}



void PathGraph::createEdges()
{
    PathGraph& pathGraph = *this;

    // Loop over all PathGraph vertices.
    BGL_FORALL_VERTICES(u0, pathGraph, PathGraph) {

        // Locate the segment and the corresponding Vertex.
        const Segment segment0 = pathGraph[u0].segment;
        const auto it0 = graph.vertexMap.find(segment0);
        SHASTA2_ASSERT(it0 != graph.vertexMap.end());
        const vertex_descriptor v0 = it0->second;
        const Vertex& vertex0 = graph[v0];

        // Loop over the RandomPathInfos of this vertex in both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            for(const Vertex::RandomPathInfo& randomPathInfo: vertex0.randomPathInfos[direction]) {
                const Graph::vertex_descriptor v1 = randomPathInfo.v;
                const uint64_t count = randomPathInfo.count;

                if(count < Graph::pathCountThreshold1) {
                    continue;
                }

                // Locate the corresponding PathGraphVertex.
                const Segment segment1 = graph[v1].segment;
                const auto it1 = vertexMap.find(segment1);
                SHASTA2_ASSERT(it1 != vertexMap.end());
                const vertex_descriptor u1 = it1->second;

                if(u1 == u0) {
                    continue;
                }

                // Locate the PathGraph edge u0->u1 (if direction==0)
                // or u1->u0 (if direction==1), creating it if necessary.
                edge_descriptor e;
                bool edgeExists = false;
                if(direction == 0) {
                    tie(e, edgeExists) = boost::edge(u0, u1, pathGraph);
                    if(not edgeExists) {
                        tie(e, edgeExists) = boost::add_edge(u0, u1, pathGraph);
                    }
                } else {
                    tie(e, edgeExists) = boost::edge(u1, u0, pathGraph);
                    if(not edgeExists) {
                        tie(e, edgeExists) = boost::add_edge(u1, u0, pathGraph);
                    }
                }
                SHASTA2_ASSERT(edgeExists);
                PathGraphEdge& edge = pathGraph[e];

                // Store the number of random paths found.
                edge.randomPathCount[direction] = count;
            }
        }
    }


}



// Find assembly paths on all edges, then remove the
// edges for which an assembly path could not be found.
void PathGraph::findAssemblyPathsOnEdges()
{
    PathGraph& pathGraph = *this;

    // Vector to contain the edges where we did not find an assembly path.
    // These edges will be removed.
    vector<edge_descriptor> failedEdges;

    // Loop over all edges of the PathGraph.
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        PathGraphEdge& edge = pathGraph[e];

        // Locate the corresponding vertices.
        const vertex_descriptor u0 = source(e, pathGraph);
        const vertex_descriptor u1 = target(e, pathGraph);

        // Find the corresponding Graph vertices.
        const Segment segment0 = pathGraph[u0].segment;
        const Segment segment1 = pathGraph[u1].segment;
        const auto it0 = graph.vertexMap.find(segment0);
        SHASTA2_ASSERT(it0 != graph.vertexMap.end());
        const Graph::vertex_descriptor v0 = it0->second;
        const auto it1 = graph.vertexMap.find(segment1);
        SHASTA2_ASSERT(it1 != graph.vertexMap.end());
        const Graph::vertex_descriptor v1 = it1->second;

        // Find the assembly path between these two grap vertices.
        const bool success = graph.findAssemblyPath(v0, v1, edge.assemblyPath);
        if(not success) {
            failedEdges.push_back(e);
        }
    }

    // Remove the failed edges.
    for(const edge_descriptor e: failedEdges) {
        boost::remove_edge(e, pathGraph);
    }
}



void PathGraph::writeGraphviz(const string& name) const
{
    ofstream dot("PathGraph-" + name + ".dot");
    writeGraphviz(dot);
}



void PathGraph::writeGraphviz(ostream& dot) const
{
    const PathGraph& pathGraph = *this;

    dot << "digraph PathGraph {\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        const Segment segment = pathGraph[v].segment;
        const AssemblyGraphEdge& assemblyGraphEdge = graph.assemblyGraph[segment];
        dot << assemblyGraphEdge.id <<
            " ["
            "label=\"" << assemblyGraphEdge.id <<
            "\\n" << assemblyGraphEdge.length() <<
            "\""
            "]"
            ";\n";
    }



    // Edges.
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const PathGraphEdge& edge = pathGraph[e];

        const PathGraph::vertex_descriptor v0 = source(e, pathGraph);
        const PathGraph::vertex_descriptor v1 = target(e, pathGraph);

        const Segment segment0 = pathGraph[v0].segment;
        const Segment segment1 = pathGraph[v1].segment;


        dot <<
            graph.assemblyGraph[segment0].id << "->" <<
            graph.assemblyGraph[segment1].id <<
            " [label=\"" <<
            edge.randomPathCount[0] << "/" <<
            edge.randomPathCount[1] <<
            "\""
            "];\n";
    }



    dot << "}\n";

}



vector< shared_ptr<PathGraph> > PathGraph::findConnectedComponents()
{
    PathGraph& pathGraph = *this;

    // Map vertices to integers.
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert({v, vertexIndex++});
    }
    const uint64_t n = vertexIndex;

    // Use DisjointSets to compute connected components.
    DisjointSets disjointSets(n);
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const vertex_descriptor v0 = source(e, pathGraph);
        const vertex_descriptor v1 = target(e, pathGraph);
        const uint64_t i0 = vertexIndexMap[v0];
        const uint64_t i1 = vertexIndexMap[v1];
        disjointSets.unionSet(i0, i1);
    }

    // Get the vertex indexes for the non-trivial connected components.
    vector< vector<uint64_t> > componentsVertexIndexes;
    disjointSets.gatherComponents(2, componentsVertexIndexes);
    cout << "Found " << componentsVertexIndexes.size() <<
        " non-trivial connected components of the PathGraph." << endl;



    // Now create the PathGraphs for each of the non-trivial components.
    vector< shared_ptr<PathGraph> > components;
    for(const vector<uint64_t>& componentVertexIndexes: componentsVertexIndexes) {

        // Create the PathGraph for this component.
        const shared_ptr<PathGraph> componentPointer = make_shared<PathGraph>(graph);
        PathGraph& component = *componentPointer;
        components.push_back(componentPointer);

        // Add the vertices, copying them form the PathGraph.
        std::map<vertex_descriptor, vertex_descriptor> componentVertexMap;
        for(const uint64_t i: componentVertexIndexes) {
            const vertex_descriptor v = vertexTable[i];
            const vertex_descriptor u = add_vertex(pathGraph[v], component);
            componentVertexMap.insert({v, u});
        }

        // Add the edges.
        // Use v for vertices of the PathGraph and u for vertices of the component.
        for(const uint64_t i0: componentVertexIndexes) {
            const vertex_descriptor v0 = vertexTable[i0];
            const vertex_descriptor u0 = componentVertexMap[v0];
            BGL_FORALL_OUTEDGES(v0, e, pathGraph, PathGraph) {
                const vertex_descriptor v1 = target(e, pathGraph);
                const vertex_descriptor u1 = componentVertexMap[v1];
                add_edge(u0, u1, pathGraph[e], component);
            }
        }
    }



    return components;
}



// This computes an assembly path, assuming it is working
// on a PathGraph with a single connected component.
void PathGraph::findAssemblyPath(vector<Segment>& assemblyPath)
{
    const bool debug = true;
    assemblyPath.clear();
    PathGraph& component = *this;

    if(debug) {
        cout << "Working on a PathGraph component with " << num_vertices(component) <<
            " vertices and " << num_edges(component) << " edges." << endl;
        cout << "The vertices correspond to segments";
        BGL_FORALL_VERTICES(u, component, PathGraph) {
            cout << " " << segmentId(u);
        }
        cout << endl;
    }

    // Find the edges of longest path.
    // This can throw if this component has cycles.
    vector<edge_descriptor> longestPathEdges;
    try {
        longestPath(component, longestPathEdges);
    } catch(std::exception&) {
        // We can do better.
        if(debug) {
            cout << "No assembly path created for this component because of cycles." << endl;
        }
        return;
    }
    if(debug) {
        cout << "The longest path in the PathGraph has " << longestPathEdges.size() << " edges." << endl;
    }
    SHASTA2_ASSERT(not longestPathEdges.empty());

    // Find the vertices of the longest path.
    vector<vertex_descriptor> longestPathVertices;
    const edge_descriptor firstPathEdge = longestPathEdges.front();
    const vertex_descriptor firstPathVertex = source(firstPathEdge, component);
    longestPathVertices.push_back(firstPathVertex);
    for(const edge_descriptor e: longestPathEdges) {
        const vertex_descriptor v = target(e, component);
        longestPathVertices.push_back(v);
    }
    SHASTA2_ASSERT(longestPathVertices.size() > 1);
    if(debug) {
        cout << "The longest path has " << longestPathVertices.size() << " vertices: " << endl;
        for(const vertex_descriptor u: longestPathVertices) {
            cout << " " << segmentId(u);
        }
        cout << endl;
    }



    // Construct the path by looping over the edges of the longest path.
    for(const edge_descriptor e: longestPathEdges) {
        const vertex_descriptor u0 = source(e, component);
        const Segment segment0 = component[u0].segment;

        // Access the assembly path portion between these two vertices.
        vector<Segment>& internalPortionOfAssemblyPath = component[e].assemblyPath;

        // Update the assembly path.
        // Don't include the last segment to avoid duplications.
        assemblyPath.push_back(segment0);
        std::ranges::copy(internalPortionOfAssemblyPath, back_inserter(assemblyPath));
    }

    // Add the final Segment to the assembly path.
    const PathGraph::edge_descriptor e = longestPathEdges.back();
    const PathGraph::vertex_descriptor u = target(e, component);
    assemblyPath.push_back(component[u].segment);

}



// Find an assembly path between two long vertices.
// This uses the randomPathsMap to locate usable short vertices.
// The Path does not include the segments corresponding to v0 and v1.
// This can fail, in which case it returns false and an empty assembly path.
bool Graph::findAssemblyPath(
    vertex_descriptor v0,
    vertex_descriptor v1,
    vector<Segment>& assemblyPath) const
{
    const bool debug = true;
    const Graph& graph = *this;
    assemblyPath.clear();

    const Vertex& vertex0 = graph[v0];
    const Vertex& vertex1 = graph[v1];
    SHASTA2_ASSERT(vertex0.isLongVertex);
    SHASTA2_ASSERT(vertex1.isLongVertex);

    // Gather the vertices of the Graph subgraph that can be used to assemble this portion.
    std::set<Graph::vertex_descriptor, Graph::OrderById> subgraphVertices(graph.orderById);
    subgraphVertices.insert(v0);
    subgraphVertices.insert(v1);
    const auto it = randomPathsMap.find({v0, v1});
    if(it != randomPathsMap.end()) {
        const vector<vertex_descriptor>& shortVertices = it->second;
        for(const vertex_descriptor v: shortVertices) {
            subgraphVertices.insert(v);
        }
    }

    if(debug) {
        cout << "Findint a possible assembly path portion between " <<
            segmentId(v0) << " and " << segmentId(v1) << endl;
        cout << "The following segments are usable to assemble this portion:" << endl;
        for(const Graph::vertex_descriptor v: subgraphVertices) {
            cout << " " << graph.segmentId(v);
        }
        cout << endl;

        // Also write a csv file that can be loaded in Bandage to visualize these segments.
        ofstream csv(
            "ReadFollowing-Subgraph-" +
            to_string(segmentId(v0)) + "-" +
            to_string(segmentId(v1)) + ".csv");
        csv << "Id,Length,Color,\n";
        for(const Graph::vertex_descriptor v: subgraphVertices) {
            csv << segmentId(v) << ",";
            csv << graph[v].length << ",";
            if(v == v0) {
                csv << "Red,";
            } else if(v == v1) {
                csv << "Green,";
            } else {
                csv << "Blue,";
            }
            csv << "\n";
        }
    }

    // Create a Subgraph of the Graph that uses only these segments.
    Subgraph subgraph(graph, subgraphVertices);
    if(debug) {
        cout << "The subgraph for this portion of the assembly path has " <<
            num_vertices(subgraph) << " vertices and " <<
            num_edges(subgraph) << " edges." << endl;
        subgraph.writeGraphviz(
            "ReadFollowing-Subgraph-" +
            to_string(segmentId(v0)) + "-" +
            to_string(segmentId(v1)) + ".dot");
    }

    // Find the longest path in the subgraph.
    // Remove cycles if necessary.
    vector<Subgraph::edge_descriptor> longestSubgraphPathEdges;
    try {
        longestPath(subgraph, longestSubgraphPathEdges);
    } catch(boost::not_a_dag&) {
        subgraph.makeAcyclic();
        longestPath(subgraph, longestSubgraphPathEdges);
        if(debug) {
            cout << "The subgraph contained cycles, and they were removed." << endl;
        }
    }
    if(longestSubgraphPathEdges.empty()) {
        if(debug) {
            cout << "Graph::findAssemblyPath failed: the longest path is empty." << endl;
        }
        return false;
    }


    // Gather the vertices of the longest path.
    vector<Subgraph::vertex_descriptor> longestSubgraphPathVertices;
    const Subgraph::edge_descriptor firstEdge = longestSubgraphPathEdges.front();
    longestSubgraphPathVertices.push_back(source(firstEdge, subgraph));
    for(const Subgraph::edge_descriptor e: longestSubgraphPathEdges) {
        longestSubgraphPathVertices.push_back(target(e, subgraph));
    }
    if(subgraph[longestSubgraphPathVertices.front()].segment != vertex0.segment) {
        if(debug) {
            cout << "Graph::findAssemblyPath failed: the longest path does not begin at v0." << endl;
        }
        return false;
    }
    if(subgraph[longestSubgraphPathVertices.back()].segment != vertex1.segment) {
        if(debug) {
            cout << "Graph::findAssemblyPath failed: the longest path does not end at v1." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "The longest path in the subgraph has " <<
            longestSubgraphPathVertices.size() << " vertices." << endl;

        // Also write a csv file that can be loaded in Bandage to visualize the longest path.
        ofstream csv(
            "ReadFollowing-Subgraph-LongestPath-" +
            to_string(segmentId(v0)) + "-" +
            to_string(segmentId(v1)) + ".csv");
        csv << "Id,Position,Length,Color,\n";
        for(uint64_t i=0; i<longestSubgraphPathVertices.size(); i++) {
            const Subgraph::vertex_descriptor w = longestSubgraphPathVertices[i];
            const Segment segment = subgraph[w].segment;
            const double H = double(i) / (3. * double(longestSubgraphPathVertices.size() - 1));
            const double S = 1.;
            const double L = 0.6;
            const string color = hslToRgbString(H, S, L);
            csv << graph.assemblyGraph[segment].id << ",";
            csv << i << ",";
            csv << graph.assemblyGraph[segment].length() << ",";
            csv << color << ",";
            csv << "\n";
        }
    }

    // Construct the assembly path, without including segment0 and segment1.
    for(uint64_t i=1; i<longestSubgraphPathVertices.size()-1; i++) {
        const Subgraph::vertex_descriptor w = longestSubgraphPathVertices[i];
        const Segment segment = subgraph[w].segment;
        assemblyPath.push_back(segment);
    }

    return true;
}



uint64_t PathGraph::segmentId(vertex_descriptor v) const
{
    const PathGraph& pathGraph = *this;
    const Segment segment = pathGraph[v].segment;
    return graph.assemblyGraph[segment].id;
}



Subgraph::Subgraph(
    const Graph& graph,
    const std::set<Graph::vertex_descriptor, Graph::OrderById>& vertices) :
    graph(graph)
{
    Subgraph& subgraph = *this;

    // Create the vertices;
    std::map<Graph::vertex_descriptor, Subgraph::vertex_descriptor> vertexMap;
    for(const Graph::vertex_descriptor v: vertices) {
        vertex_descriptor w = boost::add_vertex(subgraph);
        SubgraphVertex& subgraphVertex = subgraph[w];
        subgraphVertex.segment = graph[v].segment;
        vertexMap.insert({v, w});
    }

    // Create the edges.
    for(const Graph::vertex_descriptor v0: vertices) {
        const Subgraph::vertex_descriptor w0 = vertexMap[v0];
        BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
            const Graph::vertex_descriptor v1 = target(e, graph);
            const auto it1 = vertexMap.find(v1);
            if(it1 != vertexMap.end()) {
                const Subgraph::vertex_descriptor w1 = it1->second;
                auto[we, ignore] = add_edge(w0, w1, subgraph);
                subgraph[we].correctedJaccard = graph[e].segmentPairInformation.correctedJaccard;
            }
        }
    }

}



void Subgraph::writeGraphviz(const string& name) const
{
    const Subgraph& subgraph = *this;

    ofstream dot(name);
    dot << "digraph ReadFollowingSubgraph {\n";

    BGL_FORALL_VERTICES(w, subgraph, Subgraph) {
        dot << segmentId(w) << ";\n";
    }

    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        const vertex_descriptor w0 = source(e, subgraph);
        const vertex_descriptor w1 = target(e, subgraph);
        dot <<
            segmentId(w0) << "->" <<
            segmentId(w1) << ";\n";
    }

    dot << "}\n";
}



uint64_t Subgraph::segmentId(vertex_descriptor w) const
{
    const Subgraph& subgraph = *this;
    const SubgraphVertex& vertex = subgraph[w];
    const Segment segment = vertex.segment;
    return graph.assemblyGraph[segment].id;
}



void Subgraph::makeAcyclic()
{
    Subgraph& subgraph = *this;

    // Gather the edges with their correctedJaccard.
    class EdgeInfo {
    public:
        edge_descriptor e;
        double correctedJaccard;
        bool operator<(const EdgeInfo& that) const {
            return correctedJaccard > that.correctedJaccard;
        }
    };
    vector<EdgeInfo> edgeInfos;
    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        const double correctedJaccard = subgraph[e].correctedJaccard;
        edgeInfos.push_back(EdgeInfo({e, correctedJaccard}));
    }

    // Sort them by decreasing corrected jaccard.
    sort(edgeInfos.begin(), edgeInfos.end());

    // Do approximateTopologicalOrdewr with edges in this order.
    vector<edge_descriptor> sortedEdges;
    for(const EdgeInfo& edgeInfo: edgeInfos) {
        sortedEdges.push_back(edgeInfo.e);
    }
    approximateTopologicalSort(subgraph, sortedEdges);

    // Remove the non-DAG edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        if(not subgraph[e].isDagEdge) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, subgraph);
    }

}
