// Shasta2.
#include "AnchorSimilarityGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include <queue>



// Construct the AnchorSimilarityGraph from the completeAnchorGraph.
// Only include edges with at least the specified minCoverage.
AnchorSimilarityGraph::AnchorSimilarityGraph(
    const Anchors& anchors,
    const AnchorGraph& completeAnchorGraph) :
    MappedMemoryOwner(anchors)
{
    AnchorSimilarityGraph& anchorSimilarityGraph = *this;

    createVertices(anchors);
    createEdges(anchors, completeAnchorGraph);

    cout << "The anchor similarity graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Count entrances and exits.
    uint64_t entranceCount = 0;
    uint64_t exitCount = 0;
    BGL_FORALL_VERTICES(anchorId, anchorSimilarityGraph, AnchorSimilarityGraph) {
        if(in_degree(anchorId, anchorSimilarityGraph) == 0) {
            ++entranceCount;
        }
        if(out_degree(anchorId, anchorSimilarityGraph) == 0) {
            ++exitCount;
        }
    }
    cout << "The anchor similarity graph has " << entranceCount <<
        " entrances and " << exitCount << " exits." << endl;
}



// Constructor from binary data.
AnchorSimilarityGraph::AnchorSimilarityGraph(
    const MappedMemoryOwner& mappedMemoryOwner,
    const string& name) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    load(name);
}

void AnchorSimilarityGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AnchorSimilarityGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AnchorSimilarityGraph::save(const string& name) const
{
    // If not using persistent binary data, do nothing.
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    // First save to a string.
    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    // Now save the string to binary data.
    MemoryMapped::Vector<char> data;
    data.createNew(largeDataName(name), largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AnchorSimilarityGraph::load(const string& name)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        data.accessExistingReadOnly(largeDataName(name));
    } catch (std::exception&) {
        throw runtime_error(name + " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error(string("Error reading " + name + ": ") + e.what());
    }
}



// Create the vertices, one for each AnchorId.
// In the AnchorSimilarityGraph, vertex_descriptors are AnchorIds.
void AnchorSimilarityGraph::createVertices(const Anchors& anchors)
{
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(*this);
    }
}



// Create all the edges.
void AnchorSimilarityGraph::createEdges(
    const Anchors& anchors,
    const AnchorGraph& completeAnchorGraph)
{
    using Graph = AnchorSimilarityGraph;
    Graph& graph = *this;

    vector<uint8_t> color(anchors.size(), 0);
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
        createEdges(anchors, completeAnchorGraph, anchorId, color);
    }

    // Make sure that if an edge exists its reverse complement also exists.
    // It's not clear why this is necessary.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const AnchorId anchorIdA = source(e, graph);
        const AnchorId anchorIdB = target(e, graph);
        const AnchorId anchorIdARc = reverseComplementAnchorId(anchorIdA);
        const AnchorId anchorIdBRc = reverseComplementAnchorId(anchorIdB);
        auto [eRc, edgeExists] = boost::edge(anchorIdBRc, anchorIdARc, graph);
        if(not edgeExists) {
            add_edge(anchorIdBRc, anchorIdARc, graph[e], graph);
        }
    }
}



// Create the edges with source anchorIdA.
// This uses a forward BFS in the completeAnchorGraph starting at anchorIdA.
void AnchorSimilarityGraph::createEdges(
    const Anchors& anchors,
    const AnchorGraph& completeAnchorGraph,
    AnchorId anchorIdA,
    vector<uint8_t>& color)
{
    AnchorSimilarityGraph& anchorSimilarityGraph = *this;
    const bool debug = false;
    if(debug) {
        cout << "Creating AnchorSimilarityGraph edges with source " << anchorIdToString(anchorIdA) << endl;
    }

    // Initialize the BFS.
    color[anchorIdA] = 1;
    vector<AnchorId> visited;       // Pass as argument instead to reduce memory allocation activity.
    visited.push_back(anchorIdA);
    std::queue<AnchorId> q;         // Pass as argument instead to reduce memory allocation activity.
    q.push(anchorIdA);



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue an AnchorId.
        const AnchorId anchorId0 = q.front();
        q.pop();
        if(debug) {
            cout << "Dequeued " << anchorIdToString(anchorId0) << endl;
        }

        // Loop over its out-edges in the completeAnchorGraph.
        BGL_FORALL_OUTEDGES(anchorId0, e, completeAnchorGraph, AnchorGraph) {
            const AnchorId anchorId1 = target(e, completeAnchorGraph);
            if(debug) {
                cout << "Found " << anchorIdToString(anchorId1) << endl;
            }

            // If we already visited anchorId1, skip it.
            if(color[anchorId1] == 1) {
                continue;
            }

            // Mark it as visited.
            color[anchorId1] = 1;
            visited.push_back(anchorId1);

            // If no common reads, skip it entirely.
            auto [commonCount, nonPositiveOffsetFound] = anchors.countCommonWithFlag(anchorIdA, anchorId1);
            if(commonCount == 0) {
                continue;
            }

            // Enqueue it.
            q.push(anchorId1);

            // If it satisfies our requirements, add an edge anchorIdA->anchorId1
            // to the AnchorSimilarityGraph.
            if((commonCount >=  minCommonCount) and (not nonPositiveOffsetFound)) {
                AnchorPairInfo info;
                anchors.analyzeAnchorPair(anchorIdA, anchorId1, info);
                SHASTA2_ASSERT(info.commonPositiveOffset == commonCount);
                SHASTA2_ASSERT(info.commonNonPositiveOffset == 0);
                SHASTA2_ASSERT(info.offsetInBases < 1000000000);
                const uint64_t missingCount = info.missingCount();
                const double logP = a * double(commonCount) - b * double(missingCount);
                if(logP >= minLogP) {
                    const double weight = std::pow(10., -0.1 * logP);
                    if(debug) {
                        cout << "Added edge " << anchorIdToString(anchorIdA) << " " << anchorIdToString(anchorId1) << " " << info.offsetInBases << endl;
                    }
                    auto [ignore, edgeExists] = boost::edge(anchorIdA, anchorId1, anchorSimilarityGraph);
                    SHASTA2_ASSERT(not edgeExists);
                    add_edge(anchorIdA, anchorId1, AnchorSimilarityGraphEdge(weight, info.offsetInBases), anchorSimilarityGraph);
                }
            }
         }
    }


    // Reset the colors.
    for(const AnchorId anchorId: visited) {
        color[anchorId] = 0;
    }
    visited.clear();
}




// Compute a shortest path tree starting at teh given AnchorId.
void AnchorSimilarityGraph::shortestPaths(AnchorId anchorIdA) const
{
    using namespace boost;

    using Graph = AnchorSimilarityGraph;
    const Graph& graph = *this;

    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
    std::map<vertex_descriptor, double> distanceMap;
    std::map<vertex_descriptor, uint64_t> linearDistanceMap;
    std::map<vertex_descriptor, uint64_t> baseDistanceMap;
    linearDistanceMap.insert(make_pair(anchorIdA, 0));
    baseDistanceMap.insert(make_pair(anchorIdA, 0));

    class DijkstraVisitor : public boost::dijkstra_visitor<> {
    public:
        void examine_vertex(AnchorId anchorId, const Graph& graph)
        {
            const AnchorId predecessor = predecessorMapPointer->at(anchorId);
            if(predecessor == anchorId) {
                return;
            }
            linearDistanceMapPointer->insert(make_pair(anchorId, linearDistanceMapPointer->at(predecessor) + 1));
            auto[e, edgeExists] = boost::edge(predecessor, anchorId, graph);
            SHASTA2_ASSERT(edgeExists);
            const uint64_t baseOffset = graph[e].baseOffset;
            SHASTA2_ASSERT(baseOffset < 1000000000UL);
            baseDistanceMapPointer->insert(make_pair(anchorId, baseDistanceMapPointer->at(predecessor) + baseOffset));
        }
        std::map<vertex_descriptor, vertex_descriptor>* predecessorMapPointer;
        std::map<vertex_descriptor, double>* distanceMapPointer;
        std::map<vertex_descriptor, uint64_t>* linearDistanceMapPointer;
        std::map<vertex_descriptor, uint64_t>* baseDistanceMapPointer;
    };
    DijkstraVisitor dijkstraVisitor;
    dijkstraVisitor.predecessorMapPointer = &predecessorMap;
    dijkstraVisitor.distanceMapPointer = &distanceMap;
    dijkstraVisitor.linearDistanceMapPointer = &linearDistanceMap;
    dijkstraVisitor.baseDistanceMapPointer = &baseDistanceMap;


    cout << timestamp << "dijkstra_shortest_paths begins.\n";
    dijkstra_shortest_paths(graph, anchorIdA,
       weight_map(get(&AnchorSimilarityGraphEdge::weight, graph)).
       predecessor_map(make_assoc_property_map(predecessorMap)).
       distance_map(make_assoc_property_map(distanceMap)).
       visitor(dijkstraVisitor));
    cout << timestamp << "dijkstra_shortest_paths ends." << endl;

    {
        ofstream csv("DistanceMap.csv");
        csv << "AnchorId,Linear distance,Weight distance,Base distance,Weight per base,Predecessor,\n";
        for(const auto& [anchorId, weightDistance]: distanceMap) {
            if(weightDistance != std::numeric_limits<double>::max()) {
                const uint64_t linearDistance = linearDistanceMap.at(anchorId);
                const uint64_t baseDistance = baseDistanceMap.at(anchorId);
                const double weightPerBase = weightDistance / double(baseDistance);
                csv << anchorIdToString(anchorId) << "," <<
                    linearDistance << "," <<
                    weightDistance << "," <<
                    baseDistance << "," <<
                    weightPerBase << "," <<
                    anchorIdToString(predecessorMap.at(anchorId)) << "," <<
                    "\n";
            }
        }
    }

    cout << "AnchorSimilarityGraph::shortestPaths ends." << endl;
}



void AnchorSimilarityGraph::shortestPathsFast(
    AnchorId anchorId,
    vector<AnchorId>& predecessorMap,
    vector<double>& distanceMap,
    vector<boost::default_color_type>& colorMap,
    vector<AnchorId>& accessibleVertices
    ) const
{
    using namespace boost;
    const AnchorSimilarityGraph& graph = *this;
    accessibleVertices.clear();

    class DijkstraVisitor : public dijkstra_visitor<> {
    public:
        DijkstraVisitor(
            vector<AnchorId>& accessibleVertices,
            vector<AnchorId>& predecessorMap) :
            accessibleVerticesPointer(&accessibleVertices),
            predecessorMapPointer(&predecessorMap)
        {}
        void examine_vertex(AnchorId anchorId, const AnchorSimilarityGraph&)
        {
            accessibleVerticesPointer->push_back(anchorId);
        }
        /*
        void examine_edge(edge_descriptor e, const AnchorSimilarityGraph& graph)
        {
            const AnchorId anchorId0 = source(e, graph);
            const AnchorId anchorId1 = target(e, graph);
            cout << "examine_edge " << anchorIdToString(anchorId0) << " " << anchorIdToString(anchorId1) << endl;
        }
        */
        vector<AnchorId>* accessibleVerticesPointer;
        vector<AnchorId>* predecessorMapPointer;
    };
    DijkstraVisitor dijkstraVisitor(accessibleVertices, predecessorMap);

    // Create the shortest path tree starting at anchorId.
    // THE DOCUMENTATION for dijkstra_shortest_paths_no_init
    // FORGETS TO MENTION THAT YOU HAVE TO SET TO 0
    // THE DISTANCE OF THE STARTING VERTEX IN THE distanceMap.
    distanceMap[anchorId] = 0;
    dijkstra_shortest_paths_no_init(
        graph,
        anchorId,
        make_iterator_property_map(predecessorMap.begin(),  boost::get(vertex_index, graph)),
        make_iterator_property_map(distanceMap.begin(),  boost::get(vertex_index, graph)),
        boost::get(&AnchorSimilarityGraphEdge::weight, graph),
        boost::get(vertex_index, graph),
        std::less<double>(),
        std::plus<double>(),
        0.,
        dijkstraVisitor,
        make_iterator_property_map(colorMap.begin(), boost::get(vertex_index, graph))
        );

}



void AnchorSimilarityGraph::shortestPathsFast(
    AnchorId anchorId,
    const Anchors& anchors) const
{
    // Create the work areas required by the low level function.
    vector<AnchorId> predecessorMap(anchors.size());
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
        predecessorMap[anchorId] = anchorId;
    }
    vector<double> distanceMap(anchors.size(), std::numeric_limits<double>::max());
    vector<boost::default_color_type> colorMap(
        anchors.size(),
        boost::default_color_type::white_color);

    vector<AnchorId> accessibleVertices;

    const auto t0 = steady_clock::now();
    shortestPathsFast(anchorId, predecessorMap, distanceMap, colorMap, accessibleVertices);
    const auto t1 = steady_clock::now();
    cout << "Shortest path tree computation took " << seconds(t1-t0) << " s." << endl;

    // Reset the work areas.
    for(const AnchorId anchorId: accessibleVertices) {
        predecessorMap[anchorId] = anchorId;
        distanceMap[anchorId] = std::numeric_limits<double>::max();
        colorMap[anchorId] = boost::default_color_type::white_color;
    }

#if 0
    // Check that the work areas were reset correctly.
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
        SHASTA2_ASSERT(predecessorMap[anchorId] == anchorId);
        SHASTA2_ASSERT(distanceMap[anchorId] == std::numeric_limits<double>::max());
        SHASTA2_ASSERT(colorMap[anchorId] == boost::default_color_type::white_color);
    }
#endif
}



void AnchorSimilarityGraph::allShortestPaths(const Anchors& anchors)
{
    using Graph = AnchorSimilarityGraph;
    Graph& graph = *this;

    const bool debug = false;

    // Make sure all edges are not flagged as shortest path edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        graph[e].isShortestPathEdge = false;
    }

    // Create the work areas required for shortest path trees.
    vector<AnchorId> predecessorMap(anchors.size());
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
        predecessorMap[anchorId] = anchorId;
    }
    vector<double> distanceMap(anchors.size(), std::numeric_limits<double>::max());
    vector<boost::default_color_type> colorMap(
        anchors.size(),
        boost::default_color_type::white_color);

    vector<AnchorId> accessibleVertices;



    // Loop over all entrances.
    BGL_FORALL_VERTICES(anchorIdA, graph, Graph) {
        if((anchorIdA % 1000) == 0) {
            cout << timestamp << anchorIdA << "/" << anchors.size() << endl;
        }
        if(in_degree(anchorIdA, graph) > 0) {
            continue;
        }

        // Check that the work areas are set correctly.
        // Remove this when done debugging.
        for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
            SHASTA2_ASSERT(predecessorMap[anchorId] == anchorId);
            SHASTA2_ASSERT(distanceMap[anchorId] == std::numeric_limits<double>::max());
            SHASTA2_ASSERT(colorMap[anchorId] == boost::default_color_type::white_color);
        }

        // Find the shortest path tree.
        shortestPathsFast(anchorIdA, predecessorMap, distanceMap, colorMap, accessibleVertices);

        // DEBUG: CHECK THE predecessorMap.
        {
            ofstream csv("PredecessorMap.csv");
            for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
                csv << anchorIdToString(anchorId) << ",";
                csv << anchorIdToString(predecessorMap[anchorId]) << ",\n";
            }
        }

        // Flag all edges of the shortest path.
        for(const AnchorId anchorIdB: accessibleVertices) {
            if(anchorIdB == anchorIdA) {
                continue;
            }
            const AnchorId anchorIdC = predecessorMap[anchorIdB];
            if(debug) {
                cout << "Predecessor of  " << anchorIdToString(anchorIdB) <<
                    " is " << anchorIdToString(anchorIdC) << endl;
            }
            SHASTA2_ASSERT(anchorIdC != anchorIdB);

            auto [e, edgeExists] = boost::edge(anchorIdC, anchorIdB, graph);
            SHASTA2_ASSERT(edgeExists);
            graph[e].isShortestPathEdge = true;
        }

        // Reset the work areas.
        for(const AnchorId anchorId: accessibleVertices) {
            predecessorMap[anchorId] = anchorId;
            distanceMap[anchorId] = std::numeric_limits<double>::max();
            colorMap[anchorId] = boost::default_color_type::white_color;
        }

    }

    // Also flag as shortest path edges the reverse complements
    // of the edges already flagged as shortest path edges.
    vector<edge_descriptor> shortestPathEdges;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].isShortestPathEdge) {
            shortestPathEdges.push_back(e);
        }
    }
    for(const edge_descriptor e: shortestPathEdges) {
        const AnchorId anchorIdA = source(e, graph);
        const AnchorId anchorIdB = target(e, graph);
        const AnchorId anchorIdARc = reverseComplementAnchorId(anchorIdA);
        const AnchorId anchorIdBRc = reverseComplementAnchorId(anchorIdB);
        auto [eRc, edgeExists] = boost::edge(anchorIdBRc, anchorIdARc, graph);
        if(not edgeExists) {
            cout << "Assertion failed caused by:" << endl;
            cout << anchorIdToString(anchorIdA) << endl;
            cout << anchorIdToString(anchorIdB) << endl;
            cout << anchorIdToString(anchorIdARc) << endl;
            cout << anchorIdToString(anchorIdBRc) << endl;
        }
        SHASTA2_ASSERT(edgeExists);
        graph[eRc].isShortestPathEdge = true;
    }



    // Count the edges we flagged.
    uint64_t shortesPathEdgeCount = 0;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].isShortestPathEdge) {
            ++shortesPathEdgeCount;
        }
    }
    cout << "Flagged " << shortesPathEdgeCount <<
        " edges as shortest path edges out of " << num_edges(graph) << " total." << endl;

    writeGraphviz("AnchorSimilarityGraph.dot");

}



// Graphviz output only includes the edges flgged as shortest path edges.
void AnchorSimilarityGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AnchorSimilarityGraph::writeGraphviz(ostream& dot) const
{
    using Graph = AnchorSimilarityGraph;
    const Graph& graph = *this;

    dot << "digraph AnchorSimilarityGraph { \n";

    BGL_FORALL_VERTICES(anchorId, graph, Graph) {
        dot << "\"" << anchorIdToString(anchorId) << "\"";

        if((in_degree(anchorId, graph) == 0) or (out_degree(anchorId, graph) == 0)) {
            dot << "[style=filled color=red]";
        }

        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        if(graph[e].isShortestPathEdge) {
            const AnchorId anchorId0 = source(e, graph);
            const AnchorId anchorId1 = target(e, graph);
            dot <<
                "\"" << anchorIdToString(anchorId0) << "\"->\"" <<
                anchorIdToString(anchorId1) << "\";\n";
        }
    }
    dot << "}\n";
}
