// Shasta2.
#include "AnchorSimilarityGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "LocalAssembly4.hpp"
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



#if 0
// Compute a shortest path tree starting at the given AnchorId.
void AnchorSimilarityGraph::createShortestPathTree(AnchorId anchorIdA) const
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
#endif



void AnchorSimilarityGraph::createShortestPathTree(
    AnchorId anchorId,
    ShortestPathTreeWorkAreas& workAreas) const
{
    using namespace boost;
    const AnchorSimilarityGraph& graph = *this;

    // THIS CHECK CAN BE EXPENSIVE. ONLY USE FOR DEBUGGING.
    // workAreas.check();


    // The DijkstraVisitor records the accessible vertices.
    class DijkstraVisitor : public dijkstra_visitor<> {
    public:
        DijkstraVisitor(vector<AnchorId>& accessibleVertices) :
            accessibleVertices(accessibleVertices)
        {}
        void discover_vertex(AnchorId anchorId, const AnchorSimilarityGraph&)
        {
            accessibleVertices.push_back(anchorId);
        }
        vector<AnchorId>& accessibleVertices;
    };
    DijkstraVisitor dijkstraVisitor(workAreas.accessibleVertices);

    // Create the shortest path tree starting at anchorId.
    // THE DOCUMENTATION for dijkstra_shortest_paths_no_init
    // FORGETS TO MENTION THAT YOU HAVE TO SET TO 0
    // THE DISTANCE OF THE STARTING VERTEX IN THE distanceMap.
    workAreas.data[anchorId].distance = 0;
    workAreas.accessibleVertices.push_back(anchorId);
    dijkstra_shortest_paths_no_init(
        graph,
        anchorId,
        workAreas.predecessorMap,
        workAreas.distanceMap,
        boost::get(&AnchorSimilarityGraphEdge::weight, graph),
        boost::get(vertex_index, graph),
        std::less<double>(),
        std::plus<double>(),
        0.,
        dijkstraVisitor,
        workAreas.colorMap
        );

}


AnchorSimilarityGraph::ShortestPathTreeWorkAreas::ShortestPathTreeWorkAreas(uint64_t anchorCount) :
    data(anchorCount)
{
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        data[anchorId].predecessor = anchorId;
    }
}



void AnchorSimilarityGraph::ShortestPathTreeWorkAreas::reset()
{
    for(const AnchorId anchorId: accessibleVertices) {
        Data& d = data[anchorId];
        d.predecessor = anchorId;
        d.distance = std::numeric_limits<double>::max();
        d.color = boost::default_color_type::white_color;
    }
    accessibleVertices.clear();
}



void AnchorSimilarityGraph::ShortestPathTreeWorkAreas::check() const
{
    for(AnchorId anchorId=0; anchorId<data.size(); anchorId++) {
        const Data& d = data[anchorId];
        SHASTA2_ASSERT(d.predecessor == anchorId);
        SHASTA2_ASSERT(d.distance == std::numeric_limits<double>::max());
        SHASTA2_ASSERT(d.color == boost::default_color_type::white_color);
    }
}



void AnchorSimilarityGraph::createShortestPathTree(
    AnchorId anchorId,
    const Anchors& anchors) const
{
    ShortestPathTreeWorkAreas workAreas(anchors.size());
    createShortestPathTree(anchorId, workAreas);
    const ShortestPathTree tree(*this, anchorId, workAreas);
    workAreas.reset();

    cout << "The shortest path tree has " << num_vertices(tree) <<
        " vertices and " << num_edges(tree) << " edges." << endl;
    cout << "The maximum path length in this tree is " << tree.maximumPathLength() << endl;
    SHASTA2_ASSERT(num_edges(tree) == num_vertices(tree) - 1);
    tree.writeGraphviz("ShortestPathTree.dot");

    // Interactive loop.
    while(true) {
        cout << "Enter final AnchorId for assembly path:" << endl;
        string endAnchorIdString;
        cin >> endAnchorIdString;
        const AnchorId endAnchorId = anchorIdFromString(endAnchorIdString);

        vector<AnchorId> path;
        tree.findPathAnchorIds(endAnchorId, path);
        cout << "Path ending at " << anchorIdToString(endAnchorId) << ":";
        for(const AnchorId anchorId: path) {
            cout << " " << anchorIdToString(anchorId);
        }
        cout << endl;

        // Assemble the path.
        vector<shasta2::Base> sequence;
        for(uint64_t i1=1; i1<path.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const AnchorId anchorId0 = path[i0];
            const AnchorId anchorId1 = path[i1];
            const AnchorPair anchorPair(anchors, anchorId0, anchorId1, false);
            ostream html(0);
            vector<OrientedReadId> additionalOrientedReadIds;
            const LocalAssembly4 localAssembly(anchors, 5000, html, false, anchorPair, additionalOrientedReadIds);

            std::ranges::copy(localAssembly.sequence, back_inserter(sequence));
        }
        ofstream fasta("ShortestPathTree.fasta");
        fasta << ">ShortestPathTree\n";
        std::ranges::copy(sequence, ostream_iterator<shasta2::Base>(fasta));
    }


#if 0
    // Do a local assembly for each edge of the tree.
    ofstream fasta("ShortestPathTree.fasta");
    BGL_FORALL_EDGES(e, tree, ShortestPathTree) {
        const ShortestPathTree::vertex_descriptor v0 = source(e, tree);
        const ShortestPathTree::vertex_descriptor v1 = target(e, tree);
        const AnchorId anchorId0 = tree[v0].anchorId;
        const AnchorId anchorId1 = tree[v1].anchorId;
        AnchorPair anchorPair(anchors, anchorId0, anchorId1, false);

        cout << anchorIdToString(anchorId0) << " " << anchorIdToString(anchorId1) <<
            " " << anchorPair.orientedReadIds.size() << endl;

        ostream html(0);
        vector<OrientedReadId> additionalOrientedReadIds;
        const LocalAssembly4 localAssembly(anchors, 5000, html, false, anchorPair, additionalOrientedReadIds);

        fasta << ">" << anchorIdToString(anchorId0) << "_" << anchorIdToString(anchorId1) << "\n";
        std::ranges::copy(localAssembly.sequence, ostream_iterator<shasta2::Base>(fasta));
        fasta << "\n";

    }
#endif
}



void AnchorSimilarityGraph::flagShortestPathEdges(const Anchors& anchors)
{
    using Graph = AnchorSimilarityGraph;
    Graph& graph = *this;

    const bool debug = false;

    // Make sure all edges are not flagged as shortest path edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        graph[e].isShortestPathEdge = false;
    }

    // Create the work areas required for shortest path trees.
    ShortestPathTreeWorkAreas workAreas(anchors.size());



    // Loop over all entrances.
    BGL_FORALL_VERTICES(anchorIdA, graph, Graph) {
        if(in_degree(anchorIdA, graph) > 0) {
            continue;
        }

        // Find the shortest path tree.
        createShortestPathTree(anchorIdA, workAreas);

        // Flag all edges of the shortest path.
        for(const AnchorId anchorIdB: workAreas.accessibleVertices) {
            if(anchorIdB == anchorIdA) {
                continue;
            }
            const AnchorId anchorIdC = workAreas.data[anchorIdB].predecessor;
            if(debug) {
                cout << "Predecessor of  " << anchorIdToString(anchorIdB) <<
                    " is " << anchorIdToString(anchorIdC) << endl;
            }
            if(anchorIdC != anchorIdB) {
                auto [e, edgeExists] = boost::edge(anchorIdC, anchorIdB, graph);
                SHASTA2_ASSERT(edgeExists);
                graph[e].isShortestPathEdge = true;
            }
        }

        // Reset the work areas.
        workAreas.reset();
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



AnchorSimilarityGraph::ShortestPathTree::ShortestPathTree(
    const AnchorSimilarityGraph& graph,
    AnchorId rootAnchorId,
    const ShortestPathTreeWorkAreas& workAreas) :
    rootAnchorId(rootAnchorId)
{
    using Tree = ShortestPathTree;
    Tree& tree = *this;

    // Add the vertices.
    for(const AnchorId anchorId: workAreas.accessibleVertices) {
        if(not vertexMap.contains(anchorId)) {
            const vertex_descriptor v = add_vertex(ShortestPathTreeVertex(anchorId), tree);
            vertexMap.insert(make_pair(anchorId, v));
        }
    }

    // Add the edges.
    BGL_FORALL_VERTICES(v1, tree, Tree) {
        const AnchorId anchorId1 = tree[v1].anchorId;
        const AnchorId anchorId0 = workAreas.data[anchorId1].predecessor;
        if(anchorId0 != anchorId1) {
            const vertex_descriptor v0 = vertexMap.at(anchorId0);

            // Locate the corresponding edge in the AnchorSimilarityGraph.
            auto [e, edgeExists] = boost::edge(anchorId0, anchorId1, graph);
            SHASTA2_ASSERT(edgeExists);

            // Extract the logP.
            const double weight = graph[e].weight;
            const double logP = -10. * std::log10(weight);

            // Add the edge to the tree.
            add_edge(v0, v1, ShortestPathTreeEdge(logP), tree);
        }
    }

    computeDistancesToRoot();
    computeLongestDistancesToLeaf();
}



void AnchorSimilarityGraph::ShortestPathTree::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AnchorSimilarityGraph::ShortestPathTree::writeGraphviz(ostream& dot) const
{
    using Tree = ShortestPathTree;
    const Tree& tree = *this;

    dot << "digraph ShortestPathTree {\n";

    BGL_FORALL_VERTICES(v, tree, Tree) {
        dot << "\"" << anchorIdToString(tree[v].anchorId) << "\"";
        dot << "[label=\"" <<
            anchorIdToString(tree[v].anchorId) <<
            "\\n" << tree[v].distanceToRoot <<
            "\\n" << tree[v].longestDistanceToLeaf <<
            "\"]";
        dot << ";\n";
    }

    dot << std::fixed << std::setprecision(1);
    BGL_FORALL_EDGES(e, tree, Tree) {
        const vertex_descriptor v0 = source(e, tree);
        const vertex_descriptor v1 = target(e, tree);
        dot << "\"" << anchorIdToString(tree[v0].anchorId) << "\"->\"";
        dot << anchorIdToString(tree[v1].anchorId) << "\"";
        dot << "[label=\"" << tree[e].logP << "\"]";
        dot << ";\n";
    }

    dot << "}\n";

}



void AnchorSimilarityGraph::ShortestPathTree::computeDistancesToRoot()
{
    using Tree = ShortestPathTree;
    Tree& tree = *this;

    const vertex_descriptor root = vertexMap.at(rootAnchorId);
    tree[root].distanceToRoot = 0;

    std::queue<vertex_descriptor> q;
    q.push(root);

    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        const uint64_t distance0 = tree[v0].distanceToRoot;
        const uint64_t distance1 = distance0 + 1;

        BGL_FORALL_OUTEDGES(v0, e, tree, Tree) {
            const vertex_descriptor v1 = target(e, tree);
            q.push(v1);
            tree[v1].distanceToRoot = distance1;
        }
    }
}



uint64_t AnchorSimilarityGraph::ShortestPathTree::maximumPathLength() const
{
    using Tree = ShortestPathTree;
    const Tree& tree = *this;

    uint64_t length = 0;

    BGL_FORALL_VERTICES(v, tree, Tree) {
        length = max(length, tree[v].distanceToRoot);
    }
    return length;
}



void AnchorSimilarityGraph::ShortestPathTree::computeLongestDistancesToLeaf()
{
    using Tree = ShortestPathTree;
    Tree& tree = *this;

    // First set to 0 the longestDistanceToLeaf for all leafs.
    BGL_FORALL_VERTICES(v, tree, Tree) {
        if(out_degree(v, tree) == 0) {
            tree[v].longestDistanceToLeaf = 0;
        }
    }

    // Then loop over non-leaf vertices in order of decreasing distance from root.
    vector < vector<vertex_descriptor> > verticesByDistanceFromRoot;
    BGL_FORALL_VERTICES(v, tree, Tree) {
        const uint64_t distanceFromRoot = tree[v].distanceToRoot;
        if(verticesByDistanceFromRoot.size() <= distanceFromRoot) {
            verticesByDistanceFromRoot.resize(distanceFromRoot + 1);
        }
        verticesByDistanceFromRoot[distanceFromRoot].push_back(v);
    }
    for(auto it=verticesByDistanceFromRoot.rbegin(); it!=verticesByDistanceFromRoot.rend(); ++it) {
        const vector<vertex_descriptor>& verticesAtThisDistance = *it;
        for(const vertex_descriptor v0: verticesAtThisDistance) {
            if(out_degree(v0, tree) == 0) {
                continue;
            }

            // The longest distance to leaf of v0 is the maximum of
            // the longest distance to leaf of its children, plus 1.
            uint64_t maxChildrenDistance = 0;
            BGL_FORALL_OUTEDGES(v0, e, tree, Tree) {
                const vertex_descriptor v1 = target(e, tree);
                const uint64_t childDistance = tree[v1].longestDistanceToLeaf;
                SHASTA2_ASSERT(childDistance != invalid<uint64_t>);
                maxChildrenDistance = max(maxChildrenDistance, childDistance);
            }
            tree[v0].longestDistanceToLeaf = maxChildrenDistance + 1;
        }
    }
}



// Find the sequence of vertices of a path starting at root
// and ending at the specified vertex.
void AnchorSimilarityGraph::ShortestPathTree::findPath(
    vertex_descriptor v,
    vector<vertex_descriptor>& path) const
{
    using Tree = ShortestPathTree;
    const Tree& tree = *this;

    const vertex_descriptor root = vertexMap.at(rootAnchorId);

    path.clear();
    while(true) {
        path.push_back(v);
        if(v == root) {
            break;
        }

        auto [it, ignore] = in_edges(v, tree);
        const edge_descriptor e = *it;
        v = source(e, tree);
    }

    std::ranges::reverse(path);
}



// Find the sequence of AnchorIds of a path starting at root
// and ending at the specified AnchorId.
void AnchorSimilarityGraph::ShortestPathTree::findPathAnchorIds(
    AnchorId anchorId,
    vector<AnchorId>& anchorIds) const
{
    using Tree = ShortestPathTree;
    const Tree& tree = *this;

    vector<vertex_descriptor> path;
    findPath(vertexMap.at(anchorId), path);

    anchorIds.clear();
    for(const vertex_descriptor v: path) {
        anchorIds.push_back(tree[v].anchorId);
    }
}
