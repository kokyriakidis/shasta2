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
#include "fstream.hpp"
#include <queue>



// Construct the AnchorSimilarityGraph from the completeAnchorGraph.
// Only include edges with at least the specified minCoverage.
AnchorSimilarityGraph::AnchorSimilarityGraph(
    const Anchors& anchors,
    const AnchorGraph& completeAnchorGraph) :
    MappedMemoryOwner(anchors)
{
    createVertices(anchors);
    createEdges(anchors, completeAnchorGraph);

    cout << "The anchor similarity graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
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
    vector<uint8_t> color(anchors.size(), 0);
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId++) {
        createEdges(anchors, completeAnchorGraph, anchorId, color);
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
                        cout << anchorIdToString(anchorIdA) << " " << anchorIdToString(anchorId1) << " " << info.offsetInBases << endl;
                    }
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

    cout << "AnchorSimilarityGraph::shortestPaths begins for " <<
        anchorIdToString(anchorIdA) << endl;

    cout << "The anchor similarity graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

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



    dijkstra_shortest_paths(graph, anchorIdA,
       weight_map(get(&AnchorSimilarityGraphEdge::weight, graph)).
       predecessor_map(make_assoc_property_map(predecessorMap)).
       distance_map(make_assoc_property_map(distanceMap)).
       visitor(dijkstraVisitor));

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
