// Shasta2.
#include "AnchorSimilarityGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>

// Standard library.
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
        cout << "Creating AnchorSimilarityGraph with source " << anchorIdToString(anchorIdA) << endl;
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

            // Check anchorId1 against anchorIdA.
            AnchorPairInfo info;
            anchors.analyzeAnchorPair(anchorIdA, anchorId1, info);

            // If no common reads, skip it entirely.
            if(info.common == 0) {
                continue;
            }

            // Enqueue it.
            q.push(anchorId1);

            // If it satisfies our requirements, add an edge anchorIdA->anchorId1
            // to the AnchorSimilarityGraph.
            const uint64_t common = info.common;
            const uint64_t missing = info.missingCount();
            const double logP = a * double(common) - b * double(missing);
            if((common >= minCommonCount) and (logP >= minLogP)) {
                const double weight = std::pow(10., -0.1 * logP);
                add_edge(anchorIdA, anchorId1, AnchorSimilarityGraphEdge(weight), anchorSimilarityGraph);
            }
         }
    }


    // Reset the colors.
    for(const AnchorId anchorId: visited) {
        color[anchorId] = 0;
    }
    visited.clear();
}

