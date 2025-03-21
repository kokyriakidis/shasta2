// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "findLinearChains.hpp"
#include "Tangle.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>

// Standard library.
#include <fstream.hpp>
#include <queue>
#include <set>



AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    const AssemblerOptions::AssemblyGraphOptions& options) :
    MappedMemoryOwner(anchors),
    anchors(anchors)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains in the AnchorGraph.
    vector< vector<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(anchorGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the anchor graph." << endl;

    // Get the initial and final AnchorId of all chains.
    // These generate an AssemblyGraph vertex each, after deduplication.
    vector<AnchorId> vertexAnchorIds;
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        vertexAnchorIds.push_back(anchorIdA);
        vertexAnchorIds.push_back(anchorIdB);
    }
    deduplicate(vertexAnchorIds);

    // Create the vertices.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const AnchorId anchorId: vertexAnchorIds) {
        const vertex_descriptor v = add_vertex(AssemblyGraphVertex(anchorId), assemblyGraph);
        vertexMap.insert(make_pair(anchorId, v));
    }



    // Create the edges.
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        const vertex_descriptor vA = vertexMap[anchorIdA];
        const vertex_descriptor vB = vertexMap[anchorIdB];

        edge_descriptor e;
        tie(e, ignore) = add_edge(vA, vB, assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.id = nextEdgeId++;

        for(const AnchorGraph::edge_descriptor e: chain) {
            edge.emplace_back(anchorGraph[e]);
        }
    }

    cout << "The initial assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("A");

    transitiveReduction(
        options.transitiveReductionThreshold,
        options.transitiveReductionA,
        options.transitiveReductionB);
    compress();

    cout << "After transitive reduction, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("B");

    TrivialDetangler detangler;

    detangleVertices(detangler);
    compress();

    cout << "After trivial vertex detangling, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("C");

    detangleEdges(detangler);
    compress();

    cout << "After trivial edge detangling, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("D");
}



// Compute the sum of offsets of all the AssemblyGraphSteps of an AssemblyGraphEdge.
void AssemblyGraphEdge::getOffsets(
    const Anchors& anchors,
    uint32_t& averageOffset,
    uint32_t& minOffset,
    uint32_t& maxOffset) const
{
    averageOffset = 0;
    minOffset = 0;
    maxOffset = 0;
    for(const AssemblyGraphStep& step: *this) {
        uint32_t stepAverageOffset;
        uint32_t stepMinOffset;
        uint32_t stepMaxOffset;
        step.getOffsets(anchors, stepAverageOffset, stepMinOffset, stepMaxOffset);
        averageOffset += stepAverageOffset;
        minOffset += stepMinOffset;
        maxOffset += stepMaxOffset;
    }
}



uint32_t AssemblyGraphEdge::getAverageOffset(const Anchors& anchors) const
{
    uint32_t averageOffset = 0;
    for(const AssemblyGraphStep& step: *this) {
        averageOffset += step.getAverageOffset(anchors);
    }
    return averageOffset;
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length in bases.
        gfa << "LN:i:" << edge.getAverageOffset(anchors) << "\n";
    }

    // For each vertex, write links between each pair of incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;
                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }
}



void AssemblyGraph::write(const string& name) const
{
    save(name);
    writeGfa("AssemblyGraph-" + name + ".gfa");
    writeGraphviz("AssemblyGraph-" + name + ".dot");
    writeSegments("AssemblyGraphSegments-" + name + ".csv");
    writeSegmentDetails("AssemblyGraphSegmentsDetails-" + name + ".csv");
}



void AssemblyGraph::writeSegments(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Id,AnchorIdA,AnchorIdB,StepCount,AverageOffset,MinOffset,MaxOffset\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        uint32_t averageOffset;
        uint32_t minOffset;
        uint32_t maxOffset;
        edge.getOffsets(anchors, averageOffset, minOffset, maxOffset);

        csv <<
            edge.id << "," <<
            anchorIdToString(edge.anchorIdA()) << "," <<
            anchorIdToString(edge.anchorIdB()) << "," <<
            edge.size() << "," <<
            averageOffset << "," <<
            minOffset << "," <<
            maxOffset << "\n";
    }

}


void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream dot(fileName);
    dot << "digraph AssemblyGraph {\n";

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        dot << "\"" << anchorIdToString(vertex.anchorId) << "\";\n";
    }

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
        const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];

        dot <<
            "\"" << anchorIdToString(vertex0.anchorId) << "\""
            "->"
            "\"" << anchorIdToString(vertex1.anchorId) << "\""
            " [label=\"" << edge.id << "\"]"
            ";\n";

    }

    dot << "}\n";
}



void AssemblyGraph::writeSegmentDetails(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Id,Position,AnchorIdA,AnchorIdB,AverageOffset,MinOffset,MaxOffset\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        for(uint64_t i=0; i<edge.size(); i++) {
            const AssemblyGraphStep& step = edge[i];

            uint32_t averageOffset;
            uint32_t minOffset;
            uint32_t maxOffset;
            step.getOffsets(anchors, averageOffset, minOffset, maxOffset);

            csv <<
                edge.id << "," <<
                i << "," <<
                anchorIdToString(step.anchorIdA) << "," <<
                anchorIdToString(step.anchorIdB) << "," <<
                averageOffset << "," <<
                minOffset << "," <<
                maxOffset << "\n";
        }
    }

}



// An v0->v1 edge with length length01 is removed if:
// - length01 is no more than threshold1.
// - A BFS that starts at v0 and does not use edge v0->v1
//   reaches v1 with a BFS offset no more than a  + b * length01.
// Most of the edges removed in this way are caused by errors.
// Even if a correct one is removed at this edge, it can still
// be assembled correctly by the LocalAssembly, is phasing is
// successful in this region.
void AssemblyGraph::transitiveReduction(
    uint64_t threshold,
    uint64_t a,
    uint64_t b)
{
    AssemblyGraph& assemblyGraph = *this;

    // Work areas for the BFS are defined here to reduce nemory allocation activity.
    std::queue<vertex_descriptor> q;
    vector<vertex_descriptor> verticesEncountered;

    // Loop over all edges. We will be removing edges while iterating,
    // so we have to avoid invalidating iterators.
    edge_iterator it, itEnd;
    tie(it, itEnd) = edges(assemblyGraph);
    while(it != itEnd) {
        const edge_descriptor e = *it;
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Increment the iterator here, before possibly removing this edge.
        ++it;

        if(edge.length(anchors) > threshold) {
            continue;
        }

        // This edge is a candidate for removal.
        // Initialize the BFS.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        SHASTA_ASSERT(q.empty());
        SHASTA_ASSERT(verticesEncountered.empty());
        q.push(v0);
        verticesEncountered.push_back(v0);
        assemblyGraph[v0].bfsDistance = 0;
        const uint64_t threshold2 = a + b * edge.length(anchors);

        // Main BFS loop.
        bool removeEdge = false;
        while(not q.empty()) {
            const vertex_descriptor vA = q.front();
            q.pop();
            const uint64_t distanceA = assemblyGraph[vA].bfsDistance;

            // Loop over out-edges of vA.
            BGL_FORALL_OUTEDGES(vA, eAB, assemblyGraph, AssemblyGraph) {

                // Don't use edge e in the BFS.
                if(eAB == e) {
                    continue;
                }

                const vertex_descriptor vB = target(eAB, assemblyGraph);

                // If we already encountered vB, do nothing.
                if(assemblyGraph[vB].bfsDistance != invalid<uint64_t>) {
                    continue;
                }

                // If we got too far, don't use this in the BFS.
                const uint64_t distanceB = distanceA + assemblyGraph[eAB].length(anchors);
                if(distanceB > threshold2) {
                    continue;
                }

                // If vB is v1, edge e should be removed.
                if(vB == v1) {
                    removeEdge = true;
                    break;
                }

                // Update the BFS.
                q.push(vB);
                verticesEncountered.push_back(vB);
                assemblyGraph[vB].bfsDistance = distanceB;

            }

            if(removeEdge) {
                break;
            }
        }

        if(removeEdge) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Cleanup the BFS.
        for(const vertex_descriptor v: verticesEncountered) {
            assemblyGraph[v].bfsDistance = invalid<uint64_t>;
        }
        while(not q.empty()) {
            q.pop();
        }
        verticesEncountered.clear();
    }
}



void AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains. All the edges in a linear chain of
    // more than one edge are combined into a single edge.
    vector< vector<edge_descriptor> > linearChains;
    findLinearChains(assemblyGraph, 0, linearChains);
    cout << "Compress found " << linearChains.size() << " linear chains." << endl;

    for(const vector<edge_descriptor>& linearChain: linearChains) {
        if(linearChain.size() == 1) {
            continue;
        }

        const vertex_descriptor v0 = source(linearChain.front(), assemblyGraph);
        const vertex_descriptor v1 = target(linearChain.back(), assemblyGraph);

        // Create a new edge consisting of just this linear chain.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, assemblyGraph);
        AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
        newEdge.id = nextEdgeId++;

        // Concatenate the old edges to create the new one,
        // then remove the old edges.
        for(const edge_descriptor eOld: linearChain) {
            AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
            copy(oldEdge.begin(), oldEdge.end(), back_inserter(newEdge));
            boost::remove_edge(eOld, assemblyGraph);
        }

        // Now we can remove the middle vertices in the linear chain.
        for(uint64_t i=1; i<linearChain.size(); i++) {
            const vertex_descriptor v = source(linearChain[i], assemblyGraph);
            remove_vertex(v, assemblyGraph);
        }
    }
}



void AssemblyGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AssemblyGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AssemblyGraph::save(const string& stage) const
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
    const string name = largeDataName("AssemblyGraph-" + stage);
    MemoryMapped::Vector<char> data;
    data.createNew(name, largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AssemblyGraph::load(const string& assemblyStage)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        const string name = largeDataName("AssemblyGraph-" + assemblyStage);
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error("Assembly graph at stage " + assemblyStage +
            " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error("Error reading assembly graph at stage " + assemblyStage +
            ": " + e.what());
    }
}



// Deserialize.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const string& assemblyStage) :
    MappedMemoryOwner(anchors),
    anchors(anchors)
{
    load(assemblyStage);
}



void AssemblyGraph::detangleVertices(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather all the vertices with in_degree and out_degree greater than 1.
    vector<vertex_descriptor> detanglingCandidates;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if(
            (in_degree(v, assemblyGraph) > 1) and
            (out_degree(v, assemblyGraph) > 1)) {
            detanglingCandidates.push_back(v);
        }
    }
    cout << "Found " << detanglingCandidates.size() <<
        " tangle vertices." << endl;

    uint64_t successCount = 0;
    for(const vertex_descriptor v: detanglingCandidates) {
        Tangle tangle(assemblyGraph, v);
        const bool success = detangler(tangle);
        if(success) {
            ++successCount;
        }
    }
    cout << "Detangling successful for " << successCount << " tangle vertices." << endl;
}



// When detangling edges, we must be careful because successful
// detangling operations can remove edges.
void AssemblyGraph::detangleEdges(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather all the edges with in_degree and out_degree greater than 1.
    std::set<edge_descriptor> detanglingCandidates;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        if(
            (in_degree(v0, assemblyGraph)  + in_degree(v1, assemblyGraph) > 1) and
            (out_degree(v0, assemblyGraph) + out_degree(v1, assemblyGraph) > 1)) {
            detanglingCandidates.insert(e);
        }
    }
    cout << "Found " << detanglingCandidates.size() <<
        " tangle edges out of " << num_edges(assemblyGraph) << " total edges." << endl;

    uint64_t attemptCount = 0;
    uint64_t successCount = 0;
    while(not detanglingCandidates.empty()) {
        auto it = detanglingCandidates.begin();
        const edge_descriptor e = *it;
        detanglingCandidates.erase(it);

        ++attemptCount;
        Tangle tangle(assemblyGraph, e);
        const bool success = detangler(tangle);
        if(success) {
            ++successCount;
            for(const edge_descriptor e: tangle.removedEdges) {
                detanglingCandidates.erase(e);
            }
        }
    }
    cout << "Attempted detangling for " << attemptCount << " tangle edges." << endl;
    cout << "Detangling was successful for " << successCount << " tangle edges." << endl;

}
