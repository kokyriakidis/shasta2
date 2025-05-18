// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "findLinearChains.hpp"
#include "html.hpp"
#include "inducedSubgraphIsomorphisms.hpp"
#include "LocalAssembly.hpp"
#include "LocalAssembly1.hpp"
#include "LocalAssembly2.hpp"
#include "performanceLog.hpp"
#include "PermutationDetangler.hpp"
#include "Tangle.hpp"
#include "TrivialDetangler.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>

// Standard library.
#include <chrono.hpp>
#include <cstdlib>
#include <fstream.hpp>
#include <queue>
#include <set>

// Explicit instantiationn.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph>;



AssemblyGraph::AssemblyGraph(
    const AssemblerOptions& assemblerOptions,
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    uint64_t threadCount) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    assemblerOptions(assemblerOptions),
    anchors(anchors)
{
    AssemblyGraph& assemblyGraph = *this;
    const auto& assemblyGraphOptions = assemblerOptions.assemblyGraphOptions;
    performanceLog << "AssemblyGraph creation begins." << endl;


    // Find linear chains in the AnchorGraph.
    vector< vector<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(anchorGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the anchor graph." << endl;

    // Get the initial and final AnchorId of all chains.
    // These generate an AssemblyGraph vertex each, after deduplication.
    performanceLog << timestamp << "Finding linear chains in the AnchorGraph." << endl;
    vector<AnchorId> vertexAnchorIds;
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorPair.anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorPair.anchorIdB;
        vertexAnchorIds.push_back(anchorIdA);
        vertexAnchorIds.push_back(anchorIdB);
    }
    deduplicate(vertexAnchorIds);

    // Create the vertices.
    performanceLog << timestamp << "Creating AssemblyGraph vertices." << endl;
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const AnchorId anchorId: vertexAnchorIds) {
        const vertex_descriptor v = add_vertex(AssemblyGraphVertex(anchorId), assemblyGraph);
        vertexMap.insert(make_pair(anchorId, v));
    }



    // Create the edges.
    performanceLog << timestamp << "Creating AssemblyGraph edges." << endl;
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorPair.anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorPair.anchorIdB;
        const vertex_descriptor vA = vertexMap[anchorIdA];
        const vertex_descriptor vB = vertexMap[anchorIdB];

        edge_descriptor e;
        tie(e, ignore) = add_edge(vA, vB, assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.id = nextEdgeId++;

        edge.push_back(anchorIdA);
        for(const AnchorGraph::edge_descriptor e: chain) {
            const AnchorGraph::vertex_descriptor v = target(e, anchorGraph);
            const AnchorId anchorId = v;    // In the AnchorGraph, vertex_descriptors are AnchorIds.
            edge.push_back(anchorId);
        }
    }

    cout << "The initial assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("A");

    transitiveReduction(
        assemblyGraphOptions.transitiveReductionThreshold,
        assemblyGraphOptions.transitiveReductionA,
        assemblyGraphOptions.transitiveReductionB);
    compress();

    cout << "After transitive reduction, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("B");

    createTangleTemplates();

    PermutationDetangler detangler(
        assemblyGraphOptions.minCommonCoverage,
        assemblyGraphOptions.detangleEpsilon,
        assemblyGraphOptions.detangleMaxLogP,
        assemblyGraphOptions.detangleMinLogPDelta);


    for(uint64_t iteration=0; iteration<3; iteration++) {
        cout << "***** Detangle iteration " << iteration << " begins." << endl;

        detangleVertices(detangler);
        compress();

        cout << "After vertex detangling, the assembly graph has " << num_vertices(assemblyGraph) <<
            " vertices and " << num_edges(assemblyGraph) << " edges." << endl;

        detangleEdges(detangler);
        compress();

        cout << "After edge detangling, the assembly graph has " << num_vertices(assemblyGraph) <<
            " vertices and " << num_edges(assemblyGraph) << " edges." << endl;

        cleanupTrivialBubbles();
        compress();
        cout << "After removal of trivial bubbles, the assembly graph has " << num_vertices(assemblyGraph) <<
            " vertices and " << num_edges(assemblyGraph) << " edges." << endl;

        for(uint64_t i=0; i<tangleTemplates.size(); i++) {
            const TangleTemplate& tangleTemplate = tangleTemplates[i];
            performanceLog << timestamp << "Begin detangling template " << i << endl;
            detangle(tangleTemplate, detangler);
            performanceLog << timestamp << "End detangling template " << i << endl;
            compress();
            cout << "After detangling induced subgraphs for Tangle template " << i <<
                ", the assembly graph has " << num_vertices(assemblyGraph) <<
                " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
        }
    }
    write("Y");

    assembleAll(threadCount);
    write("Z");
    writeFasta("AssemblyGraph-Z.fasta");


    performanceLog << timestamp << "AssemblyGraph creation ends." << endl;
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    using shasta::Base;
    vector<Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(gfa));
            gfa << "\tLN:i:" << sequence.size() << "\n";
        } else {
            gfa << "*\t";
            gfa << "LN:i:" << edge.length(anchors) << "\n";
        }
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



void AssemblyGraph::writeFasta(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream fasta(fileName);

    using shasta::Base;
    vector<Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(edge.wasAssembled);
        edge.getSequence(sequence);

        fasta << ">" << edge.id << " " << sequence.size() << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
        fasta << "\n";
    }
}



void AssemblyGraph::write(const string& name) const
{
    performanceLog << timestamp << "Begin AssemblyGraph::save " << name << endl;
    save(name);
    writeGfa("AssemblyGraph-" + name + ".gfa");
    writeGraphviz("AssemblyGraph-" + name + ".dot");
    writeSegments("AssemblyGraphSegments-" + name + ".csv");
    writeSegmentDetails("AssemblyGraphSegmentsDetails-" + name + ".csv");
    performanceLog << timestamp << "End AssemblyGraph::save " << name << endl;
}



void AssemblyGraph::writeSegments(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Id,AnchorIdA,AnchorIdB,AnchorCount,Length\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        csv <<
            edge.id << "," <<
            anchorIdToString(edge.front()) << "," <<
            anchorIdToString(edge.back()) << "," <<
            edge.size() << "," <<
            edge.length(anchors) <<  "\n";
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
    csv << "Id,Position,AnchorIdA,AnchorIdB,Common,Base offset,\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        for(uint64_t i=1; i<edge.size(); i++) {
            const AnchorId anchorIdA = edge[i - 1];
            const AnchorId anchorIdB = edge[i];

            uint64_t baseOffset;
            const uint64_t commonCount = anchors.countCommon(anchorIdA, anchorIdB, baseOffset);

            csv <<
                edge.id << "," <<
                i - 1 << "," <<
                anchorIdToString(anchorIdA) << "," <<
                anchorIdToString(anchorIdB) << "," <<
                commonCount << "," <<
                baseOffset << "\n";
        }
    }

}



// An v0->v1 edge with length length01 is removed if:
// - It has no internal anchors.
// - length01 is no more than threshold1.
// - A BFS that starts at v0 and only uses edges with minimum common count
//   greater than the common count of the v0->v1 edge
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
    performanceLog << timestamp << "AssemblyGraph::transitiveReduction begins." << endl;

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

        const bool debug = false;
        if(debug) cout << "Transitive reduction for " << edge.id << endl;

        // Increment the iterator here, before possibly removing this edge.
        ++it;

        // If the edge has internal anchors, exclude it from transitive reduction.
        if(edge.size() > 2) {
            if(debug) cout << "Skipped due to internal anchors." << endl;
            continue;
        }

        if(edge.length(anchors) > threshold) {
            if(debug) cout << " Skipped due to length." << endl;
            continue;
        }

        // This edge is a candidate for removal.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        const AnchorId anchorId1 = assemblyGraph[v1].anchorId;
        const uint64_t edgeCommonCount = minCommonCountOnEdgeAdjacent(e);

        if(debug) cout << "Edge " << edge.id << " is " <<
            anchorIdToString(anchorId0) << "->" << anchorIdToString(anchorId1) <<
            ", commonCount " << edgeCommonCount << endl;

        // Initialize the BFS.
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

            if(debug) cout << "Dequeued " << anchorIdToString(assemblyGraph[vA].anchorId) << " at distance " << distanceA <<
                ", out=degree " << out_degree(vA, assemblyGraph) << endl;

            // Loop over out-edges of vA.
            BGL_FORALL_OUTEDGES(vA, eAB, assemblyGraph, AssemblyGraph) {
                if(debug) cout << "Processing out-edge " << assemblyGraph[eAB].id << endl;

                // Don't use edge e in the BFS.
                if(eAB == e) {
                    if(debug) cout << "Skipped because it is the edge we are currently processing." << endl;
                    continue;
                }

                // If the min common count is too small, don't consider this edge.
                if(minCommonCountOnEdgeAdjacent(eAB)< edgeCommonCount) {
                    if(debug) cout << "Skipped because of minCommonCount " << minCommonCountOnEdgeAdjacent(eAB) << endl;
                    continue;
                }

                const vertex_descriptor vB = target(eAB, assemblyGraph);

                // If we already encountered vB, do nothing.
                if(assemblyGraph[vB].bfsDistance != invalid<uint64_t>) {
                    if(debug) cout << "Skipped because already encountered." << endl;
                    continue;
                }

                // If we got too far, don't use this in the BFS.
                const uint64_t distanceB = distanceA + assemblyGraph[eAB].length(anchors);
                if(distanceB > threshold2) {
                    if(debug) cout << "Skipped because too far." << endl;
                    continue;
                }

                // If vB is v1, edge e should be removed.
                if(vB == v1) {
                    removeEdge = true;
                    if(debug) cout << "Edge will be removed." << endl;
                    break;
                }

                // Update the BFS.
                if(debug) cout << "Updating the BFS." << endl;
                q.push(vB);
                verticesEncountered.push_back(vB);
                assemblyGraph[vB].bfsDistance = distanceB;
                if(debug) cout << "Done processing out-edge " << assemblyGraph[eAB].id << endl;

            }

            if(removeEdge) {
                break;
            }
        }


        if(removeEdge) {
            boost::remove_edge(e, assemblyGraph);
            if(debug) cout << "Removed." << endl;
        }  else {
            if(debug) cout << "Not removed." << endl;
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

    performanceLog << timestamp << "AssemblyGraph::transitiveReduction ends." << endl;
}



void AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;

    performanceLog << timestamp << "AssemblyGraph::compress begins." << endl;

    // Find linear chains. All the edges in a linear chain of
    // more than one edge are combined into a single edge.
    vector< vector<edge_descriptor> > linearChains;
    findLinearChains(assemblyGraph, 0, linearChains);

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
        newEdge.push_back(assemblyGraph[v0].anchorId);
        for(const edge_descriptor eOld: linearChain) {
            AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
            copy(oldEdge.begin() + 1, oldEdge.end(), back_inserter(newEdge));
            boost::remove_edge(eOld, assemblyGraph);
        }

        // Now we can remove the middle vertices in the linear chain.
        for(uint64_t i=1; i<linearChain.size(); i++) {
            const vertex_descriptor v = source(linearChain[i], assemblyGraph);
            remove_vertex(v, assemblyGraph);
        }

        SHASTA_ASSERT(newEdge.front() == assemblyGraph[v0].anchorId);
        SHASTA_ASSERT(newEdge.back() == assemblyGraph[v1].anchorId);
    }

    performanceLog << timestamp << "AssemblyGraph::compress ends." << endl;
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
    const AssemblerOptions& assemblerOptions,
    const Anchors& anchors,
    const string& assemblyStage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    assemblerOptions(assemblerOptions),
    anchors(anchors)
{
    load(assemblyStage);
}



void AssemblyGraph::detangleVertices(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    performanceLog << timestamp << "AssemblyGraph::detangleVertices begins." << endl;

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
    performanceLog << timestamp << "AssemblyGraph::detangleVertices ends." << endl;

    cout << "Detangling successful for " << successCount << " tangle vertices." << endl;
}



// When detangling edges, we must be careful because successful
// detangling operations can remove edges.
void AssemblyGraph::detangleEdges(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    performanceLog << timestamp << "AssemblyGraph::detangleEdges begins." << endl;

    // Gather all the edges with more than one entrance and exit.
    std::set<edge_descriptor> detanglingCandidates;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const uint64_t entranceCount = (in_degree (v0, assemblyGraph) + in_degree (v1, assemblyGraph)) - 1;
        const uint64_t exitCount     = (out_degree(v0, assemblyGraph) + out_degree(v1, assemblyGraph)) - 1;
        if((entranceCount > 1) and (exitCount > 1)) {
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
        vector<uint64_t> debugList = {43951};
        tangle.debug = find(debugList.begin(), debugList.end(), assemblyGraph[e].id) != debugList.end();
        detangler.debug = tangle.debug;
        if(tangle.debug) {
            cout << "Detangling edge " << assemblyGraph[e].id << endl;
        }
        const bool success = detangler(tangle);
        if(success) {
            ++successCount;
            for(const edge_descriptor e: tangle.removedEdges) {
                detanglingCandidates.erase(e);
            }
        }
    }

    performanceLog << timestamp << "AssemblyGraph::detangleEdges ends." << endl;

    cout << "Attempted detangling for " << attemptCount << " tangle edges." << endl;
    cout << "Detangling was successful for " << successCount << " tangle edges." << endl;

}



// The length of an AssemblyGraphEdge is the estimated length of its sequence,
// equal to the sum of the base offsets for adjacent pairs of AnchorIds in this edge.
uint64_t AssemblyGraphEdge::length(const Anchors& anchors) const
{
    SHASTA_ASSERT(size() > 1);
    const AssemblyGraphEdge& edge = *this;

    uint64_t sumBaseOffset = 0;
    for(uint64_t i=1; i<size(); i++) {
        const AnchorId anchorIdA = edge[i - 1];
        const AnchorId anchorIdB = edge[i];
        uint64_t baseOffset = 0;
        anchors.countCommon(anchorIdA, anchorIdB, baseOffset);
        sumBaseOffset += baseOffset;
    }
    return sumBaseOffset;
}



// For each set of parallel edges with identical AnchorId sequences,
// keep only one.
void AssemblyGraph::cleanupTrivialBubbles()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        while(true) {

            // Look for a pair of out-edges of v0 with identical AnchorId sequences.
            // Every time we find one, we remove one of the two and we start the iteration
            // from scratch. This way we don't invalidate edge iterators.
            out_edge_iterator itBegin, itEnd;
            tie(itBegin, itEnd) = out_edges(v0, assemblyGraph);
            bool found = false;
            for(auto it0=itBegin; it0!=itEnd; it0++) {
                const edge_descriptor e0 = *it0;
                const vector<AnchorId>& anchorSequence0 = assemblyGraph[e0];
                auto it1 = it0;
                ++it1;
                for(; it1!=itEnd; ++it1) {
                    const edge_descriptor e1 = *it1;
                    const vector<AnchorId>& anchorSequence1 = assemblyGraph[e1];
                    if(anchorSequence1 == anchorSequence0) {
                        boost::remove_edge(e1, assemblyGraph);
                        found = true;
                        break;
                    }
                }
                if(found) {
                    break;
                }
            }

            if(not found) {
                break;
            }
        }
    }
}



void AssemblyGraph::createTangleTemplates()
{

    tangleTemplates.emplace_back(4);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
    }

    tangleTemplates.emplace_back(3);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        tangleTemplates.push_back(reverse(g));
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
    }

    tangleTemplates.emplace_back(8);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
        add_edge(5, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
    }

    // Skip four bubbles.
    tangleTemplates.emplace_back(10);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
        add_edge(5, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
        add_edge(7, 8, g);
        add_edge(7, 8, g);
        add_edge(8, 9, g);
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        tangleTemplates.push_back(reverse(g));
    }

    tangleTemplates.emplace_back(8);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        add_edge(4, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
    }


    tangleTemplates.emplace_back(5);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(1, 4, g);
        add_edge(2, 3, g);
        add_edge(2, 4, g);
        tangleTemplates.push_back(reverse(g));
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 3, g);
        add_edge(2, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
    }

    for(uint64_t i=0; i<tangleTemplates.size(); i++) {
        const TangleTemplate& tangleTemplate = tangleTemplates[i];
        const string dotFileName = "TangleTemplate-" + to_string(i) + ".dot";
        ofstream dot(dotFileName);
        writeGraphviz(dot, tangleTemplate);
        dot.close();
        std::system(("dot -O -T svg -Nshape=rectangle " + dotFileName).c_str());
    }
}



// This reverses all of the edges of the TangleTemplate.
// We can't use boost::reverse_graph because that creates
// a graph of a different type.
AssemblyGraph::TangleTemplate AssemblyGraph::reverse(const TangleTemplate& x)
{
    TangleTemplate y(num_vertices(x));

    BGL_FORALL_EDGES(e, x, TangleTemplate) {
        add_edge(target(e, x), source(e, x), y);
    }
    return y;
}



void AssemblyGraph::writeGraphviz(ostream& s, const TangleTemplate& g)
{
    s << "digraph TangleTemplate {\n";

    BGL_FORALL_VERTICES(v, g, TangleTemplate) {
        s << v << ";\n";
    }

    BGL_FORALL_EDGES(e, g, TangleTemplate) {
        const auto v0 = source(e, g);
        const auto v1 = target(e, g);
        s << v0 << "->" << v1 << ";\n";
    }

    s << "}\n";
}



void AssemblyGraph::detangle(
    const TangleTemplate& tangleTemplate,
    Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    vector< vector<vertex_descriptor> > isomorphisms;
    inducedSubgraphIsomorphisms(assemblyGraph, tangleTemplate, isomorphisms);

    cout << "Found " << isomorphisms.size() << " instances for the requested Tangle template." << endl;

    // ofstream html("Tangles.html");
    // writeHtmlBegin(html, "Tangles");

    uint64_t successCount = 0;
    std::set<vertex_descriptor> removedVertices;
    for(const vector<vertex_descriptor>& isomorphism: isomorphisms) {
        /*
        cout << "Working on the following tangle:";
        for(const vertex_descriptor v: isomorphism) {
            cout << " " << anchorIdToString(assemblyGraph[v].anchorId);
        }
        cout << endl;
        */

        // If this contains any vertices that were already removed, skip it.
        bool skip = false;
        for(const vertex_descriptor v: isomorphism) {
            if(removedVertices.contains(v)) {
                skip = true;
                break;
            }
        }
        if(skip) {
            continue;
        }


        Tangle tangle(assemblyGraph, isomorphism);
        // tangle.tangleMatrix.writeHtml(assemblyGraph, html);

        const bool success = detangler(tangle);
        if(success) {
            for(const vertex_descriptor v: isomorphism) {
                removedVertices.insert(v);
            }
            ++successCount;
        }

    }

    // writeHtmlEnd(html);

}



// For a given edge, this returns the minimum common count
// for pairs of adjacent anchors in the edge.
uint64_t AssemblyGraph::minCommonCountOnEdge(edge_descriptor e) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphEdge& edge = assemblyGraph[e];
    SHASTA_ASSERT(edge.size() > 1);

    uint64_t minCommonCount = std::numeric_limits<uint64_t>::max();
    for(uint64_t i=1; i<edge.size(); i++) {
        const AnchorId anchorId0 = edge[i-1];
        const AnchorId anchorId1 = edge[i];
        const uint64_t commonCount = anchors.countCommon(anchorId0, anchorId1);
        minCommonCount = min(minCommonCount, commonCount);
    }

    return minCommonCount;
}



// Same, but only counting journey offsets equal to 1.
uint64_t AssemblyGraph::minCommonCountOnEdgeAdjacent(edge_descriptor e) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphEdge& edge = assemblyGraph[e];
    SHASTA_ASSERT(edge.size() > 1);

    uint64_t minCommonCount = std::numeric_limits<uint64_t>::max();
    for(uint64_t i=1; i<edge.size(); i++) {
        const AnchorId anchorId0 = edge[i-1];
        const AnchorId anchorId1 = edge[i];
        const AnchorPair anchorPair(anchors, anchorId0, anchorId1, true);
        const uint64_t commonCount = anchorPair.size();
        minCommonCount = min(minCommonCount, commonCount);
    }

    return minCommonCount;

}



// Assemble sequence for all edges.
void AssemblyGraph::assembleAll(uint64_t threadCount)
{
    cout << timestamp << "Sequence assembly begins." << endl;
    const AssemblyGraph& assemblyGraph = *this;
    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble(threadCount);
    cout << timestamp << "Sequence assembly ends." << endl;
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the edgeStepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph::assemble(uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    edgeStepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t stepCount = edge.size() - 1;
        edge.sequences.resize(stepCount);
        for(uint64_t i=0; i<stepCount; i++) {
            edgeStepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(edgeStepsToBeAssembled.size(), batchCount);
    runThreads(&AssemblyGraph::assembleThreadFunction, threadCount);

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }
}


void AssemblyGraph::assembleThreadFunction(uint64_t threadId)
{
    AssemblyGraph& assemblyGraph = *this;

    ofstream out("Assemble-Thread-" + to_string(threadId) + ".txt");

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {
            if((j % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << j << "/" << edgeStepsToBeAssembled.size() << endl;
            }

            const auto& p = edgeStepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            AssemblyGraphEdge& edge = assemblyGraph[e];
            SHASTA_ASSERT(i < edge.sequences.size());
            out << "Begin " << edge.id << "," << i << endl;
            const auto t0 = steady_clock::now();
            assembleStep(e, i);
            const auto t1 = steady_clock::now();
            out << "End " << edge.id << "," << i << "," << seconds(t1-t0) << endl;
        }
    }
}



// Assemble sequence for the specified edge.
void AssemblyGraph::assemble(edge_descriptor e, uint64_t threadCount)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble(threadCount);
}



// Assemble sequence for step i of the specified edge.
// The sequences vector for the edge must have already been sized to the correct length.
// This is the lowest level sequence assembly functions and is not multithreaded.
// It runs a LocalAssembly between the appropriate pir of adjacent anchors in the
// AssemblyGraphEdge.
void AssemblyGraph::assembleStep(edge_descriptor, uint64_t)
{
    // Use AssemblyGraph2 instead.
    SHASTA_ASSERT(0);

#if 0
    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphEdge& edge = assemblyGraph[e];

    // Get the AnchorIds for this step.
    SHASTA_ASSERT(i + 1 < edge.size());
    const AnchorId anchorId0 = edge[i];
    const AnchorId anchorId1 = edge[i+1];

    // Access the vector where we want to store assembled sequence of this step.
    SHASTA_ASSERT(i < edge.sequences.size());
    vector<shasta::Base>& sequence = edge.sequences[i];

    // Run the LocalAssembly.
    ofstream html;  // Not open, so no html output takes place.
#if 0
    LocalAssemblyDisplayOptions localAssemblyDisplayOptions(html);
    const LocalAssembly localAssembly(
        anchors.k, anchors.reads, anchors.markers, anchors,
        anchorId0, anchorId1,
        0, localAssemblyDisplayOptions, assemblerOptions.localAssemblyOptions,
        false, false);
#else

    SHASTA_ASSERT(0);
    /*
    LocalAssembly2 localAssembly(anchors, anchorId0, anchorId1, false,
        assemblerOptions.localAssemblyOptions.maxAbpoaLength,
        assemblerOptions.aDrift,
        assemblerOptions.bDrift,
        html, false);
    */
#endif
    // localAssembly.getSequence(sequence);
#endif
}



uint64_t AssemblyGraphEdge::sequenceLength() const
{
    SHASTA_ASSERT(wasAssembled);
    uint64_t length = 0;
    for(const vector<Base>& sequence: sequences) {
        length += sequence.size();
    }
    return length;
}



void AssemblyGraphEdge::getSequence(vector<Base>& edgeSequence) const
{
    SHASTA_ASSERT(wasAssembled);
    edgeSequence.clear();

    for(const vector<Base>& sequence: sequences) {
        copy(sequence.begin(), sequence.end(), back_inserter(edgeSequence));
    }
}
