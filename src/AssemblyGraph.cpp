// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "findLinearChains.hpp"
#include "html.hpp"
#include "inducedSubgraphIsomorphisms.hpp"
#include "performanceLog.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
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
    const AssemblerOptions::AssemblyGraphOptions& options ) :
    MappedMemoryOwner(anchors),
    anchors(anchors)
{
    AssemblyGraph& assemblyGraph = *this;
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
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
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
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
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
        options.transitiveReductionThreshold,
        options.transitiveReductionA,
        options.transitiveReductionB);
    compress();

    cout << "After transitive reduction, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("B");

    fillInducedSubgraphTemplates();

    TrivialDetangler detangler(options.minCommonCoverage);

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

    cleanupTrivialBubbles();
    compress();
    cout << "After removal of trivial bubbles, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("E");

    for(uint64_t i=0; i<inducedSubgraphTemplates.size(); i++) {
        const InducedSubgraph& inducedSubgraph = inducedSubgraphTemplates[i];
        write("X" + to_string(i));
        detangleInducedSubgraphs(inducedSubgraph, detangler);
        compress();
        cout << "After detangling induced subgraphs for template " << i <<
            ", the assembly graph has " << num_vertices(assemblyGraph) <<
            " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    }
    write("F");

    performanceLog << timestamp << "AssemblyGraph creation ends." << endl;
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
        gfa << "LN:i:" << edge.length(anchors) << "\n";
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
        // vector<uint64_t> debugList = {110089,110193,75121,105461,105503};
        // tangle.debug = find(debugList.begin(), debugList.end(), assemblyGraph[e].id) != debugList.end();
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



void AssemblyGraph::fillInducedSubgraphTemplates()
{
    inducedSubgraphTemplates.push_back(InducedSubgraph(4));
    {
        InducedSubgraph& g = inducedSubgraphTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
    }

    inducedSubgraphTemplates.push_back(InducedSubgraph(6));
    {
        InducedSubgraph& g = inducedSubgraphTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
    }

    for(uint64_t i=0; i<inducedSubgraphTemplates.size(); i++) {
        const InducedSubgraph& inducedSubgraph = inducedSubgraphTemplates[i];
        ofstream dot("InducedSubgraphTemplate-" + to_string(i) + ".dot");
        writeGraphviz(dot, inducedSubgraph);
    }
}



void AssemblyGraph::writeGraphviz(ostream& s, const InducedSubgraph& g)
{
    s << "digraph InducedSubgraphTemplate {\n";

    BGL_FORALL_VERTICES(v, g, InducedSubgraph) {
        s << v << ";\n";
    }

    BGL_FORALL_EDGES(e, g, InducedSubgraph) {
        const auto v0 = source(e, g);
        const auto v1 = target(e, g);
        s << v0 << "->" << v1 << ";\n";
    }

    s << "}\n";
}



void AssemblyGraph::detangleInducedSubgraphs(
    const InducedSubgraph& inducedSubgraphTemplate,
    Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    vector< vector<vertex_descriptor> > isomorphisms;
    inducedSubgraphIsomorphisms(assemblyGraph, inducedSubgraphTemplate, isomorphisms);

    cout << "Found " << isomorphisms.size() << " induced subgraphs for the requested template." << endl;

    // ofstream html("InducedSubgraphTangles.html");
    // writeHtmlBegin(html, "InducedSubgraphTangles");

    uint64_t successCount = 0;
    std::set<vertex_descriptor> removedVertices;
    for(const vector<vertex_descriptor>& isomorphism: isomorphisms) {
        /*
        cout << "Working on the following induced subgraph:";
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
