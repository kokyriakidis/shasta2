// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "findReachableVertices.hpp"
#include "ReadId.hpp"
#include "timestamp.hpp"
#include "tmpDirectory.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AnchorGraph>;



// Construct the AnchorGraph from the Journeys.
// Only include edges with at least the specified minCoverage.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AnchorGraph>(*this)
{
    AnchorGraph& anchorGraph = *this;

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    // Loop over possible source vertices to create edges.
    nextEdgeId = 0;
    vector<AnchorPair> anchorPairs;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, 0, anchorPairs);
        for(const AnchorPair& anchorPair: anchorPairs) {
            if(anchorPair.size() >= minEdgeCoverage) {
                const uint64_t offset = anchorPair.getAverageOffset(anchors);
                edge_descriptor e;
                tie(e, ignore) = add_edge(anchorIdA, anchorPair.anchorIdB,
                    AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                anchorGraph[e].useForAssembly = true;
            }
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}



// Constructor from binary data.
AnchorGraph::AnchorGraph(const MappedMemoryOwner& mappedMemoryOwner, const string& name) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<AnchorGraph>(*this)
{
    load(name);
}



void AnchorGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AnchorGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AnchorGraph::save(const string& name) const
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



void AnchorGraph::load(const string& name)
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



void AnchorGraph::transitiveReduction(
    uint64_t transitiveReductionMaxEdgeCoverage,
    uint64_t transitiveReductionMaxDistance)
{
    AnchorGraph& anchorGraph = *this;
    cout << "AnchorGraph transitive reduction begins." << endl;

    // Initially make sure all edges are flag as "useForAssembly".
    // The transitive reduction process sets useForAssembly to false
    // for edges removed by transitive reduction.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        anchorGraph[e].useForAssembly = true;
    }

    // Loop over edge coverage.
    // At each iteration we only consider edges with this coverage.
    vector<edge_descriptor> edgesToProcess;
    vector<edge_descriptor> edgesToRemove;
    for(uint64_t edgeCoverage=1; edgeCoverage<=transitiveReductionMaxEdgeCoverage; edgeCoverage++) {

        // Gather edges with this coverage.
        edgesToProcess.clear();
        BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
            if(anchorGraph[e].coverage() == edgeCoverage) {
                edgesToProcess.push_back(e);
            }
        }

        // If there are none, there is nothing to do.
        if(edgesToProcess.empty()) {
            continue;
        }

        // Loop over all edges with this coverage.
        // This can be multithreaded.
        edgesToRemove.clear();
        for(const edge_descriptor e: edgesToProcess) {
            if(transitiveReductionCanRemove(e, transitiveReductionMaxDistance)) {
                edgesToRemove.push_back(e);
            }
        }

        // Turn off the useForAssembly flag for edges removed at this iteration over coverage.
        for(const edge_descriptor e: edgesToRemove) {
            anchorGraph[e].useForAssembly = false;
        }
        cout << "Edge coverage " << edgeCoverage <<
            ": processed " << edgesToProcess.size() <<
            " edges and flagged " << edgesToRemove.size() << endl;
    }
    cout << "AnchorGraph transitive reduction ends." << endl;

    uint64_t useForAssemblyCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(anchorGraph[e].useForAssembly) {
            ++useForAssemblyCount;
        }
    }
    cout << useForAssemblyCount << " flagged for use in assembly out of " <<
        num_edges(anchorGraph) << " total." << endl;

}



bool AnchorGraph::transitiveReductionCanRemove(
    edge_descriptor e,
    uint64_t transitiveReductionMaxDistance) const
{
    const AnchorGraph& anchorGraph = *this;
    const uint64_t edgeCoverage = anchorGraph[e].coverage();

    const vertex_descriptor v0 = source(e, anchorGraph);
    const vertex_descriptor v1 = target(e, anchorGraph);

    const bool debug = ((anchorIdToString(v0) == "45549+") and (anchorIdToString(v1) == "78505-"));

    // Do a forward BFS starting at v0, using edges
    // still marked as "use for assembly"
    // with coverage greater than edgeCoverage
    // and with maximum distance (number of edges)
    // equal to transitiveReductionMaxDistance.
    // If we encounter v1, return true.
    std::queue<vertex_descriptor> q;
    q.push(v0);

    // A map to store vertices already encountered and their distance from v0.
    std::map<vertex_descriptor, uint64_t> m;
    m.insert(make_pair(v0, 0));



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor vA = q.front();
        q.pop();
        const auto itA = m.find(vA);
        SHASTA2_ASSERT(itA != m.end());
        const uint64_t distanceA = itA->second;
        const uint64_t distanceB = distanceA + 1;

        // Loop over its out-edges still marked as useForAssembly
        // and with sufficient coverage.
        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            const AnchorGraphEdge& edgeAB = anchorGraph[eAB];
            if(not edgeAB.useForAssembly) {
                continue;
            }

            // Only use edges with higher coverage for the BFS,
            if(edgeAB.coverage() <= edgeCoverage) {
                continue;
            }

            // If we reached v1, return true;
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(vB == v1) {
                if(debug) {
                    cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
                        " flagged by transitive reduction." << endl;
                }
                return true;
            }

            // If we already encountered vB, don't do anything.
            if(m.contains(vB)) {
                continue;
            }

            if(distanceB < transitiveReductionMaxDistance) {
                q.push(vB);
                m.insert(make_pair(vB, distanceB));
            }
        }
    }

    // If getting here we did not encounter v1 in the BFS loop.
    if(debug) {
        cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
            " not flagged by transitive reduction." << endl;
    }
    return false;
}



// Constructor to create an anchor similarity graph,
// in which an edge between two anchors is created if the
// oriented read compositions of the two anchors are
// sufficiently similar.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    const UseSimilarity&) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AnchorGraph>(*this)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double a = 1.;
    const double b = 1.;
    const double minLogP = 0.;
    const uint64_t m = 100;

    AnchorGraph& anchorGraph = *this;

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }


    class EdgeCandidate {
    public:
        AnchorPair anchorPair;
        AnchorPairInfo anchorPairInfo;
        uint64_t offset;
        double logP = invalid<double>;
        EdgeCandidate(const Anchors& anchors, AnchorId anchorIdA, AnchorId anchorIdB) :
            anchorPair(anchors, anchorIdA, anchorIdB, false)
        {
            offset = anchorPair.getAverageOffset(anchors);
            anchors.analyzeAnchorPair(anchorIdA, anchorIdB, anchorPairInfo);
        }
        bool operator<(const EdgeCandidate& that) const
        {
            return logP > that.logP;
        }
    };

    // Loop over possible source vertices to create edges.
    vector<AnchorId> anchorIds;
    vector<uint64_t> count;
    AnchorPairInfo anchorPairInfo;
    vector<EdgeCandidate> edgeCandidates;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        const Anchor anchorA = anchors[anchorIdA];

        const bool debug = anchorIdToString(anchorIdA) == "436-";

        // Loop over OrientedReadIds in this anchor.
        // For each OrientedReadId, gather the AnchorIds
        // in the journey portion beginning here.
        anchorIds.clear();
        for(const AnchorMarkerInfo& anchorMarkerInfo: anchorA) {
            const OrientedReadId orientedReadId = anchorMarkerInfo.orientedReadId;
            const Journey journey = journeys[orientedReadId];
            const uint32_t positionInJourneyA = anchorMarkerInfo.positionInJourney;

            // Gather the AnchorIds seen by this OrientedRead after anchorIdA.
            for(uint32_t positionInJourney=positionInJourneyA+1; positionInJourney<journey.size(); positionInJourney++) {
                anchorIds.push_back(journey[positionInJourney]);
            }
        }
        deduplicateAndCountWithThreshold(anchorIds, count, minEdgeCoverage);

        // For each of these AnchorIds that meets out requirements, create an EdgeCandidate.
        edgeCandidates.clear();
        for(const AnchorId anchorIdB: anchorIds) {
            EdgeCandidate& edgeCandidate = edgeCandidates.emplace_back(anchors, anchorIdA, anchorIdB);
            const AnchorPairInfo& info = edgeCandidate.anchorPairInfo;
            if(debug) {
                cout << anchorIdToString(anchorIdA) << " " << anchorIdToString(anchorIdB) << " " << info.common << endl;
            }
            if(info.common < minEdgeCoverage) {
                if(debug) {
                    cout << "Skipped due to low commonCount." << endl;
                }
                edgeCandidates.pop_back();
                continue;
            }
            const uint64_t missing = info.onlyA + info.onlyB - info.onlyAShort - info.onlyBShort;
            edgeCandidate.logP = a * double(info.common) - b * double(missing);
            if(debug) {
                cout << "Common " << info.common << ", missing " << missing << ", logP " << edgeCandidate.logP << endl;
            }
            if(edgeCandidate.logP < minLogP)  {
                edgeCandidates.pop_back();
                continue;
            }
        }
        if(debug) {
            cout << edgeCandidates.size() << " edge candidates for " << anchorIdToString(anchorIdA) << endl;
        }

        // Sort the edge candidates by decreasing logP.
        // Keep the m ones with the greatest logP.
        sort(edgeCandidates.begin(), edgeCandidates.end());
        for(uint64_t i=0; i<min(m, edgeCandidates.size()); i++) {
            const EdgeCandidate& edgeCandidate = edgeCandidates[i];
            edge_descriptor e;
            tie(e, ignore) = add_edge(edgeCandidate.anchorPair.anchorIdA, edgeCandidate.anchorPair.anchorIdB,
                AnchorGraphEdge(edgeCandidate.anchorPair, edgeCandidate.offset, nextEdgeId++), anchorGraph);
            anchorGraph[e].useForAssembly = true;
        }
    }

    cout << "The anchor similarity graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}



// This uses read following in the complete AnchorGraph
// to create the AnchorGraph to be used for assembly.
// This is meant to be used with strict anchor generation,
// where most anchors correspond to a single copy.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const AnchorGraph& completeAnchorGraph) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AnchorGraph>(*this)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = 3;

    // The AnchorIds of the edge we find.
    // Indexed by edgeAnchorIds[direction][anchorId0] .
    array<vector< vector<AnchorId> >, 2> edgeAnchorIds;
    edgeAnchorIds[0].resize(anchors.size());
    edgeAnchorIds[1].resize(anchors.size());

    // Loop over all anchors.
    vector<Subgraph::vertex_descriptor> path;
    BGL_FORALL_VERTICES(anchorIdStart, completeAnchorGraph, AnchorGraph) {

        // Loop over both directions (0=forward, 1=backward).
        for(uint64_t direction=0; direction<2; direction++) {

            // Create a Subgraph starting at anchorIdStart and moving in this direction.
            Subgraph subgraph(anchors, completeAnchorGraph, anchorIdStart, direction, minCommonCount);
            subgraph.removeCycles();
            subgraph.transitiveReduction();
            subgraph.pruneMultipleExits();

            // Create the dominator tree.
            const Subgraph dominatorTree(
                subgraph,
                AnchorGraph::Subgraph::DominatorTree(),
                anchors);

            // Walk up from the exit.
            subgraph.walkUp(dominatorTree, path);

            for(uint64_t i1=1; i1<path.size(); i1++) {
                const uint64_t i0 = i1 - 1;
                const Subgraph::vertex_descriptor v0 = path[i0];
                const Subgraph::vertex_descriptor v1 = path[i1];
                const AnchorId anchorId0 = dominatorTree[v0].anchorId;
                const AnchorId anchorId1 = dominatorTree[v1].anchorId;
                edgeAnchorIds[direction][anchorId0].push_back(anchorId1);
            }
        }
    }

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(*this);
    }

    // Create the edges. An edge is added if it was found in both directions.
    nextEdgeId = 0;
    vector<AnchorId> goodAnchorIds;
    for(AnchorId anchorIdA=0; anchorIdA<anchors.size(); anchorIdA++) {
        for(uint64_t direction=0; direction<2; direction++) {
            vector<AnchorId>& anchorIdsB = edgeAnchorIds[direction][anchorIdA];
            deduplicate(anchorIdsB);
            /*
            cout << anchorIdToString(anchorIdA) << " " << direction << ":";
            for(const AnchorId anchorIdB: anchorIdsB) {
                cout << " " << anchorIdToString(anchorIdB);
            }
            cout << endl;
            */
        }

        // Get the ones that were found in both directions.
        goodAnchorIds.clear();
        std::set_intersection(
            edgeAnchorIds[0][anchorIdA].begin(), edgeAnchorIds[0][anchorIdA].end(),
            edgeAnchorIds[1][anchorIdA].begin(), edgeAnchorIds[1][anchorIdA].end(),
            back_inserter(goodAnchorIds));

        /*
        cout << anchorIdToString(anchorIdA) << ":";
        for(const AnchorId anchorIdB: goodAnchorIds) {
            cout << " " << anchorIdToString(anchorIdB);
        }
        cout << endl;
        */

        // Generate these edges.
        for(const AnchorId anchorIdB: goodAnchorIds) {
            const AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);
            const uint64_t offset = anchorPair.getAverageOffset(anchors);
            edge_descriptor e;
            tie(e, ignore) = add_edge(anchorIdA, anchorIdB,
                AnchorGraphEdge(anchorPair, offset, nextEdgeId++), *this);
            (*this)[e].useForAssembly = true;
        }
    }

    cout << "The complete anchor graph has " <<
        num_vertices(completeAnchorGraph) << " vertices and " <<
        num_edges(completeAnchorGraph) << " edges." << endl;
    cout << "The anchor graph has " <<
        num_vertices(*this) << " vertices and " <<
        num_edges(*this) << " edges." << endl;
}



// Do a BFS starting at anchorIdStart and moving in the specified direcrtion,
// disregarding vertices that have no common OrientedReadIds with anchorIdStart.
AnchorGraph::Subgraph::Subgraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    AnchorId anchorIdStart,
    uint64_t direction,
    uint64_t minCommonCount) :
    direction(direction)
{
    const bool debug = false;
    if(debug) {
        cout << "Creating AnchorGraph::Subgraph for " << anchorIdToString(anchorIdStart) <<
            " direction " << direction << endl;
    }
    Subgraph& subgraph = *this;

    // Initialize the BFS.
    std::queue<AnchorId> q;
    q.push(anchorIdStart);
    std::set<AnchorId> visited;
    visited.insert(anchorIdStart);
    vStart = add_vertex(SubgraphVertex(anchorIdStart, anchors[anchorIdStart].size()), subgraph);
    vertexMap.insert(make_pair(anchorIdStart, vStart));

    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue an AnchorId.
        const AnchorId anchorId0 = q.front();
        if(debug) {
            cout << "Dequeued " << anchorIdToString(anchorId0) << endl;
        }
        q.pop();
        const vertex_descriptor v0 = vertexMap.at(anchorId0);

        vector<AnchorId> next;
        if(direction == 0) {
            BGL_FORALL_OUTEDGES(anchorId0, e, anchorGraph, AnchorGraph) {
                next.push_back(target(e, anchorGraph));
            }
        } else {
            BGL_FORALL_INEDGES(anchorId0, e, anchorGraph, AnchorGraph) {
                next.push_back(source(e, anchorGraph));
            }
        }

        for(const AnchorId anchorId1: next) {
            if(debug) {
                cout << "Found " << anchorIdToString(anchorId1) << endl;
            }
            const uint64_t commonCount =
                (direction == 0) ?
                anchors.countCommon(anchorIdStart, anchorId1) :
                anchors.countCommon(anchorId1, anchorIdStart);

            if(commonCount < minCommonCount) {
                if(debug) {
                    cout << "Not enough common oriented reads, discarded." << endl;
                }
                continue;
            } else {
                if(debug) {
                    cout << commonCount << " common oriented reads." << endl;
                }
            }

            auto it1 = vertexMap.find(anchorId1);
            if(it1 == vertexMap.end()) {
                tie(it1, ignore) = vertexMap.insert(make_pair(anchorId1,
                    add_vertex(SubgraphVertex(anchorId1, commonCount), subgraph)));
            }
            const vertex_descriptor v1 = it1->second;


            // Get vertex descriptors consistent with the direction.
            vertex_descriptor u0 = v0;
            vertex_descriptor u1 = v1;
            if(direction == 1) {
                std::swap(u0, u1);
            }

            // Get the corresponding AnchorGraph edge.
            AnchorGraph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(subgraph[u0].anchorId, subgraph[u1].anchorId, anchorGraph);
            SHASTA2_ASSERT(edgeExists);
            const AnchorGraphEdge anchorGraphEdge = anchorGraph[e];

            add_edge(u0, u1, SubgraphEdge(anchorGraphEdge.anchorPair.size()), subgraph);
            if(debug) {
                cout << "Added edge " << anchorIdToString(anchorId0) << " " << anchorIdToString(anchorId1) << endl;
            }

            if(not visited.contains(anchorId1)) {
                q.push(anchorId1);
                visited.insert(anchorId1);
                if(debug) {
                    cout << "Enqueued " << anchorIdToString(anchorId1) << endl;
                }
                SHASTA2_ASSERT(subgraph[v1].anchorId == anchorId1);
             }
        }

    }
}


// This does an approximate topological sort, then removes
// edges not flagged as DAG edges.
void AnchorGraph::Subgraph::removeCycles()
{
    Subgraph& subgraph = *this;

    // Add the edges in order of decreasing coverage.
    vector<pair<edge_descriptor, uint64_t> > edgesWithCoverage;
    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        const uint64_t coverage = subgraph[e].coverage;
        edgesWithCoverage.push_back(make_pair(e, coverage));
    }
    sort(edgesWithCoverage.begin(), edgesWithCoverage.end(), OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());

    vector<edge_descriptor> sortedEdges;
    for(const auto& [e, ignore]: edgesWithCoverage) {
        sortedEdges.push_back(e);
    }

    shasta2::approximateTopologicalSort(subgraph, sortedEdges);

    // Remove edges not flagged as DAG edges.
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



void AnchorGraph::Subgraph::transitiveReduction()
{
    shasta2::transitiveReduction(*this);
}




void AnchorGraph::Subgraph::writeGraphviz(const string& fileName, const Anchors& anchors) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, anchors);
}



void AnchorGraph::Subgraph::writeGraphviz(ostream& dot, const Anchors& anchors) const
{
    const Subgraph& subgraph = *this;
    const uint64_t startAnchorIdCoverage = anchors[subgraph[vStart].anchorId].size();

    // For better display, write the vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > sortedVertices;
    BGL_FORALL_VERTICES(v, subgraph, Subgraph) {
        if((v != vStart) and (out_degree(v, subgraph) == 0) and (in_degree(v, subgraph) == 0)) {
            continue;
        }
        const SubgraphVertex& vertex = subgraph[v];
        sortedVertices.push_back(make_pair(v, vertex.rank));
    }
    sort(sortedVertices.begin(), sortedVertices.end(), OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());




    dot << "digraph S {\n";
    for(const auto& [v, ignore]: sortedVertices) {
        const SubgraphVertex& vertex = subgraph[v];
        const AnchorId anchorId = vertex.anchorId;
        dot << "\"" << anchorIdToString(anchorId) << "\"";
        dot << "[";
        dot << "label=\"" << anchorIdToString(subgraph[v].anchorId) <<
            "\\nCoverage " << anchors[anchorId].size() <<
            "\\nCommon " << vertex.commonCount;
        if(isValid(vertex.rank)) {
            dot << "\\nRank " << vertex.rank;
        }
        dot << "\"";

        string color;
        if(v == vStart) {
            color = "LightBlue";
        } else {
            if(vertex.commonCount == 1) {
                color = "Grey";
            } else if(vertex.commonCount == 2) {
                color = "LightGrey";
            } else {
                const double hue = double(vertex.commonCount) / double(startAnchorIdCoverage);
                color = "\"" + to_string(hue / 3.) + " 1 1\"";
            }
        }
        dot << " style=filled fillcolor=" << color;

        dot << "]";
        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        const SubgraphEdge& edge = subgraph[e];
        const vertex_descriptor v0 = source(e, subgraph);
        const vertex_descriptor v1 = target(e, subgraph);
        dot <<
            "\"" << anchorIdToString(subgraph[v0].anchorId) << "\""
            "->"
            "\"" << anchorIdToString(subgraph[v1].anchorId) << "\""
            "["
            "label=\"" << subgraph[e].coverage << "\""
            " penwidth=" << std::fixed << std::setprecision(2) << 0.3 * double(edge.coverage);
        if(not edge.isDagEdge) {
            dot << " color=red";
        }
        dot <<
            "]"
            ";\n";
    }
    dot << "}\n";

}



void AnchorGraph::Subgraph::writeHtml(ostream& html, const Anchors& anchors) const
{
    html << "The subgraph has " << num_vertices(*this) << " vertices and " <<
        num_edges(*this) << " edges.<br>";

    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, anchors);


    // Display it in html in svg format.
    const double timeout = 120.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    try {
        graphvizToHtml(dotFileName, "dot", timeout, options, html);
    } catch(const std::exception& exception) {
        html << "<br>Error during graph layout: " << exception.what();
    }

}



// Recursively prune leafs with commonCount less than minLeafCommonCount.
void AnchorGraph::Subgraph::prune(uint64_t minLeafCommonCount)
{
    const bool debug = false;
    Subgraph& subgraph = *this;

    // Gather the low coverage leafs that currently exist.
    std::set<vertex_descriptor> leaves;
    BGL_FORALL_VERTICES(v, subgraph, Subgraph) {
        if(subgraph[v].commonCount < minLeafCommonCount) {
            const bool isLeaf =
                ((direction == 0) and (out_degree(v, subgraph) == 0))
                or
                ((direction == 1) and (in_degree(v, subgraph) == 0));
            if(isLeaf) {
                leaves.insert(v);
            }
        }
    }

    if(debug) {
        cout << "Initial leaves:";
        for(const vertex_descriptor v: leaves) {
            cout << " " << anchorIdToString(subgraph[v].anchorId);
        }
    }
    cout << endl;



    // Remove leaves, updating the set of currently existing leaves.
    while(not leaves.empty()) {
        const auto it = leaves.begin();
        const vertex_descriptor v0 = *it;
        leaves.erase(it);

        if(debug) {
            cout << "Working on leaf " << anchorIdToString(subgraph[v0].anchorId) << endl;
        }

        // Update the set before removing it.
        if(direction == 0) {
            BGL_FORALL_INEDGES(v0, e, subgraph, Subgraph) {
                const vertex_descriptor v1 = source(e, subgraph);
                if(subgraph[v1].commonCount < minLeafCommonCount) {
                    const uint64_t currentOutDegree = out_degree(v1, subgraph);
                    const uint64_t outDegreeAfterRemoval = currentOutDegree - 1;
                    if(outDegreeAfterRemoval == 0) {
                        leaves.insert(v1);
                        if(debug) {
                            cout << "Adding " << anchorIdToString(subgraph[v1].anchorId) << " to set of leaves." << endl;
                        }
                    }
                }
            }
        } else {
            BGL_FORALL_OUTEDGES(v0, e, subgraph, Subgraph) {
                const vertex_descriptor v1 = target(e, subgraph);
                if(subgraph[v1].commonCount < minLeafCommonCount) {
                    const uint64_t currentInDegree = in_degree(v1, subgraph);
                    const uint64_t inDegreeAfterRemoval = currentInDegree - 1;
                    if(inDegreeAfterRemoval == 0) {
                        leaves.insert(v1);
                        if(debug) {
                            cout << "Adding " << anchorIdToString(subgraph[v1].anchorId) << " to set of leaves." << endl;
                        }
                    }
                }
            }

        }

        // Now we can remove this low coverage leaf.
        // We can't remove the vertex because that would invalidate all vertex descriptors
        // (because we use vecS for the Subgraph).
        if(debug) {
            cout << "Removing " << anchorIdToString(subgraph[v0].anchorId) << endl;
        }
        boost::clear_vertex(v0, subgraph);
    }
}



void AnchorGraph::Subgraph::writeFasta(
    const vector<vertex_descriptor>& path,
    const string& fileName,
    const Anchors& anchors) const
{
    ofstream fasta(fileName);
    writeFasta(path, fasta, anchors);
}



void AnchorGraph::Subgraph::writeFastaHtml(
    const vector<vertex_descriptor>& path,
    ostream& html,
    const Anchors& anchors) const
{
    html << "<pre>";
    writeFasta(path, html, anchors);
    html << "</pre>";
}



void AnchorGraph::Subgraph::writeFasta(
    const vector<vertex_descriptor>& path,
    ostream& fasta,
    const Anchors& anchors) const
{
    using shasta2::Base;
    const Subgraph& subgraph = *this;

    for(uint64_t i=0; i<path.size(); i++) {
        const vertex_descriptor v = path[i];
        const AnchorId anchorId = subgraph[v].anchorId;
        const vector<Base> kmerSequence = anchors.anchorKmerSequence(anchorId);
        fasta << ">" << i << "-" << anchorIdToString(anchorId) << "\n";
        std::ranges::copy(kmerSequence, ostream_iterator<Base>(fasta));
        fasta << "\n";
    }
}



// Get the AnchorIds of the linear chain on vertices containing vStart.
// They are returned in order.
void AnchorGraph::Subgraph::getLinearPortion(vector<AnchorId>& anchorIds) const
{
    const Subgraph& subgraph = *this;
    anchorIds.clear();

    vertex_descriptor v = vStart;

    if(direction == 0) {
        while(true) {
            const uint64_t outDegree = out_degree(v, subgraph);
            if(outDegree > 1) {
                break;
            }
            anchorIds.push_back(subgraph[v].anchorId);
            if(outDegree == 0) {
                break;
            }
            SHASTA2_ASSERT(outDegree == 1);
            auto [it, ignore] = out_edges(v, subgraph);
            const edge_descriptor e = *it;
            v = target(e, subgraph);
        }
    } else if(direction == 1) {
        while(true) {
            const uint64_t inDegree = in_degree(v, subgraph);
            if(inDegree > 1) {
                break;
            }
            anchorIds.push_back(subgraph[v].anchorId);
            if(inDegree == 0) {
                break;
            }
            SHASTA2_ASSERT(inDegree == 1);
            auto [it, ignore] = in_edges(v, subgraph);
            const edge_descriptor e = *it;
            v = source(e, subgraph);
        }
        std::ranges::reverse(anchorIds);
    } else {
        SHASTA2_ASSERT(0);
    }
}



// Constructor the dominator tree (in the given direction).
AnchorGraph::Subgraph::Subgraph(
    const Subgraph& subgraph,
    const DominatorTree&,
    const Anchors& anchors)
{
    Subgraph& dominatorTree = *this;

    // Create vertices of the dominator tree.
    BGL_FORALL_VERTICES(v, subgraph, Subgraph) {
        const SubgraphVertex& vertex = subgraph[v];
        const vertex_descriptor u = add_vertex(vertex, dominatorTree);
        vertexMap[vertex.anchorId] = u;
    }
    vStart = vertexMap.at(subgraph[subgraph.vStart].anchorId);
    direction = subgraph.direction;



    if(direction == 0) {

        // Compute the dominator tree of the input subgraph.
        std::map<vertex_descriptor, vertex_descriptor> dominatorTreeMap;
            shasta2::lengauer_tarjan_dominator_tree(subgraph, subgraph.vStart,
                boost::make_assoc_property_map(dominatorTreeMap));

        // Create edges of the dominator tree.
        for(auto& [v1, v0]: dominatorTreeMap) {
            const AnchorId anchorId0 = subgraph[v0].anchorId;
            const AnchorId anchorId1 = subgraph[v1].anchorId;
            const vertex_descriptor u0 = vertexMap.at(anchorId0);
            const vertex_descriptor u1 = vertexMap.at(anchorId1);
            const AnchorPair anchorPair(anchors, anchorId0, anchorId1, false);
            add_edge(u0, u1, SubgraphEdge(anchorPair.size(), true), dominatorTree);
        }
    } else {

        // Compute the dominator tree of the reverse_graph of the input subgraph.
        using ReverseGraph = boost::reverse_graph<Subgraph>;
        const ReverseGraph reverseGraph(subgraph);
        std::map<vertex_descriptor, vertex_descriptor> dominatorTreeMap;
            shasta2::lengauer_tarjan_dominator_tree(reverseGraph, subgraph.vStart,
                boost::make_assoc_property_map(dominatorTreeMap));

        // Create edges of the dominator tree.
        for(auto& [v0, v1]: dominatorTreeMap) {
            const AnchorId anchorId0 = subgraph[v0].anchorId;
            const AnchorId anchorId1 = subgraph[v1].anchorId;
            const vertex_descriptor u0 = vertexMap.at(anchorId0);
            const vertex_descriptor u1 = vertexMap.at(anchorId1);
            const AnchorPair anchorPair(anchors, anchorId0, anchorId1, false);
            add_edge(u0, u1, SubgraphEdge(anchorPair.size(), true), dominatorTree);
        }
    }
}



// Find exits, ignoring isolated verties.
void AnchorGraph::Subgraph::findExits(vector<vertex_descriptor>& exits) const
{
    const Subgraph& subgraph = *this;

    exits.clear();
    BGL_FORALL_VERTICES(v, subgraph, Subgraph) {
        const uint64_t outDegree = out_degree(v, subgraph);
        const uint64_t inDegree = in_degree(v, subgraph);
        if((v != vStart) and (outDegree == 0) and (inDegree == 0)) {
            continue;
        }
        const bool isExit =
            ((direction == 0) and (outDegree == 0))
            or
            ((direction == 1) and (inDegree == 0));
        if(isExit) {
            exits.push_back(v);
        }
    }

}



// If there are multiple exits, keep only vertices that are
// reachable (backward) from all exits.
void AnchorGraph::Subgraph::pruneMultipleExits()
{
    Subgraph& subgraph = *this;


    while(true) {
        vector<vertex_descriptor> exits;
        findExits(exits);

        if(exits.size() < 2) {
            break;
        }

        // Find how many times each vertex is reachable from one of the exits.
        std::set<vertex_descriptor> reachableVertices;
        std::map<vertex_descriptor, uint64_t> reachCount;
        for(const vertex_descriptor exit: exits) {
            findReachableVertices(subgraph, exit, 1 - direction, reachableVertices);
            for(const vertex_descriptor v: reachableVertices) {
                const auto it = reachCount.find(v);
                if(it == reachCount.end()) {
                    reachCount.insert(make_pair(v, 1));
                } else {
                    ++it->second;
                }
            }
        }

        // Remove the vertices that are not reachable from all the exits.
        vector<vertex_descriptor> verticesToBeRemoved;
        BGL_FORALL_VERTICES(v, subgraph, Subgraph) {
            const auto it = reachCount.find(v);
            if((it == reachCount.end()) or (it->second != exits.size())) {
                verticesToBeRemoved.push_back(v);
            }
        }
        for(const vertex_descriptor v: verticesToBeRemoved) {
            // We can't remove it because the Subgraph uses vecS,
            clear_vertex(v, subgraph);
        }
    }
}



// Walk up the dominator tree.
// This will assert if not called on the dominator tree.
void AnchorGraph::Subgraph::walkUp(
    vertex_descriptor v,
    vector<vertex_descriptor>& path) const
{
    const Subgraph& dominatorTree = *this;
    path.clear();

    // Keep going until we reach vStart.
    while(true) {

        // Add v to the path.
        path.push_back(v);

        // If we reached vStart, stop here.
        if(v == vStart) {
            break;
        }

        // Move back.
        if(direction == 0) {
            auto[begin, end] = in_edges(v, dominatorTree);
            SHASTA2_ASSERT(begin != end);
            const edge_descriptor e = *begin++;
            SHASTA2_ASSERT(begin == end);
            v = source(e, dominatorTree);
        } else {
            auto[begin, end] = out_edges(v, dominatorTree);
            SHASTA2_ASSERT(begin != end);
            const edge_descriptor e = *begin++;
            SHASTA2_ASSERT(begin == end);
            v = target(e, dominatorTree);
        }
    }

    if(direction == 0) {
        std::ranges::reverse(path);
    }
}




// Walk up the dominator tree.
// Note this returns a path in the dominator tree.
void AnchorGraph::Subgraph::walkUp(
    const Subgraph& dominatorTree,
    vector<vertex_descriptor>& path) const
{
    const Subgraph& subgraph = *this;

    // Find the exit in the Subgraph.
    vector<vertex_descriptor> exits;
    subgraph.findExits(exits);
    SHASTA2_ASSERT(exits.size() == 1);
    const vertex_descriptor exit = exits.front();

    // Walk up the dominator tree starting at the corresponding vertex.
    const AnchorId exitAnchorId = subgraph[exit].anchorId;
    dominatorTree.walkUp(dominatorTree.vertexMap.at(exitAnchorId), path);

}
