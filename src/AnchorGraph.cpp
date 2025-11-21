// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "ReadLengthDistribution.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AnchorGraph>;



// Simple generation of edges.
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



// Constructor that uses read following.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    uint64_t minContinueReadFollowingCount,
    double aDrift,
    double bDrift,
    uint64_t /* threadCount */) :
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

    // Create the edges using read following.
    createEdges3(anchors, journeys, minEdgeCoverage, minContinueReadFollowingCount, aDrift, bDrift);

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

}



void AnchorGraph::createEdges1(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    double aDrift,
    double bDrift)
{
    AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = anchors.size();

    // Work vectors for the loop below.
    vector<AnchorPair> anchorPairs;
    vector<AnchorPair> newAnchorPairs;
    vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
    vector<uint64_t> offsets;

    // Create the edges using read following.
    for(uint64_t direction=0; direction<2; direction++) {
        for(AnchorId anchorId0=0; anchorId0<anchorCount; anchorId0++) {
            const AnchorId anchorId1 = anchors.readFollowing(journeys, anchorId0, direction, minEdgeCoverage, aDrift, bDrift);
            if(anchorId1 == invalid<AnchorId>) {
                continue;
            }

            AnchorId anchorIdA = anchorId0;
            AnchorId anchorIdB = anchorId1;
            if(direction == 1) {
                std::swap(anchorIdA, anchorIdB);
            }

            bool edgeExists = false;
            tie(ignore, edgeExists) = boost::edge(anchorIdA, anchorIdB, anchorGraph);
            if(edgeExists) {
                continue;
            }

            // Create an AnchorPair with the common oriented reads between
            // anchorIdA and anchorIdB.
            const AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);



            // If the AnchorPair is consistent, generate a single edge.
            if(anchorPair.isConsistent(anchors, aDrift, bDrift, positions, offsets)) {
                const uint64_t offset = anchorPair.getAverageOffset(anchors);
                edge_descriptor e;
                tie(e, ignore) = add_edge(anchorIdA, anchorIdB,
                            AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                anchorGraph[e].useForAssembly = true;
            }



            // If the AnchorPair is not consistent, split it and generate one or more edges.
            else {

                // We have to split this AnchorPair into consistent AnchorPairs,
                // then generate a new edge for each (a set of parallel edges).
                anchorPair.splitByOffsets(anchors, aDrift, bDrift, newAnchorPairs);

                for(const AnchorPair& anchorPair: newAnchorPairs) {
                    if(anchorPair.orientedReadIds.size() >= minEdgeCoverage) {
                        const uint64_t offset = anchorPair.getAverageOffset(anchors);
                        edge_descriptor e;
                        tie(e, ignore) = add_edge(anchorIdA, anchorPair.anchorIdB,
                            AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                        anchorGraph[e].useForAssembly = true;
                        anchorGraph[e].isParallelEdge = true;
                    }
                }

            }
        }
    }

}



void AnchorGraph::createEdges2(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    double aDrift,
    double bDrift)
{
    AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = anchors.size();

    // Create the edges using read following.
    AnchorPair anchorPair;
    uint32_t offset;
    for(uint64_t direction=0; direction<2; direction++) {
        for(AnchorId anchorId0=0; anchorId0<anchorCount; anchorId0++) {

            // Do read following.
            // If successful, this finds a consistent AnchorPair
            // so we don't have to worry about splitting it.
            const bool found = anchors.readFollowing(
                journeys, anchorId0, direction,
                minEdgeCoverage, aDrift, bDrift,
                anchorPair, offset);

            // Not successful, do nothing.
            if(not found) {
                continue;
            }

            // If the edge already exists, do nothing.
            bool edgeExists = false;
            tie(ignore, edgeExists) = boost::edge(anchorPair.anchorIdA, anchorPair.anchorIdB, anchorGraph);
            if(edgeExists) {
                continue;
            }

            // We don't already have this edge. Add it.
            edge_descriptor e;
            tie(e, ignore) = add_edge(anchorPair.anchorIdA, anchorPair.anchorIdB,
                AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
            anchorGraph[e].useForAssembly = true;
        }
    }
}



void AnchorGraph::createEdges3(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    uint64_t minContinueReadFollowingCount,
    double aDrift,
    double bDrift)
{
    AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = anchors.size();

    // Create the edges using read following.
    vector<pair<AnchorPair, uint32_t> > anchorPairs;    // AnchorPairs with offsets.
    for(uint64_t direction=0; direction<2; direction++) {
        for(AnchorId anchorId0=0; anchorId0<anchorCount; anchorId0++) {

            // Do read following.
            anchors.readFollowing(
                journeys, anchorId0, direction,
                minEdgeCoverage, minContinueReadFollowingCount, aDrift, bDrift,
                anchorPairs);

            // Loop over the AnchorPairs that were found.
            for(const auto& p: anchorPairs) {
                const AnchorPair& anchorPair = p.first;
                const uint32_t offset = p.second;

                // If the edge already exists, do nothing.
                bool edgeExists = false;
                tie(ignore, edgeExists) = boost::edge(anchorPair.anchorIdA, anchorPair.anchorIdB, anchorGraph);
                if(edgeExists) {
                    continue;
                }

                // We don't already have this edge. Add it.
                edge_descriptor e;
                tie(e, ignore) = add_edge(anchorPair.anchorIdA, anchorPair.anchorIdB,
                    AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                anchorGraph[e].useForAssembly = true;
            }
        }
    }
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



// Dijkstra search.
// This performs a shortest path search starting at the specified AnchorId
// and stops when it finds a consistent AnchorPair that
// satisfies the coverage criteria and can therefore be used to create a new edge.
// If direction is 0, this returns an AnchorPair with AnchorIdA = startAnchorId;
// If direction is 1, this returns an AnchorPair with AnchorIdB = startAnchorId;
bool AnchorGraph::search(
    uint64_t direction,
    AnchorId startAnchorId,
    const Anchors& anchors,
    const ReadLengthDistribution&  readLengthDistribution,
    double aDrift,
    double bDrift,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    uint64_t maxDistance,
    AnchorPair& anchorPairArgument,
    uint64_t& offsetArgument
    ) const
{
    if(direction == 0) {
        return searchForward(
            startAnchorId,
            anchors,
            readLengthDistribution,
            aDrift, bDrift,
            minEdgeCoverageNear, minEdgeCoverageFar,
            maxDistance,
            anchorPairArgument, offsetArgument
            );
    } else if(direction == 1) {
        return searchBackward(
            startAnchorId,
            anchors,
            readLengthDistribution,
            aDrift, bDrift,
            minEdgeCoverageNear, minEdgeCoverageFar,
            maxDistance,
            anchorPairArgument, offsetArgument
            );
    } else {
        SHASTA2_ASSERT(0);
    }
}



bool AnchorGraph::searchForward(
    AnchorId startAnchorId,
    const Anchors& anchors,
    const ReadLengthDistribution&  readLengthDistribution,
    double aDrift,
    double bDrift,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    uint64_t maxDistance,
    AnchorPair& anchorPairArgument,
    uint64_t& offsetArgument
    ) const
{
    using boost::multi_index_container;
    using boost::multi_index::indexed_by;
    using boost::multi_index::member;
    using boost::multi_index::ordered_unique;
    using boost::multi_index::ordered_non_unique;

    const AnchorGraph& anchorGraph = *this;

    const bool debug = false;

    if(debug) {
        cout << "AnchorGraph::search called for " << anchorIdToString(startAnchorId) << endl;
    }



    // A container with the vertices (AnchorIds) encountered so far and their distance in the search.
    class VertexInfo {
    public:
        AnchorId anchorId;
        uint64_t distance;
    };
    using VerticesEncountered = multi_index_container<VertexInfo, indexed_by<
            ordered_unique<member <VertexInfo, AnchorId, &VertexInfo::anchorId> >,
            ordered_non_unique< member<VertexInfo, uint64_t ,&VertexInfo::distance> >
        > >;

    // Nomenclature consistent with the Wikipedia article
    // https://en.wikipedia.org/wiki/Dijkstra's_algorithm#Algorithm
    VerticesEncountered unvisited;
    auto& indexByAnchorId = unvisited.get<0>();
    auto& indexByDistance = unvisited.get<1>();

    std::set<AnchorId> visited;
    unvisited.insert({startAnchorId, 0});

    vector<AnchorPair> newAnchorPairs;

    // Main loop.
    while(not unvisited.empty()) {

        // Get the unvisited with the smallest distance.
        auto it0 = indexByDistance.begin();
        const VertexInfo& info0 = *it0;
        const AnchorId anchorId0 = info0.anchorId;
        const uint64_t distance0 = info0.distance;
        indexByDistance.erase(it0);

        if(debug) {
            cout << "Dequeued " << anchorIdToString(anchorId0) << " at distance " << distance0 << endl;
        }

        // See if this is the AnchorPair we want.
        if(anchorId0 != startAnchorId) {
            AnchorPair anchorPair(anchors, startAnchorId, anchorId0, false);

            // Split it to make it consistent.
            anchorPair.splitByOffsets(anchors, aDrift, bDrift, newAnchorPairs);
            const AnchorPair& consistentAnchorPair = newAnchorPairs.front();
            const uint64_t offset = consistentAnchorPair.getAverageOffset(anchors);
            const uint64_t coverage = consistentAnchorPair.orientedReadIds.size();

            if(coverage <= 1) {
                continue;
            }

            // Compute the edge coverage threshold that applies for this offset.
            const uint64_t bin = offset / readLengthDistribution.binWidth;
            const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
            uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
            edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);

            if(debug) {
                cout << "    Offset from start " << offset << ", coverage " << coverage <<
                    ", coverage threshold " << edgeCoverageThreshold << endl;
            }

            if(coverage >= edgeCoverageThreshold) {
                anchorPairArgument = anchorPair;
                offsetArgument = offset;
                return true;
            }

       }

        BGL_FORALL_OUTEDGES(anchorId0, e, anchorGraph, AnchorGraph) {
            const AnchorId anchorId1 = target(e, anchorGraph);
            const uint64_t distance1 = distance0 + anchorGraph[e].offset;

            if(distance1 > maxDistance) {
                continue;
            }

            if(debug) {
                cout << "    Found " << anchorIdToString(anchorId1) << " at distance " << distance1 << endl;
            }


            const auto it1 = indexByAnchorId.find(anchorId1);
            if(it1 == indexByAnchorId.end()) {
                unvisited.insert({anchorId1, distance1});

                if(debug) {
                    cout << "    Added " << anchorIdToString(anchorId1) << " at distance " << distance1 << endl;
                }
            } else {
                const uint64_t oldDistance1 = it1->distance;
                if(distance1 < oldDistance1) {
                    indexByAnchorId.replace(it1, {anchorId1, distance1});
                    if(debug) {
                        cout << "    Replaced " << anchorIdToString(anchorId1) <<
                            ", old distance " << oldDistance1 <<
                            ", new distance " << distance1 << endl;
                    }
                } else {
                    if(debug) {
                        cout << "    Did not replace " << anchorIdToString(anchorId1) <<
                            ", old distance " << oldDistance1 <<
                            ", new distance " << distance1 << endl;
                    }
                }
            }
        }



    }


    return false;
}



bool AnchorGraph::searchBackward(
    AnchorId startAnchorId,
    const Anchors& anchors,
    const ReadLengthDistribution&  readLengthDistribution,
    double aDrift,
    double bDrift,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    uint64_t maxDistance,
    AnchorPair& anchorPairArgument,
    uint64_t& offsetArgument
    ) const
{
    using boost::multi_index_container;
    using boost::multi_index::indexed_by;
    using boost::multi_index::member;
    using boost::multi_index::ordered_unique;
    using boost::multi_index::ordered_non_unique;

    const AnchorGraph& anchorGraph = *this;

    const bool debug = false; // anchorIdToString(startAnchorId) == "198725-";

    if(debug) {
        cout << "AnchorGraph::search called for " << anchorIdToString(startAnchorId) << endl;
    }



    // A container with the vertices (AnchorIds) encountered so far and their distance in the search.
    class VertexInfo {
    public:
        AnchorId anchorId;
        uint64_t distance;
    };
    using VerticesEncountered = multi_index_container<VertexInfo, indexed_by<
            ordered_unique<member <VertexInfo, AnchorId, &VertexInfo::anchorId> >,
            ordered_non_unique< member<VertexInfo, uint64_t ,&VertexInfo::distance> >
        > >;

    // Nomenclature consistent with the Wikipedia article
    // https://en.wikipedia.org/wiki/Dijkstra's_algorithm#Algorithm
    VerticesEncountered unvisited;
    auto& indexByAnchorId = unvisited.get<0>();
    auto& indexByDistance = unvisited.get<1>();

    std::set<AnchorId> visited;
    unvisited.insert({startAnchorId, 0});

    vector<AnchorPair> newAnchorPairs;

    // Main loop.
    while(not unvisited.empty()) {

        // Get the unvisited with the smallest distance.
        auto it0 = indexByDistance.begin();
        const VertexInfo& info0 = *it0;
        const AnchorId anchorId0 = info0.anchorId;
        const uint64_t distance0 = info0.distance;
        indexByDistance.erase(it0);

        if(debug) {
            cout << "Dequeued " << anchorIdToString(anchorId0) << " at distance " << distance0 << endl;
        }

        // See if this is the AnchorPair we want.
        if(anchorId0 != startAnchorId) {
            AnchorPair anchorPair(anchors, anchorId0, startAnchorId, false);

            // Split it to make it consistent.
            anchorPair.splitByOffsets(anchors, aDrift, bDrift, newAnchorPairs);
            const AnchorPair& consistentAnchorPair = newAnchorPairs.front();
            const uint64_t offset = consistentAnchorPair.getAverageOffset(anchors);
            const uint64_t coverage = consistentAnchorPair.orientedReadIds.size();
            if(coverage <= 1) {
                continue;
            }

            // Compute the edge coverage threshold that applies for this offset.
            const uint64_t bin = offset / readLengthDistribution.binWidth;
            const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
            uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
            edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);

            if(debug) {
                cout << "    Offset from start " << offset << ", coverage " << coverage <<
                    ", coverage threshold " << edgeCoverageThreshold << endl;
            }

            if(coverage >= edgeCoverageThreshold) {
                anchorPairArgument = anchorPair;
                offsetArgument = offset;
                return true;
            }

       }

        BGL_FORALL_INEDGES(anchorId0, e, anchorGraph, AnchorGraph) {
            const AnchorId anchorId1 = source(e, anchorGraph);
            const uint64_t distance1 = distance0 + anchorGraph[e].offset;
            if(distance1 > maxDistance) {
                continue;
            }



            const auto it1 = indexByAnchorId.find(anchorId1);
            if(it1 == indexByAnchorId.end()) {
                unvisited.insert({anchorId1, distance1});

                if(debug) {
                    cout << "    Added " << anchorIdToString(anchorId1) << " at distance " << distance1 << endl;
                }
            } else {
                const uint64_t oldDistance1 = it1->distance;
                if(distance1 < oldDistance1) {
                    indexByAnchorId.replace(it1, {anchorId1, distance1});
                    if(debug) {
                        cout << "    Replaced " << anchorIdToString(anchorId1) <<
                            ", old distance " << oldDistance1 <<
                            ", new distance " << distance1 << endl;
                    }
                } else {
                    if(debug) {
                        cout << "    Did not replace " << anchorIdToString(anchorId1) <<
                            ", old distance " << oldDistance1 <<
                            ", new distance " << distance1 << endl;
                    }
                }
            }
        }



    }


    return false;
}



// Eliminate dead ends where possible, using shortest path searches.
void AnchorGraph::handleDeadEnds(
    const Anchors& anchors,
    const ReadLengthDistribution& readLengthDistribution,
    double aDrift,
    double bDrift,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    uint64_t maxDistance,
    uint64_t threadCount)
{
    AnchorGraph& anchorGraph = *this;
    performanceLog << timestamp << "AnchorGraph::handleDeadEnds begins." << endl;

    // Fill in the HandleDeadEndsData that need to be visible by all threads.
    HandleDeadEndsData& data = handleDeadEndsData;
    data.anchors = &anchors;
    data.readLengthDistribution = &readLengthDistribution;
    data.aDrift = aDrift;
    data.bDrift = bDrift;
    data.minEdgeCoverageNear = minEdgeCoverageNear;
    data.minEdgeCoverageFar = minEdgeCoverageFar;
    data.maxDistance = maxDistance;

    // Create a filtered AnchorGraph containing only the edges marked as "useForAssembly".
    class EdgePredicate {
    public:
        bool operator()(const AnchorGraph::edge_descriptor& e) const
        {
            return (*anchorGraph)[e].useForAssembly;
        }
        EdgePredicate(const AnchorGraph& anchorGraph) : anchorGraph(&anchorGraph) {}
        EdgePredicate() : anchorGraph(0) {}
        const AnchorGraph* anchorGraph;
    };
    using FilteredAnchorGraph = boost::filtered_graph<AnchorGraph, EdgePredicate>;
    FilteredAnchorGraph filteredAnchorGraph(anchorGraph, EdgePredicate(anchorGraph));

    // Fill in the dead ends in the filtered graph. Ignore isolated vertices.
    // Probably we should ignore all vertices in very small connected components instead.
    BGL_FORALL_VERTICES(anchorId0, filteredAnchorGraph, FilteredAnchorGraph) {
        const bool isForwardDeadEnd = (out_degree(anchorId0, filteredAnchorGraph) == 0);
        const bool isBackwardDeadEnd = (in_degree(anchorId0, filteredAnchorGraph) == 0);

        // If it is isolated, ignore it.
        if(isForwardDeadEnd and isBackwardDeadEnd) {
            continue;
        }

        if(isForwardDeadEnd) {
            data.deadEnds.push_back(make_pair(anchorId0, 0));
        } else if(isBackwardDeadEnd) {
            data.deadEnds.push_back(make_pair(anchorId0, 1));
        }
    }
    cout << "Found " << data.deadEnds.size() << " dead ends in the initial AnchorGraph." << endl;

    // Run the searches in parallel.
    data.threadPairs.clear();
    data.threadPairs.resize(threadCount);
    const uint64_t batchCount = 1;
    setupLoadBalancing(data.deadEnds.size(), batchCount);
    runThreads(&AnchorGraph::handleDeadEndsThreadFunction, threadCount);

    // For reproducibility, consolidate the pairs found by each thread and deduplicate.
    vector< pair<AnchorPair, uint64_t> > allPairs;
    for(const auto& v: data.threadPairs) {
        copy(v.begin(), v.end(), back_inserter(allPairs));
    }
    data.threadPairs.clear();
    class SortHelper{
    public:
        bool operator()(const pair<AnchorPair, uint64_t>& x, const pair<AnchorPair, uint64_t>& y) const
        {
            return tie(x.first.anchorIdA, x.first.anchorIdB) < tie(y.first.anchorIdA, y.first.anchorIdB);
        }
    };
    sort(allPairs.begin(), allPairs.end(), SortHelper());
    std::unique(allPairs.begin(), allPairs.end(), SortHelper());

    cout << "Found " << allPairs.size() << " candidate edges for " <<
        data.deadEnds.size() << " dead ends." << endl;


    // Add one edge for each AnchorPair, but only if at least one of the two anchors
    // is still a dead end.
    uint64_t addedCount = 0;
    for(const auto& p: allPairs) {
        const AnchorPair& anchorPair = p.first;
        const uint64_t offset = p.second;
        const AnchorId anchorIdA = anchorPair.anchorIdA;
        const AnchorId anchorIdB = anchorPair.anchorIdB;

        if((out_degree(anchorIdA, filteredAnchorGraph) == 0) or (in_degree(anchorIdB, filteredAnchorGraph) == 0)) {
            edge_descriptor e;
            tie(e, ignore) = add_edge(anchorIdA, anchorIdB,
                AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
            ++addedCount;
            anchorGraph[e].useForAssembly = true;
            anchorGraph[e].addedAtDeadEnd = true;
        }
    }
    cout << "Added " << addedCount << " new edges to fix dead ends." << endl;

}



void AnchorGraph::handleDeadEndsThreadFunction(uint64_t threadId)
{
    HandleDeadEndsData& data = handleDeadEndsData;

    // Get the search parameters
    const Anchors& anchors = *data.anchors;
    const ReadLengthDistribution& readLengthDistribution = *data.readLengthDistribution;
    const double aDrift = data.aDrift;
    const double bDrift = data.bDrift;
    const uint64_t minEdgeCoverageNear = data.minEdgeCoverageNear;
    const uint64_t minEdgeCoverageFar = data.minEdgeCoverageFar;
    const uint64_t maxDistance = data.maxDistance;

    // Vector where we will store the AnchorPairs found by this thread.
    auto& threadPairs = data.threadPairs[threadId];
    threadPairs.clear();

    // Loop over batches of dead ends assigned to this thread.
    uint64_t begin, end;
    AnchorPair anchorPair;
    while(getNextBatch(begin, end)) {

        // Loop over dead ends in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const auto& p = data.deadEnds[i];
            const AnchorId anchorId = p.first;
            const uint64_t direction = p.second;

            uint64_t offset;
            const bool found = search(direction, anchorId, anchors, readLengthDistribution,
                aDrift, bDrift,
                minEdgeCoverageNear, minEdgeCoverageFar, maxDistance,
                anchorPair, offset);

            // If successful, store it.
            if(found) {
                threadPairs.emplace_back(anchorPair, offset);
            }
        }
    }


#if 0
    // Loop over forward or backward dead ends in the filtered anchor graph,
    // ignoring isolated vertices.
    // Probably we should ignore all vertices in very small connected components instead.
    BGL_FORALL_VERTICES(anchorId0, filteredAnchorGraph, FilteredAnchorGraph) {
        const bool isForwardDeadEnd = (out_degree(anchorId0, filteredAnchorGraph) == 0);
        const bool isBackwardDeadEnd = (in_degree(anchorId0, filteredAnchorGraph) == 0);

        // If it is isolated, ignore it.
        if(isForwardDeadEnd and isBackwardDeadEnd) {
            continue;
        }


        // Forward dead end.
        if(isForwardDeadEnd) {
            if(debug) {
                cout << timestamp << "Working on forward dead end at " << anchorIdToString(anchorId0) << endl;
            }

            AnchorPair anchorPair;
            uint64_t offset;
            const bool found = searchForward(anchorId0, anchors, readLengthDistribution,
                aDrift, bDrift,
                minEdgeCoverageNear, minEdgeCoverageFar, maxDistance,
                anchorPair, offset);

            if(debug) {
                if(found) {
                    cout << "Success. Found " << anchorIdToString(anchorPair.anchorIdB) <<
                        ", coverage " << anchorPair.orientedReadIds.size() <<
                        ", offset " << offset << endl;
                } else {
                    cout << "Failure." << endl;
                }
            }

        }
    }
#endif
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

