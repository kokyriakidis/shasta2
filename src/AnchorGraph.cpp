// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
#include "ReadLengthDistribution.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>



AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage) :
    AnchorGraphBaseClass(anchors.size()),
    MappedMemoryOwner(anchors)
{
    AnchorGraph& anchorGraph = *this;

    uint64_t nextEdgeId = 0;

    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    vector<AnchorPair> anchorPairs;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, minEdgeCoverage, anchorPairs);

        for(const AnchorPair& anchorPair: anchorPairs) {
            const uint64_t offset = anchorPair.getAverageOffset(anchors);
            add_edge(anchorIdA, anchorPair.anchorIdB,
                AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

}



// Constructor that splits edges that have an AnchorPair
// with inconsistent offsets, and also does local search to
// eliminate dead ends where possible.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const ReadLengthDistribution& readLengthDistribution,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    double aDrift,
    double bDrift) :
    MappedMemoryOwner(anchors)
{
    AnchorGraph& anchorGraph = *this;

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    // Initial creation of the edges, without any limits on coverage.
    // Edges that don't have a consistent offset are split as necessary.
    uint64_t nextEdgeId = 0;

    // Work vectors for the loop below.
    vector<AnchorPair> anchorPairs;
    vector<AnchorPair> newAnchorPairs;
    vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
    vector<uint64_t> offsets;

    // Loop over possible source vertices.
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {

        // Find next AnchorIds in journeys, without coverage limitations.
        AnchorPair::createChildren(anchors, journeys, anchorIdA, 0, anchorPairs);


        // Each AnchorPair generates an edge if it is consistent
        // or more than one if it is not. Here too we don't use any coverage limitation.
        for(const AnchorPair& anchorPair: anchorPairs) {

            if(anchorPair.isConsistent(anchors, aDrift, bDrift, positions, offsets)) {

                const uint64_t offset = anchorPair.getAverageOffset(anchors);
                add_edge(anchorIdA, anchorPair.anchorIdB,
                    AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);

            } else {

                // We have to split this AnchorPair into consistent AnchorPairs,
                // then generate a new edge for each (a set of parallel edges).
                anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);

                // cout << "Generating " << newAnchorPairs.size() << " parallel edges for " <<
                //    anchorIdToString(anchorIdA) << " -> " << anchorIdToString(anchorPair.anchorIdB) << endl;


                for(const AnchorPair& anchorPair: newAnchorPairs) {
                    const uint64_t offset = anchorPair.getAverageOffset(anchors);
                    edge_descriptor e;
                    tie(e, ignore) = add_edge(anchorIdA, anchorPair.anchorIdB,
                        AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                    anchorGraph[e].isParallelEdge = true;
                }
            }
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Now all the edges must have consistent offsets.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        SHASTA_ASSERT(anchorGraph[e].anchorPair.isConsistent(anchors, aDrift, bDrift, positions, offsets));
        SHASTA_ASSERT(anchorGraph[e].offset < 100000000);
    }


    // Mark as "useForAssembly" edges with sufficient coverage.
    uint64_t useForAssemblyCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        AnchorGraphEdge& edge = anchorGraph[e];
        const uint64_t offset = edge.offset;

        // Compute the edge coverage threshold that applies for this offset.
        const uint64_t bin = offset / readLengthDistribution.binWidth;
        const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
        uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
        edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);

        // If coverage is sufficient, mark the edge as useForAssembly.
        if(edge.anchorPair.orientedReadIds.size() >= edgeCoverageThreshold) {
            edge.useForAssembly = true;
            ++useForAssemblyCount;
        }
    }
    cout << useForAssemblyCount << " AnchorGraph edges were marked to be used for assembly." << endl;

#if 0
    // Test search.
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    search(1, anchorIdFromString("421976+"),
        anchors,
        readLengthDistribution,
        aDrift, bDrift,
        minEdgeCoverageNear, minEdgeCoverageFar,
        anchorPair, offset);
    cout << "Found " << anchorIdToString(anchorPair.anchorIdA) << " -> " <<
        anchorIdToString(anchorPair.anchorIdB) <<
        ", coverage " << anchorPair.orientedReadIds.size() <<
        ", offset " << offset << endl;
#endif
}



// Constructor from binary data.
AnchorGraph::AnchorGraph(const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    load();
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



void AnchorGraph::save() const
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
    const string name = largeDataName("AnchorGraph");
    MemoryMapped::Vector<char> data;
    data.createNew(name, largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AnchorGraph::load()
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        const string name = largeDataName("AnchorGraph");
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error("AnchorGraph is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error(string("Error reading AnchorGraph: ") + e.what());
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
    AnchorPair& anchorPairArgument,
    uint64_t& offsetArgument
    ) const
{
    if(direction == 0) {
        return searchForward(
            startAnchorId,
            anchors,
            readLengthDistribution,
            aDrift,
            bDrift,
            minEdgeCoverageNear,
            minEdgeCoverageFar,
            anchorPairArgument,
            offsetArgument
            );
    } else if(direction == 1) {
        return searchBackward(
            startAnchorId,
            anchors,
            readLengthDistribution,
            aDrift,
            bDrift,
            minEdgeCoverageNear,
            minEdgeCoverageFar,
            anchorPairArgument,
            offsetArgument
            );
    } else {
        SHASTA_ASSERT(0);
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
            anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);
            const AnchorPair& consistentAnchorPair = newAnchorPairs.front();
            const uint64_t offset = consistentAnchorPair.getAverageOffset(anchors);
            const uint64_t coverage = consistentAnchorPair.orientedReadIds.size();

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
            AnchorPair anchorPair(anchors, anchorId0, startAnchorId, false);

            // Split it to make it consistent.
            anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);
            const AnchorPair& consistentAnchorPair = newAnchorPairs.front();
            const uint64_t offset = consistentAnchorPair.getAverageOffset(anchors);
            const uint64_t coverage = consistentAnchorPair.orientedReadIds.size();

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
