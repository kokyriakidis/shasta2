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
#include <boost/pending/disjoint_sets.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AnchorGraph>;




AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage) :
    AnchorGraphBaseClass(anchors.size()),
    MappedMemoryOwner(anchors),
    MultithreadedObject<AnchorGraph>(*this)
{
    AnchorGraph& anchorGraph = *this;

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
	double coverageFractionThreshold,
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



    // Mark as "lowCoverage1" edges which don't pass a first coverage test.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        AnchorGraphEdge& edge = anchorGraph[e];
        const uint64_t offset = edge.offset;

        // Compute the edge coverage threshold that applies for this offset.
        const uint64_t bin = offset / readLengthDistribution.binWidth;
        const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
        uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
        edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);

        // If coverage is low, mark the edge as "lowCoverage1".
        if(edge.anchorPair.orientedReadIds.size() < edgeCoverageThreshold) {
            edge.lowCoverage1 = true;
        }
    }


    const bool debug = false;

    // Mark as "lowCoverage2" edges which don't pass a second coverage test.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {

    	if(debug) {
            const vertex_descriptor v0 = source(e, anchorGraph);
            const vertex_descriptor v1 = target(e, anchorGraph);
            cout << "Working on edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) << endl;
    	}

    	// Gather information about this edge.
        AnchorGraphEdge& edge = anchorGraph[e];
        if(edge.lowCoverage1) {
        	if(debug) {
        		cout << "Already flagged as lowCoverage1." << endl;
        	}
        	continue;
        }
        const uint64_t offset = edge.offset;
        const uint64_t coverage = edge.anchorPair.orientedReadIds.size();

        // Compute the coverage correlation factor for this offset.
        const uint64_t bin = offset / readLengthDistribution.binWidth;
        const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;

        // Compute maximum and total coverage for the outgoing edges of the source vertex of e.
        // Don't consider edges flagged as "lowCoverage1".
        const vertex_descriptor v0 = source(e, anchorGraph);
        uint64_t maxOutCoverage0 = 0;
        uint64_t totalOutCoverage0 = 0;
        BGL_FORALL_OUTEDGES(v0, e0, anchorGraph, AnchorGraph) {
        	const AnchorGraphEdge& edge0 = anchorGraph[e0];
        	if(edge0.lowCoverage1) {
        		continue;
        	}
        	const uint64_t outCoverage0 = edge0.anchorPair.orientedReadIds.size();
        	maxOutCoverage0 = max(maxOutCoverage0, outCoverage0);
        	totalOutCoverage0 += outCoverage0;
        }

        // If this edge has maximum coverage among the outgoing edges of v0,
        // don't flag it as "lowCoverage2".
        if(coverage == maxOutCoverage0) {
        	if(debug) {
        		cout << "Has maximum out-coverage." << endl;
        	}
        	continue;
        }

        // Compute maximum and total coverage for the incoming edges of the target vertex of e.
        // Don't consider edges flagged as "lowCoverage1".
        const vertex_descriptor v1 = target(e, anchorGraph);
        uint64_t maxInCoverage1 = 0;
        uint64_t totalInCoverage1 = 0;
        BGL_FORALL_INEDGES(v1, e1, anchorGraph, AnchorGraph) {
        	const AnchorGraphEdge& edge1 = anchorGraph[e1];
        	if(edge1.lowCoverage1) {
        		continue;
        	}
        	const uint64_t inCoverage1 = edge1.anchorPair.orientedReadIds.size();
        	maxInCoverage1 = max(maxInCoverage1, inCoverage1);
        	totalInCoverage1 += inCoverage1;
        }

        // If this edge has maximum coverage among the incoming edges of v1,
        // don't flag it as "lowCoverage2".
        if(coverage == maxInCoverage1) {
        	if(debug) {
        		cout << "Has maximum in-coverage." << endl;
        	}
        	continue;
        }

        // This is not the edge with the most out-coverage of v0
        // or with the most in-coverage of v1. We will apply an additional
        // coverage test.

        // Compute the additional coverage threshold we are going to use for this edge.
        const uint64_t minCoverage2 = uint64_t(std::round(
        		coverageFractionThreshold *
				double(min(totalOutCoverage0, totalInCoverage1)) *
				coverageCorrelation
				));

        // If coverage is low, flag this edge as "lowCoverage2".
        if(coverage < minCoverage2) {
        	edge.lowCoverage2 = true;
        	if(debug) {
        		cout << "Flagged as lowCoverage2." << endl;
        	}
        }

    }


    // Now mark as "useForAssembly" the edges not flagged as "lowCoverage1" or "lowCoverage2".
    uint64_t useForAssemblyCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
    	AnchorGraphEdge& edge = anchorGraph[e];
    	if(not (edge.lowCoverage1 or edge.lowCoverage2)) {
    		edge.useForAssembly = true;
    		++useForAssemblyCount;
    	}
    	if(debug) {
            const vertex_descriptor v0 = source(e, anchorGraph);
            const vertex_descriptor v1 = target(e, anchorGraph);
            cout << "Flags for edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) << " " <<
            	edge.lowCoverage1 << " " << edge.lowCoverage2 << " " << edge.useForAssembly << endl;
    	}
    }
    cout << useForAssemblyCount << " AnchorGraph edges were marked to be used for assembly." << endl;



    // Compute connected components using only the edges flagged as "useForAssembly".
    vector<uint64_t> rank(anchorCount);
    vector<uint64_t> parent(anchorCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<anchorCount; i++) {
        disjointSets.make_set(i);
    }
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(anchorGraph[e].useForAssembly) {
            const AnchorId anchorId0 = source(e, anchorGraph);
            const AnchorId anchorId1 = target(e, anchorGraph);
            disjointSets.union_set(anchorId0, anchorId1);
        }
    }

    // Count the vertices in each connected component.
    vector<uint64_t> componentSize(anchorCount, 0);
    BGL_FORALL_VERTICES(anchorId, anchorGraph, AnchorGraph) {
        const uint64_t componentId = disjointSets.find_set(anchorId);
        ++componentSize[componentId];
    }

    // Create a histogram of component sizes.
    vector<uint64_t> histogram;
    for(uint64_t componentId=0; componentId<anchorCount; componentId++) {
        const uint64_t size = componentSize[componentId];
        if(size > 0) {
            if(size >= histogram.size()) {
                histogram.resize(size + 1, 0);
            }
            ++histogram[size];
        }
    }

    /*
    cout << "Component size histogram for the AnchorGraph before handling dead ends:" << endl;
    for(uint64_t size=1; size<histogram.size(); size++) {
        const uint64_t frequency = histogram[size];
        if(frequency) {
            cout << size << "," << frequency << endl;
        }
    }
    */

    // Mark as not to be used for assembly edges in small connected components.
    uint64_t removedCount = 0;
    const uint64_t minComponentSize = 5;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(anchorGraph[e].useForAssembly) {
            const AnchorId anchorId0 = source(e, anchorGraph);
            const AnchorId anchorId1 = target(e, anchorGraph);
            const uint64_t componentId0 = disjointSets.find_set(anchorId0);
            const uint64_t componentId1 = disjointSets.find_set(anchorId1);
            SHASTA_ASSERT(componentId0 == componentId1);
            const uint64_t size = componentSize[componentId0];
            if(size < minComponentSize) {
                anchorGraph[e].useForAssembly = false;
                anchorGraph[e].inSmallComponent = true;
                ++removedCount;
            }
        }
    }
    cout << "Found " << removedCount << " edges in small connected components." << endl;




#if 0
    const uint64_t maxDistance = 300000;

    // Test search.
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    const bool success = search(0, anchorIdFromString("2149-"),
        anchors,
        readLengthDistribution,
        aDrift, bDrift,
        minEdgeCoverageNear, minEdgeCoverageFar, maxDistance,
        anchorPair, offset);
    if(success) {
        cout << "Found " << anchorIdToString(anchorPair.anchorIdA) << " -> " <<
            anchorIdToString(anchorPair.anchorIdB) <<
            ", coverage " << anchorPair.orientedReadIds.size() <<
            ", offset " << offset << endl;
    } else {
        cout << "Search failed." << endl;
    }
    return;

    // Eliminate dead ends where possible, using shortest path searches.
    handleDeadEnds(anchors, readLengthDistribution,
        aDrift, bDrift,
        minEdgeCoverageNear, minEdgeCoverageFar, maxDistance, threadCount);
#endif
}



// Constructor from binary data.
AnchorGraph::AnchorGraph(const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<AnchorGraph>(*this)
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
            anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);
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
            anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);
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
