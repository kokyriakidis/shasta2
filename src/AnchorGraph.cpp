// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
#include "ReadLengthDistribution.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include "boost/graph/iteration_macros.hpp"
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
// with inconsistent offsets.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    double aDrift,
    double bDrift) :
    MappedMemoryOwner(anchors)
{
    // cout << "AAA " << minEdgeCoverageNear << " " << minEdgeCoverageFar << endl;
    AnchorGraph& anchorGraph = *this;

    // The edge coverage threshold is a function of estimated offset.
    // It is computed as max(edgeCoverageThresholdNear * coverageCorrelation(offset), edgeCoverageThresholdFar).
    // As small offset coverageCorrelation is 1, so for small offsets the
    // edge coverage threshold is edgeCoverageThresholdNear.
    // At large offset the edge coverage threshold is edgeCoverageThresholdFar.

    ReadLengthDistribution readLengthDistribution(anchors);

    uint64_t nextEdgeId = 0;

    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    vector<AnchorPair> anchorPairs;
    vector<AnchorPair> newAnchorPairs;

    vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
    vector<uint64_t> offsets;

    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, minEdgeCoverageFar, anchorPairs);

        // cout << "BBB " << anchorIdToString(anchorIdA) << " found " << anchorPairs.size() << " children." << endl;

        for(const AnchorPair& anchorPair: anchorPairs) {

            if(anchorPair.isConsistent(anchors, aDrift, bDrift, positions, offsets)) {

                // Compute the edge coverage threshold for this AnchorPair.
                // It is a function of the offset.
                const uint64_t offset = anchorPair.getAverageOffset(anchors);
                const uint64_t bin = offset / readLengthDistribution.binWidth;
                const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
                uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
                edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);
                // cout << "CCC " << offset << " " << edgeCoverageThreshold << " " << anchorPair.orientedReadIds.size() << endl;

                // Just generate an edge with this AnchorPair.
                // This is the most common case.
                if(anchorPair.orientedReadIds.size() >= edgeCoverageThreshold) {
                    const uint64_t offset = anchorPair.getAverageOffset(anchors);
                    add_edge(anchorIdA, anchorPair.anchorIdB,
                        AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                }

            } else {

                // We have to split this AnchorPair into consistent AnchorPairs,
                // then generate a new edge for each (a set of parallel edges).
                anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);

                for(const AnchorPair& anchorPair: newAnchorPairs) {

                    // Compute the edge coverage threshold for this AnchorPair.
                    // It is a function of the offset.
                    const uint64_t offset = anchorPair.getAverageOffset(anchors);
                    const uint64_t bin = offset / readLengthDistribution.binWidth;
                    const double coverageCorrelation = readLengthDistribution.data[bin].coverageCorrelation;
                    uint64_t edgeCoverageThreshold = uint64_t(coverageCorrelation * double(minEdgeCoverageNear));
                    edgeCoverageThreshold = max(edgeCoverageThreshold, minEdgeCoverageFar);
                    // cout << "DDD " << offset << " " << edgeCoverageThreshold << " " << anchorPair.orientedReadIds.size() << endl;

                    if(anchorPair.size() >= edgeCoverageThreshold) {
                        const uint64_t offset = anchorPair.getAverageOffset(anchors);
                        add_edge(anchorIdA, anchorPair.anchorIdB,
                            AnchorGraphEdge(anchorPair, offset, nextEdgeId++), anchorGraph);
                    }
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

}



// Constructor that splits edges that have an AnchorPair
// with inconsistent offsets, and also does local search to
// eliminate dead ends where possible.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverageNear,
    uint64_t minEdgeCoverageFar,
    double aDrift,
    double bDrift,
    FixDeadEnds) :
    MappedMemoryOwner(anchors)
{
    AnchorGraph& anchorGraph = *this;

    // The edge coverage threshold is a function of estimated offset.
    // It is computed as max(edgeCoverageThresholdNear * coverageCorrelation(offset), edgeCoverageThresholdFar).
    // As small offset coverageCorrelation is 1, so for small offsets the
    // edge coverage threshold is edgeCoverageThresholdNear.
    // At large offset the edge coverage threshold is edgeCoverageThresholdFar.
    // The coverageCorrelation is provided by the ReadLengthDistribution;
    ReadLengthDistribution readLengthDistribution(anchors);

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
