// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include "boost/graph/iteration_macros.hpp"

// Standard library.
#include "fstream.hpp"



AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage) :
    AnchorGraphBaseClass(anchors.size())
{
    AnchorGraph& anchorGraph = *this;

    const uint64_t anchorCount = anchors.size();
    vector<AnchorPair> anchorPairs;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, minEdgeCoverage, anchorPairs);

        for(const AnchorPair& anchorPair: anchorPairs) {
            add_edge(anchorIdA, anchorPair.anchorIdB, anchorPair, anchorGraph);
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
    uint64_t minEdgeCoverage,
    double aDrift,
    double bDrift)
{
    AnchorGraph& anchorGraph = *this;

    const uint64_t anchorCount = anchors.size();
    vector<AnchorPair> anchorPairs;
    vector<AnchorPair> newAnchorPairs;

    vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
    vector<uint64_t> offsets;

    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, minEdgeCoverage, anchorPairs);

        for(const AnchorPair& anchorPair: anchorPairs) {

            if(anchorPair.isConsistent(anchors, aDrift, bDrift, positions, offsets)) {

                // Just generate an edge with this AnchorPair.
                // This is the most common case.
                add_edge(anchorIdA, anchorPair.anchorIdB, anchorPair, anchorGraph);

            } else {

                // We have to split this AnchorPair into consistent AnchorPairs,
                // then generate a new edge for each (a set of parallel edges).
                anchorPair.split(anchors, aDrift, bDrift, newAnchorPairs);

                for(const AnchorPair& anchorPair: newAnchorPairs) {
                    if(anchorPair.size() >= minEdgeCoverage) {
                        add_edge(anchorIdA, anchorPair.anchorIdB, anchorPair, anchorGraph);
                    }
                }

            }
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Now all the edges must have consistent offsets.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        SHASTA_ASSERT(anchorGraph[e].isConsistent(anchors, aDrift, bDrift, positions, offsets));
    }

}



// Compute the edge journeys.
// The edge journey of an OrientedReadId is the sequence of
// AnchorGraph edges visited by the OrientedReadId.
// An OrientedReadId visits an AnchorGraph edges if it
// appears in the AnchorPairt for the edge.
// Edge journeys are indexed by OrientedReadId::getValue().
void AnchorGraph::computeEdgeJourneys(
    const Anchors& anchors,
    vector< vector<edge_descriptor> >& edgeJourneys
) const
{
    const AnchorGraph& anchorGraph = *this;
    const uint64_t orientedReadCount = anchors.markers.size();

    // Construct a vector of pairs (ordinal, edge_descriptor) for each OrientedReadId;
    vector< vector< pair<uint32_t, edge_descriptor> > > v(orientedReadCount);
    vector< pair<uint32_t, uint32_t> > ordinals;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const AnchorPair& anchorPair = anchorGraph[e];
        anchorPair.getOrdinals(anchors, ordinals);
        SHASTA_ASSERT(ordinals.size() == anchorPair.orientedReadIds.size());
        for(uint64_t i=0; i<anchorPair.orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
            const uint32_t ordinal0 = ordinals[i].first;
            v[orientedReadId.getValue()].push_back(make_pair(ordinal0, e));
        }
    }

    // Sort them by ordinals.
    for(vector< pair<uint32_t, edge_descriptor> >& u: v) {
        sort(u.begin(), u.end(), OrderPairsByFirstOnly<uint32_t, edge_descriptor>());
    }

    // Now we can create the edge journeys.
    edgeJourneys.clear();
    edgeJourneys.resize(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        vector<edge_descriptor> edgeJourney = edgeJourneys[i];
        const vector< pair<uint32_t, edge_descriptor> >& u = v[i];
        for(const auto& p: u) {
            edgeJourney.push_back(p.second);
        }
    }
 }

