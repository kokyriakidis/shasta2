// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
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

