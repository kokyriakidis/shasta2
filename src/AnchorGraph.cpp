#include "AnchorGraph.hpp"
#include "Anchor.hpp"
using namespace shasta;


AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage) :
    AnchorGraphBaseClass(anchors.size())
{
    const uint64_t anchorCount = anchors.size();
    vector<AnchorPair> anchorPairs;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        AnchorPair::createChildren(anchors, journeys, anchorIdA, minEdgeCoverage, anchorPairs);

        for(const AnchorPair& anchorPair: anchorPairs) {
            add_edge(anchorIdA, anchorPair.anchorIdB, anchorPair, *this);
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}
