// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
using namespace shasta;

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
            add_edge(anchorIdA, anchorPair.anchorIdB, anchorPair, *this);
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Write out the isolated anchors.
    uint64_t isolatedAnchorCount = 0;
    ofstream csv("IsolatedAnchors.csv");
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        if(
            (in_degree (anchorId, anchorGraph) == 0) and
            (out_degree(anchorId, anchorGraph) == 0)
            ) {
            ++ isolatedAnchorCount;
            csv << anchorIdToString(anchorId) << "\n";
        }
    }
    cout << "The anchor graph has " << isolatedAnchorCount << " isolated vertices." << endl;
}
