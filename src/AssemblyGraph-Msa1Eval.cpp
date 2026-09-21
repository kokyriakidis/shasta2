// Shasta2.
#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "algorithm.hpp"
#include "iterator.hpp"



// For the msa1 hard-region evaluation harness. Return, for every step of every
// edge, exactly the (anchorIdA, anchorIdB, orientedReadIds) that assembleStep
// would pass to LocalAssembly7 - including the additional oriented reads
// borrowed from the previous/next step. See assembleStep for the logic mirrored here.
vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string> > >
    AssemblyGraph::getAssemblyGraphSteps() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string> > > steps;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        for(uint64_t i=0; i<edge.size(); i++) {
            vector<OrientedReadId> additionalOrientedReadIds;
            if(i > 0) {
                std::ranges::copy(edge[i - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
            }
            if(i < edge.size() - 1) {
                std::ranges::copy(edge[i + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
            }
            deduplicate(additionalOrientedReadIds);

            const AnchorPair& anchorPair = edge[i].anchorPair;
            vector<OrientedReadId> orientedReadIds = additionalOrientedReadIds;
            std::ranges::copy(anchorPair.orientedReadIds, back_inserter(orientedReadIds));
            deduplicate(orientedReadIds);

            vector<string> orientedReadIdStrings;
            orientedReadIdStrings.reserve(orientedReadIds.size());
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                orientedReadIdStrings.push_back(orientedReadId.getString());
            }

            steps.push_back(make_tuple(
                uint64_t(edge.id), i,
                anchorPair.anchorIdA, anchorPair.anchorIdB,
                orientedReadIdStrings));
        }
    }

    return steps;
}
