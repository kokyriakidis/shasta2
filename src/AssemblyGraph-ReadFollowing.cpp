#include "AssemblyGraph.hpp"
#include "ReadFollowing5.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;



bool AssemblyGraph::readFollowingStrandSymmetric(
    [[maybe_unused]] uint64_t tangleId,
    [[maybe_unused]] const Tangle& tangle)
{
    /*
    if(not tangle.isSelfComplementary()) {
        const Graph graph(*this, tangleId, tangle);
    }
    */
    return false;
}
