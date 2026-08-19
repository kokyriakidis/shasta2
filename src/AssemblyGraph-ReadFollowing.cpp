#include "AssemblyGraph.hpp"
#include "ReadFollowing5.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;



bool AssemblyGraph::readFollowingStrandSymmetric(uint64_t, const Tangle&)
{
    /*
    if(not tangle.isSelfComplementary()) {
        const Graph graph(*this, tangleId, tangle);
    }
    */
    return false;
}
