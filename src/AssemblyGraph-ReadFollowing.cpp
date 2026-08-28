#include "AssemblyGraph.hpp"
#include "ReadFollowing5.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;



bool AssemblyGraph::readFollowingStrandSymmetric(
    [[maybe_unused]] uint64_t tangleId,
    [[maybe_unused]] const Tangle& tangle,
    [[maybe_unused]] ostream& html)
{

#if 0
    const Graph graph(*this, tangleId, tangle, html);
#endif

    return false;
}
