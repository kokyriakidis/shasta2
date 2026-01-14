#pragma once

#include "Superbubble.hpp"

namespace shasta2 {
    class SuperbubbleChain;
}



class shasta2::SuperbubbleChain: public vector<Superbubble> {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    uint64_t phase1(
        AssemblyGraph&,
        uint64_t superbubbleChainId);

private:

};
