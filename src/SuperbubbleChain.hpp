#pragma once

#include "Superbubble.hpp"

namespace shasta {
    class SuperbubbleChain;
}



class shasta::SuperbubbleChain: public vector<Superbubble> {
public:
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    void phase(AssemblyGraph&, uint64_t superbubbleChainId);

private:

};
