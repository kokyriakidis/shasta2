#pragma once

#include "Superbubble.hpp"

namespace shasta2 {
    class SuperbubbleChain;
}



class shasta2::SuperbubbleChain: public vector<Superbubble> {
public:
    using vertex_descriptor = AssemblyGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AssemblyGraphBaseClass::edge_descriptor;

    uint64_t phase1(
        AssemblyGraph&,
        uint64_t superbubbleChainId) const;

private:

};
