#pragma once

#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class AlignedBase;

    void poasta(

        // The input sequences to be aligned.
        // They are presented to abpoa in this order.
        const vector< vector<Base> >& sequences,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The computed alignment.
        // Each element of the vector correspond to one of the input sequences,
        // in the same order.
        // These all have the same length, which equals the length of the aligned consensus.
        vector< vector<AlignedBase> >& alignment,

        // The aligned consensus.
        vector<AlignedBase>& alignedConsensus
    );

    // This tests a direct call to poasta.
    void testPoasta1();

    // This tests a call to poasta using the wrapper.
    void testPoasta2();

}
