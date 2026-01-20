#pragma once

#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class AlignedBase;

    void abpoa(

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
        vector<AlignedBase>& alignedConsensus,

        // Consensus is always computed.
        // Alignment and alignedConsensus are only computed if this is set to true.
        bool computeAlignment

    );


    // The version with weights accepts a weight for each input sequence.
    // This can improve performance when the same sequence appears
    // more than once in a multiple sequence alignment.
    void abpoa(

        // The input sequences to be aligned.
        // They are presented to abpoa in this order.
        const vector< pair<vector<Base>, uint64_t> >& sequencesWithWeights,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The computed alignment.
        // Each element of the vector correspond to one of the input sequences,
        // in the same order.
        // These all have the same length, which equals the length of the aligned consensus.
        vector< vector<AlignedBase> >& alignment,

        // The aligned consensus.
        vector<AlignedBase>& alignedConsensus);



    void testAbpoa();
    void testAbpoaWithWeights();

}
