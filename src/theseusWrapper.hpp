#pragma once

#include "cstdint.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class AlignedBase;


    void theseus(

        // The input sequences fixed on both sides, with their coverage.
        // They are passed to theseus in this order.
        const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

        // The input sequences fixed on the left only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

        // The input sequences fixed on the right only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

        // The consensus sequence and its coverage.
        vector< pair<Base, uint64_t> >& consensus,

        // The computed alignment.
        // Each element of the vector correspond to one of the input sequences,
        // in the same order.
        // These all have the same length, which equals the length of the aligned consensus.
        vector< vector<AlignedBase> >& alignment,

        // The aligned consensus.
        vector<AlignedBase>& alignedConsensus,

        // Consensus and alignedConsensus are always computed.
        // Alignment is only computed if this set to true.
        bool computeAlignment
    );

    // This writes a file that can be used as input to pericles.
    void theseusWriteFile(

        // The input sequences fixed on both sides, with their coverage.
        // They are passed to theseus in this order.
        const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

        // The input sequences fixed on the left only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

        // The input sequences fixed on the right only, with their coverage.
        const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

        const string& fileName
    );

    void testTheseus();
}
