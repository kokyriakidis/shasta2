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
        vector<Base>& consensus,
        vector<AlignedBase>& alignedConsensus,
        vector<uint64_t>& coverage,

        // The alignment is only computed if computeAlignment is true.
        vector< vector<AlignedBase> >& alignment,
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
