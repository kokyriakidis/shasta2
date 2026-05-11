#pragma once

#include "cstdint.hpp"
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

        // The consensus sequence.
        vector<Base>& consensus,

        // The alignment is only computed if computeAlignment is true.
        vector< vector<AlignedBase> >& alignment,
        bool computeAlignment
    );

    void testTheseus();
}
