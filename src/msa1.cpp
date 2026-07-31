// Shasta2.
#include "msa1.hpp"
using namespace shasta2;

// Standard library.
#include "stdexcept.hpp"



// All arguments have the [[maybe_unused]]  attribute to prevent compiler
// warnings before implementation. They should be removed once
// the code is complete.
void shasta2::msa1(

    // The input sequences fixed on both sides, with their coverage.
    [[maybe_unused]] const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

    // The input sequences fixed on the left only, with their coverage.
    [[maybe_unused]] const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

    // The input sequences fixed on the right only, with their coverage.
    [[maybe_unused]] const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

    // The consensus sequence and its coverage.
    [[maybe_unused]] vector< pair<Base, uint64_t> >& consensus,

    // The computed alignment.
    // Each element of the vector correspond to one of the input sequences,
    // in the same order.
    // These all have the same length, which equals the length of the aligned consensus.
    [[maybe_unused]] vector< vector<AlignedBase> >& alignment,

    // The aligned consensus.
    [[maybe_unused]] vector<AlignedBase>& alignedConsensus,

    // Consensus and alignedConsensus are always computed.
    // Alignment is only computed if this set to true.
    [[maybe_unused]] bool computeAlignment
)
{
    throw runtime_error("Msa1 is not implemented.");
}
