#pragma once

/*****************************************************************

Two Base sequences are defined to be "similar" if they are likely to be
identical except for sequencing errors, that is,
if they differ only by copy numbers in repeats of short period.

areSimilarSequences returns true if its input sequences are "similar".
The precise criterion for similarity is defined by the minRepeatCount vector.

minRepeatCount[0] is ignored.

minRepeatCount[1] is the minimum homopolymer length for differences
in the length of a homopolymer run. If the length difference occurs
in a homopolymer run of at least this length,
it is considered likely to be caused by sequencing errors and
does not cause the two sequences to be considered "similar".

Similarly, minRepeatCount[p] defines the minimum length
(number of repeat units of period p) for a repeat with
with period p. If the copy number difference occurs in a repeat with p or
more units, it is considered likely to be caused by sequencing errors and
does not cause the two sequences to be considered "similar".
For example, if minRepeatCount[2] is 3, differences in lengths
of repeats with period 2 a longer than 3 units (6 bases)
do not cause the two sequences to be considered similar.

Note that if any mismatches are present, the two sequences are never considered to be similar.

areSimilarSequences uses abpoa to compute an alignment between
its input sequences. If the alignment contains any mismatches,
it returns false. Otherwise, it checks all indels
and returns false if indels that don't satisfy the above criteria
are found.

If html is open, details of the computation are output there.

*****************************************************************/

#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;

    bool areSimilarSequences(
        const vector<Base>&,
        const vector<Base>&,
        const vector<uint64_t>& minRepeatCount,
        ostream& html);
}
