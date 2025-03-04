#pragma once

/*******************************************************************************

Among all 4^k k-mers of length k, we choose a subset that we call "markers".
The markers are selected at the beginning of an assembly
and never changed, and selected in such a way that,
if (and only if) a k-mer is a marker, its reverse complement
is also a marker.

*******************************************************************************/

#include "Kmer.hpp"
#include "ReadId.hpp"
#include "Uint.hpp"

#include "utility.hpp"

namespace shasta {

    class CompressedMarker;
}



// Markers in shared memory are stored using class CompressedMarker.
class shasta::CompressedMarker {
public:

    // The position of this marker in the oriented read.
    // This limits the length of a read to 2^24=16Mib bases.
    Uint24 position;
};


