#pragma once

#include "ShortBaseSequence.hpp"

namespace shasta {

    // Type used to represent a k-mer.
    // This limits the maximum k-mer length that can be used.
    using Kmer16 = ShortBaseSequence16;
    using Kmer32 = ShortBaseSequence32;
    using Kmer64 = ShortBaseSequence64;

    using Kmer128 = ShortBaseSequence128;
    using Kmer = Kmer128;
    static_assert(Kmer::capacity == 128, "Unexpected Kmer capacity.");

}


