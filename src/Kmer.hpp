#pragma once

#include "ShortBaseSequence.hpp"

namespace shasta {

    // Type used to represent a k-mer.
    // This limits the maximum k-mer length that can be used.
    using Kmer128 = ShortBaseSequence128;
    using Kmer = Kmer128;
    static_assert(Kmer::capacity == 128, "Unexpected Kmer capacity.");

}


