#ifndef SHASTA_READ_FLAGS_HPP
#define SHASTA_READ_FLAGS_HPP

#include "cstdint.hpp"

namespace shasta {
    class ReadFlags;
}

class shasta::ReadFlags {
public:

    // Set if we have other reads with the same name.
    uint8_t isDuplicate : 1;

    // Set if this read is not to be used in the assembly
    // due to the presence of duplicates.
    // The way this is set is determined by the value of --Reads.handleDuplicates:
    uint8_t discardDueToDuplicates : 1;

    // This is set for reads that are approximate palindromic,
    // that is, are well aligned with their own reverse complement.
    uint8_t isPalindromic : 1;

    // Set if the read is marked as chimeric.
    uint8_t isChimeric : 1;

    // The strand that this read will be assembled on.
    // Only used by Mode 2 assembly.
    // Set in flagCrossStrandReadGraphEdges2.
    uint8_t strand : 1;

    ReadFlags()
    {
        static_assert(sizeof(ReadFlags) == 1, "Unexpected size of ReadFlags.");
        *reinterpret_cast<uint8_t*>(this) = 0;
    }
};

#endif
