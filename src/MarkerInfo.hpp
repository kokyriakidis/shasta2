#pragma once

// Shasta.
#include "Kmer.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"

namespace shasta2 {
    class MarkerInfo;
    class Markers;
    class Reads;
}



// Class to describe a single marker of an oriented read.
class shasta2::MarkerInfo {
public:
    OrientedReadId orientedReadId;
    uint32_t ordinal;

    // The position of the middle of the marker relative to the beginning of the oriented read.
    // This equals the position of the first base of the marker plus k/2.
    uint32_t position;

    MarkerInfo() {}
    MarkerInfo(
        OrientedReadId orientedReadId,
        uint32_t ordinal,
        uint32_t position) :
        orientedReadId(orientedReadId),
        ordinal(ordinal),
        position(position)
    {}

    // Construct the reverse complement of a MarkerInfo.
    MarkerInfo reverseComplement(
        const Reads&,
        const Markers&) const;

    // Get the Kmer corresponding to this MarkerInfo.
    Kmer getKmer(uint64_t k, const Reads&) const;
};

