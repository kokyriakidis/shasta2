#pragma once

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"

namespace shasta {
    class MarkerInfo;
    class Markers;
}



// Class to describe a single marker of an oriented read.
class shasta::MarkerInfo {
public:
    OrientedReadId orientedReadId;
    uint32_t ordinal;
    MarkerInfo() {}
    MarkerInfo(OrientedReadId orientedReadId, uint32_t ordinal) :
        orientedReadId(orientedReadId), ordinal(ordinal) {}

    // Construct the reverse complement of a MarkerInfo.
    MarkerInfo reverseComplement(const Markers&) const;
};
