#ifndef SHASTA_MARKER_ACCESS_FUNCTIONS_HPP

#include "Kmer.hpp"
#include "ReadId.hpp"

namespace shasta {

    class Markers;
    class Reads;

    // Access functions for markers Kmers and KmerIds.
    // There are similar member functions in class Assembler,
    // but these are accessible anywhere else.

    Kmer getOrientedReadMarkerKmer(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const Markers& markers
        );

    Kmer getOrientedReadMarkerKmerStrand0(
        ReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const Markers& markers
        );

    Kmer getOrientedReadMarkerKmerStrand1(
        ReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const Markers& markers
        );

    // Get the marker KmerId for an oriented read and ordinal.
    KmerId getOrientedReadMarkerKmerId(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t k,
        const Reads&,
        const Markers& markers
        );

}

#endif
