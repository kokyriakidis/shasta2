#include "MarkerInfo.hpp"
#include "Reads.hpp"
using namespace shasta2;



// Construct the reverse complement of a MarkerInfo.
MarkerInfo MarkerInfo::reverseComplement(const Reads& reads) const
{
    MarkerInfo markerInfoRc = *this;
    const uint32_t orientedReadBaseCount = uint32_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

    markerInfoRc.orientedReadId.flipStrand();
    markerInfoRc.position = orientedReadBaseCount - markerInfoRc.position;

    return markerInfoRc;
}



// Get the Kmer corresponding to this MarkerInfo.
Kmer MarkerInfo::getKmer(uint64_t k, const Reads& reads) const
{
    // The MarkerInfo contains the mid position of the marker.
    return reads.getKmer(k, orientedReadId, position - uint32_t(k/2));
}
