#include "MarkerInfo.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta2;



// Construct the reverse complement of a MarkerInfo.
MarkerInfo MarkerInfo::reverseComplement(
    const Reads& reads,
    const Markers& markers) const
{
    MarkerInfo markerInfoRc = *this;
    const uint32_t orientedReadMarkerCount = uint32_t(markers[orientedReadId.getValue()].size());
    const uint32_t orientedReadBaseCount = uint32_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

    markerInfoRc.orientedReadId.flipStrand();
    markerInfoRc.ordinal = orientedReadMarkerCount - 1 - markerInfoRc.ordinal;
    markerInfoRc.position = orientedReadBaseCount - markerInfoRc.position;

    return markerInfoRc;
}
