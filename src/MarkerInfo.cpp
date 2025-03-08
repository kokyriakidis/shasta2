#include "MarkerInfo.hpp"
#include "Markers.hpp"
using namespace shasta;



// Construct the reverse complement of a MarkerInfo.
MarkerInfo MarkerInfo::reverseComplement(const Markers& markers) const
{
    MarkerInfo markerInfoRc = *this;
    const uint32_t orientedReadMarkerCount = uint32_t(markers[orientedReadId.getValue()].size());

    markerInfoRc.orientedReadId.flipStrand();
    markerInfoRc.ordinal = orientedReadMarkerCount - 1 - markerInfoRc.ordinal;

    return markerInfoRc;
}
