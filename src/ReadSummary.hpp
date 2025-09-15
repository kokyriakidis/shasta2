#pragma once

namespace shasta {
    class ReadSummary;
}



// Summary information about a read in the assembly.
class shasta::ReadSummary {
public:

    bool isUsedForAssembly = true;

    // The marker error rate computed using all reads.
    // Reads that have this value larger than Options::maxMarkerErrorRate
    // get their isUsedForAssembly set to false and are not used for
    // assembly.
    double initialMarkerErrorRate = 0.;

    // The marker error rate computed using the reads with the
    // isUsedForassembly flag set.
    double markerErrorRate = 0.;

    // The number of bases before the first anchor.
    uint32_t initialAnchorGap;

    // The width, in bases, of the largest gap between adjacent anchors.
    uint32_t middleAnchorGap;

    // The number of bases after the last anchor.
    uint32_t finalAnchorGap;

};
