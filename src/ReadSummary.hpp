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
};
