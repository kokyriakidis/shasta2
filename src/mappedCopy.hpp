#pragma once

#include "string.hpp"

namespace shasta2 {

    // This can be used to copy a file to the huge page filesystem.
    // The regular cp command does not work (but it works to copy
    // the other way around, from the huge page filesystem).
    void mappedCopy(
        const string& inputPath,
        const string& outputPath);
}

