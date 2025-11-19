#pragma once

#include "fstream.hpp"
#include "string.hpp"

// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

namespace shasta {
    extern ofstream performanceLog;
    void openPerformanceLog(const string& fileName);
}

