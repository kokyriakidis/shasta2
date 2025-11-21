// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

#include "performanceLog.hpp"

namespace shasta2 {
    ofstream performanceLog;
}



void shasta2::openPerformanceLog(const string& fileName)
{
    performanceLog.open(fileName);
}
