// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

#include "performanceLog.hpp"
#include "memoryInformation.hpp"
#include "timestamp.hpp"

#include "iostream.hpp"

namespace shasta2 {
    ofstream performanceLog;
}



void shasta2::openPerformanceLog(const string& fileName)
{
    performanceLog.open(fileName);
}



void shasta2::writePerformanceStatistics(const string& name)
{
    performanceLog <<
        timestamp << "At " << name << ":\n" <<
        "    Peak virtual memory usage " << getPeakMemoryUsage() << "\n" <<
        flush;

}
