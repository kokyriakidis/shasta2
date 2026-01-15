#include "memoryInformation.hpp"
using namespace shasta2;

#include "fstream.hpp"
#include "string.hpp"



uint64_t shasta2::getPeakMemoryUsage() {
    uint64_t peakMemoryUsage = 0ULL;

    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        while (std::getline(procStats, line)) {
            if (string::npos == line.find("VmPeak")) {
                continue;
            }
            size_t pos = line.find(":");
            while (pos < line.size() && !isdigit(line[pos])) {
                pos++;
            }
            char* end;
            peakMemoryUsage = std::strtoull(line.c_str() + pos, &end, 10);
            // Convert from kB to bytes.
            peakMemoryUsage *= 1024;
            break;
        }
    }

    return peakMemoryUsage;
}



// Get total physical memory available, in bytes.
uint64_t shasta2::getTotalPhysicalMemory()
{
    ifstream meminfo("/proc/meminfo");
    string s;
    uint64_t memoryKb;
    meminfo >> s >> memoryKb;
    return 1024 * memoryKb;
}
