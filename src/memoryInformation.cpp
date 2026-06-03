#include "memoryInformation.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta2;

#include "fstream.hpp"
#include "iostream.hpp"
#include <sstream>
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



// Get current and peak virtual memory utilization of the current process, in bytes.
VirtualMemoryInfo shasta2::getVirtualMemoryInfo()
{
    VirtualMemoryInfo virtualMemoryInfo;

    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        uint64_t done = 0;
        while((done < 2) and std::getline(procStats, line)) {
            std::istringstream is(line);
            string name;
            string value;
            is >> name >> value;
            if(name == "VmSize:") {
                virtualMemoryInfo.current = 1024UL * std::stol(value);
                ++done;
            }
            if(name == "VmPeak:") {
                virtualMemoryInfo.peak = 1024UL * std::stol(value);
                ++done;
            }
        }
    }

    return virtualMemoryInfo;
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



void shasta2::writeMemoryStatistics(const string& name)
{
    performanceLog <<
        timestamp << "At " << name << ":\n" <<
        "    Peak virtual memory usage " << getPeakMemoryUsage() << "\n" <<
        flush;

}

