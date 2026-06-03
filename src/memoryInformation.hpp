#pragma once

#include <cstdint>
#include "string.hpp"

namespace shasta2 {

    // Get peak virtual memory utilization of the current process, in bytes.
    uint64_t getPeakMemoryUsage();

    // Get current and peak virtual memory utilization of the current process, in bytes.
    class VirtualMemoryInfo {
    public:
        uint64_t current = 0;
        uint64_t peak = 0;
    };
    VirtualMemoryInfo getVirtualMemoryInfo();

    // Get total physical memory available, in bytes.
    uint64_t getTotalPhysicalMemory();


    void writeMemoryStatistics(const string&);
}
