#pragma once

#include <cstdint>
#include "iosfwd.hpp"
#include "string.hpp"



namespace shasta2 {

    // Get peak virtual memory utilization of the current process, in bytes.
    uint64_t getPeakMemoryUsage();

    // Get total physical memory available, in bytes.
    uint64_t getTotalPhysicalMemory();



    // Get current and peak virtual memory utilization of the current process, in bytes.
    class VirtualMemoryInfo {
    public:
        uint64_t current = 0;
        uint64_t peak = 0;
    };
    VirtualMemoryInfo getVirtualMemoryInfo();



    // Get statistics about memory allocated by malloc.
    // mallinfo and mallinfo2 cannot be used as they
    // return inaccurate/incomplete information.
    class MallocInfo {
    public:
        uint64_t heapSize = 0;
        uint64_t maxHeapSize = 0;
        uint64_t mmapSize = 0;
        uint64_t heapFree = 0;
        uint64_t heapUsed() const
        {
            return heapSize - heapFree;
        }
        void write(ostream&) const;
    };
    MallocInfo getMallocInfo();



    void writeMemoryStatistics(const string&);
}
