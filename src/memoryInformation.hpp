#pragma once

#include <cstdint>

namespace shasta2 {

    // Get peak virtual memory utilization of the current process, in bytes.
    uint64_t getPeakMemoryUsage();

    // Get total physical memory available, in bytes.
    uint64_t getTotalPhysicalMemory();

}
