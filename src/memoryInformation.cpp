// Shasta2.
#include "memoryInformation.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

// Linux.
#include <malloc.h>
#include <stdio.h>

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"
#include <sstream>



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



void MallocInfo::write(ostream& s) const
{
    s << "Heap size " << heapSize << "\n";
    s << "Heap used " << heapUsed() << "\n";
    s << "Heap free " << heapFree << "\n";
    s << "Maximum heap size " << maxHeapSize << "\n";
    s << "Mmap size " << mmapSize << "\n";
    s << flush;
}



MallocInfo shasta2::getMallocInfo()
{

    // Create a FILE* for use with malloc_info.
    char* buffer = 0;
    uint64_t bufferSize = 0;
    FILE* file = ::open_memstream(&buffer, &bufferSize);
    if(not file) {
        throw runtime_error("Error in open_memstream.");
    }

    // Call malloc_info.
    const int returnCode = ::malloc_info(0, file);
    if(returnCode == -1) {
        throw runtime_error("Error in malloc_info.");
    }
    ::fclose(file);

    // Create a string to use to create the property tree.
    // We could use a string_view but c++23
    // would not support that constructor.
    const string s(buffer, bufferSize);
    std::istringstream is(s);
    // cout << s << flush;

    // Create the property tree.
    boost::property_tree::ptree propertyTree;
    boost::property_tree::read_xml(is, propertyTree);
    ::free(buffer);

    // boost::property_tree::write_info(cout, propertyTree);

    // Extract the information we want.
    MallocInfo mallocInfo;
    for(const auto& a: propertyTree) {
        for(const auto& b: a.second) {

            if(b.first == "total") {
                for(const auto& c: b.second) {
                    const string type = c.second.get<string>("type");
                    // cout << "total type " << c.second.get<string>("type") << endl;
                    // cout << "total count " << c.second.get<uint64_t>("count") << endl;
                    // cout << "total size " << c.second.get<uint64_t>("size") << endl;
                    if((type == "fast") or (type == "rest")) {
                        mallocInfo.heapFree += c.second.get<uint64_t>("size");
                    }
                    if(type == "mmap") {
                        mallocInfo.mmapSize += c.second.get<uint64_t>("size");
                    }
                }
            }

            if(b.first == "system") {
                for(const auto& c: b.second) {
                    const string type = c.second.get<string>("type");
                    // cout << "system type " << c.second.get<string>("type") << endl;
                    // cout << "system size " << c.second.get<uint64_t>("size") << endl;
                    if(type == "current") {
                        mallocInfo.heapSize = c.second.get<uint64_t>("size");
                    }
                    if(type == "max") {
                        mallocInfo.maxHeapSize = c.second.get<uint64_t>("size");
                    }
                }
            }
        }
    }

    return mallocInfo;

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
    const VirtualMemoryInfo virtualMemoryInfo = getVirtualMemoryInfo();

    performanceLog <<
        timestamp <<
        "At " << name << ":\n" <<
        "Virtual memory (current) " << virtualMemoryInfo.current << "\n" <<
        "Virtual memory (peak) " << virtualMemoryInfo.peak << "\n";
    getMallocInfo().write(performanceLog);

}

