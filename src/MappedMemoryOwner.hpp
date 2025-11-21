#pragma once

#include "cstdint.hpp"
#include "string.hpp"

namespace shasta2 {
    class MappedMemoryOwner;
}



class shasta2::MappedMemoryOwner {
public:

    string largeDataFileNamePrefix;
    uint64_t largeDataPageSize;

    // Function to construct names for binary objects.
    // The output can be passed to createNew or accessExisting
    // member functions of MemoryMapped obkects.
    string largeDataName(const string& name) const
    {
        if(largeDataFileNamePrefix.empty()) {
            return "";  // Anonymous;
        } else {
            return largeDataFileNamePrefix + name;
        }
    }

    MappedMemoryOwner() {}
    MappedMemoryOwner(const MappedMemoryOwner&) = default;

    MappedMemoryOwner(
        const string& largeDataFileNamePrefix,
        uint64_t largeDataPageSize) :
        largeDataFileNamePrefix(largeDataFileNamePrefix),
        largeDataPageSize(largeDataPageSize)
    {}

    template<class T> void createNew(T& t, const string& name)
    {
        t.createNew(largeDataName(name), largeDataPageSize);
    }
    template<class T> void accessExistingReadOnly(T& t, const string& name)
    {
        t.accessExistingReadOnly(largeDataName(name));
    }
};

