#ifndef SHASTA_MAPPED_MEMORY_OWNER_HPP
#define SHASTA_MAPPED_MEMORY_OWNER_HPP

#include "cstdint.hpp"
#include "string.hpp"

namespace shasta {
    class MappedMemoryOwner;
}



class shasta::MappedMemoryOwner {
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

    template<class T> void createNew(T& t, const string& name)
    {
        t.createNew(largeDataName(name), largeDataPageSize);
    }
    template<class T> void accessExistingReadOnly(T& t, const string& name)
    {
        t.accessExistingReadOnly(largeDataName(name));
    }
};


#endif
