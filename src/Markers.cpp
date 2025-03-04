#include "Markers.hpp"
#include "MarkerFinder.hpp"
using namespace shasta;



Markers::Markers(
    const MappedMemoryOwner& mappedMemoryOwner,
    size_t k,
    const KmerChecker& kmerChecker,
    const Reads& reads,
    size_t threadCount) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    MemoryMapped::VectorOfVectors<Marker, uint64_t>::createNew(largeDataName("Markers"), largeDataPageSize);

    MarkerFinder markerFinder(
        k,
        kmerChecker,
        reads,
        *this,
        threadCount);

}



Markers::Markers(const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    MemoryMapped::VectorOfVectors<Marker, uint64_t>::accessExistingReadOnly(largeDataName("Markers"));
}
