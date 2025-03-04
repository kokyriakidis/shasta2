// Shasta.
#include "Mode3Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Standard library.
#include "iostream.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Mode3Assembler>;



// This constructor runs the assembly.
Mode3Assembler::Mode3Assembler(
    const MappedMemoryOwner& mappedMemoryOwner,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<Marker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer,
    uint64_t /* threadCount */,
    const Mode3AssemblyOptions& options,
    bool debug) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads),
    markers(markers),
    debug(debug),
    anchorsPointer(anchorsPointer),
    options(options)
{
    SHASTA_ASSERT(anchorsPointer);
    anchors().writeCoverageHistogram();
    anchors().writeAnchorGapsByRead();
}



// This constructor just accesses binary data.
Mode3Assembler::Mode3Assembler(
    const MappedMemoryOwner& mappedMemoryOwner,
    uint64_t k,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<Marker, uint64_t>& markers,
    shared_ptr<mode3::Anchors> anchorsPointer,
    const Mode3AssemblyOptions& options) :
    MultithreadedObject<Mode3Assembler>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    k(k),
    reads(reads),
    markers(markers),
    anchorsPointer(anchorsPointer),
    options(options)
{
    SHASTA_ASSERT(anchorsPointer);
}



