#pragma once

/*******************************************************************************

Among all 4^k k-mers of length k, we choose a subset that we call "markers".
The markers are selected at the beginning of an assembly
and never changed, and selected in such a way that,
if (and only if) a k-mer is a marker, its reverse complement
is also a marker.

*******************************************************************************/

// Shasta.
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"
#include "Uint.hpp"

// Standard library.
#include "memory.hpp"

namespace shasta2 {

    class Marker;
    class Markers;

    class KmerChecker;
    class Reads;
    class OrientedReadId;
}



// Markers in shared memory are stored using class Marker.
class shasta2::Marker {
public:

    // The position of this marker in the oriented read.
    // This limits the length of a read to 2^24=16Mib bases.
    Uint24 position;

    bool operator<(const Marker& that) const
    {
        return position < that.position;
    }
};



// The markers on all oriented reads. Indexed by OrientedReadId::getValue().
class shasta2::Markers:
    public MappedMemoryOwner,
    public MultithreadedObject<Markers>,
    public MemoryMapped::VectorOfVectors<Marker, uint64_t> {
public:

    Markers(
        const MappedMemoryOwner&,
        uint64_t k,
        const shared_ptr<const Reads>,
        const shared_ptr<const KmerChecker>,
        size_t threadCount);

    Markers(
        const MappedMemoryOwner&,
        uint64_t k,
        const shared_ptr<const Reads>);

    // Access functions for markers Kmers and KmerIds.
    Kmer getKmer(
        OrientedReadId,
        uint32_t ordinal) const;

    // These are filled in by all constructors.
    size_t k;
    const Reads& reads;

private:

    // The remaining arguments are only filled in by the first constructor,
    // which constructs the Markers from scratch.
    const shared_ptr<const KmerChecker> kmerChecker;
    size_t threadCount;

    void threadFunction(size_t threadId);

    // In pass 1, we count the number of markers for each
    // read and call reads->incrementCountMultithreaded.
    // In pass 2, we store the markers.
    size_t pass;

    // Low level, private access functions.
    Kmer getKmerStrand0(
        ReadId,
        uint32_t ordinal) const;

    Kmer getKmerStrand1(
        ReadId,
        uint32_t ordinal) const;
};
