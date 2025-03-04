#include "Markers.hpp"
#include "Kmer.hpp"
#include "KmerChecker.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Markers>;



Markers::Markers(
    const MappedMemoryOwner& mappedMemoryOwner,
    size_t k,
    const shared_ptr<const KmerChecker> kmerChecker,
    const shared_ptr<const Reads> reads,
    size_t threadCount) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<Markers>(*this),
    k(k),
    kmerChecker(kmerChecker),
    reads(reads)
{
    MemoryMapped::VectorOfVectors<Marker, uint64_t>::createNew(largeDataName("Markers"), largeDataPageSize);

    // Initial message.
    const uint64_t readCount = reads->readCount();
    performanceLog << timestamp << "Finding markers in " << readCount << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }


    const size_t batchSize = 100;
    beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    pass = 1;
    runThreads(&Markers::threadFunction, threadCount);
    beginPass2();
    endPass2(false);
    setupLoadBalancing(readCount, batchSize);
    pass = 2;
    runThreads(&Markers::threadFunction, threadCount);

    unreserve();



    // Final message.
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Finding markers completed in " << tTotal << " s." << endl;
    cout << "Created " << totalSize() << " markers." << endl;

}



void Markers::threadFunction(uint64_t)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads of this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            const LongBaseSequenceView read = reads->getRead(readId);
            size_t markerCount = 0; // For this read.
            Marker* markerPointerStrand0 = 0;
            Marker* markerPointerStrand1 = 0;
            if(pass == 2) {
                markerPointerStrand0 = this->begin(OrientedReadId(readId, 0).getValue());
                markerPointerStrand1 = this->end(OrientedReadId(readId, 1).getValue()) - 1ULL;
            }

            if(read.baseCount >= k) {   // Avoid pathological case.

                // Loop over k-mers of this read.
                Kmer kmer;
                for(size_t position=0; position<k; position++) {
                    kmer.set(position, read[position]);
                }
                for(uint32_t position=0; /*The check is done later */; position++) {
                    const KmerId kmerId = KmerId(kmer.id(k));
                    if(kmerChecker->isMarker(kmerId)) {
                        // This k-mer is a marker.

                        if(pass == 1) {
                            ++markerCount;
                        } else {
                            // Strand 0.
                            markerPointerStrand0->position = position;
                            ++markerPointerStrand0;

                            // Strand 1.
                            markerPointerStrand1->position = uint32_t(read.baseCount - k - position);
                            --markerPointerStrand1;

                        }
                    }

                    if(position+k == read.baseCount) {
                        break;
                    }

                    // Update the k-mer.
                    kmer.shiftLeft();
                    kmer.set(k-1, read[position+k]);
                }
            }

            if(pass == 1) {
                incrementCount(OrientedReadId(readId, 0).getValue(), markerCount);
                incrementCount(OrientedReadId(readId, 1).getValue(), markerCount);
            } else {
                SHASTA_ASSERT(markerPointerStrand0 ==
                    this->end(OrientedReadId(readId, 0).getValue()));
                SHASTA_ASSERT(markerPointerStrand1 ==
                    this->begin(OrientedReadId(readId, 1).getValue()) - 1ULL);
            }
        }
    }

}


Markers::Markers(const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<Markers>(*this)
{
    MemoryMapped::VectorOfVectors<Marker, uint64_t>::accessExistingReadOnly(largeDataName("Markers"));
}
