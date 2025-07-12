#ifndef SHASTA_READ_LOADER_HPP
#define SHASTA_READ_LOADER_HPP

// Shasta
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "Reads.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"

namespace shasta {
    class ReadLoader;

    extern template class MultithreadedObject<ReadLoader>;
}



// Class used to load reads from a fasta file.
class shasta::ReadLoader :
    public MultithreadedObject<ReadLoader>{
public:

    // The constructor does all the work.
    ReadLoader(
        const string& fileName,
        uint64_t minReadLength,
        size_t threadCount,
        const string& dataNamePrefix,
        size_t pageSize,
        Reads& reads);

    ~ReadLoader();

private:

    // The name of the file we are processing.
    const string& fileName;

    // The minimum read length. Shorter reads are not stored.
    const uint64_t minReadLength;

    // The number of threads to be used for processing.
    // Reading is done single-threaded as there is usually no benefit
    // frm multithreaded reading.
    size_t threadCount;

    // Information that we can use to create temporary
    // memory mapped binary data structures.
    const string& dataNamePrefix;
    const size_t pageSize;

    // The data structure that the reads will be added to.
    Reads& reads;

    // Create the name to be used for a MemoryMapped object.
    string dataName(
        const string& dataName) const;
    string threadDataName(
        size_t threadId,
        const string& dataName) const;

    // Read an entire file into a buffer,
    // using threadCountForReading threads.
    int64_t fileSize;
    MemoryMapped::Vector<char> buffer;
    void allocateBuffer();
    bool readFile(bool useODirect);
    void allocateBufferAndReadFile();

    // Vectors where each thread stores the reads it found.
    // Indexed by threadId.
    vector< unique_ptr<MemoryMapped::VectorOfVectors<char, uint64_t> > > threadReadNames;
    vector< unique_ptr<LongBaseSequences> > threadReads;
    void allocatePerThreadDataStructures();
    void allocatePerThreadDataStructures(size_t threadId);

    // Store the reads computed by each thread and free
    // the per-thread data structures.
    void storeReads();

    // Functions used for fasta files.
    void processFastaFile();
    void processFastaFileThreadFunction(size_t threadId);
    // Function that returns true if a read begins
    // at this position in Fasta format.
    bool fastaReadBeginsHere(uint64_t offset) const;

    // Functions used for fastq files.
    void processFastqFile();
    void processFastqFileThreadFunction(size_t threadId);

    // Find all line ends in the file.
    void findLineEnds();
    void findLineEndsThreadFunction(size_t threadId);
    vector< vector<uint64_t> > threadLineEnds;
    vector<uint64_t> lineEnds;

};



#endif
