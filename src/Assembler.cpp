#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "KmerCheckerFactory.hpp"
#include "mode3-Anchor.hpp"
#include "MurmurHash2.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;

#include "MultithreadedObject.tpp"
template class MultithreadedObject<Assembler>;


// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefixArgument,
    bool createNew,
    size_t largeDataPageSizeArgument) :
    MultithreadedObject(*this)
{
    largeDataFileNamePrefix = largeDataFileNamePrefixArgument;

    if(createNew) {

        // Create a new assembly.
        assemblerInfo.createNew(largeDataName("Info"), largeDataPageSizeArgument);
        assemblerInfo->largeDataPageSize = largeDataPageSizeArgument;
        largeDataPageSize = largeDataPageSizeArgument;

        readsPointer = make_shared<Reads>();
        readsPointer->createNew(
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadIdsSortedByName"),
            largeDataPageSize
        );
        // cout << "Created a new assembly with page size " << largeDataPageSize << endl;

    } else {

        // Access an existing assembly.
        assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
        largeDataPageSize = assemblerInfo->largeDataPageSize;

        readsPointer = make_shared<Reads>();
        readsPointer->access(
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadIdsSortedByName")
        );

    }
    SHASTA_ASSERT(largeDataPageSize == assemblerInfo->largeDataPageSize);

    // In both cases, assemblerInfo, reads, readNames, readRepeatCounts are all open for write.

    fillServerFunctionTable();
}




// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void Assembler::assemble(
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFileNames)
{
    // Adjust the number of threads, if necessary.
    uint64_t threadCount = assemblerOptions.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Number of threads: " << threadCount << endl;

    // Add reads from the specified input files.
    addReads(
        inputFileNames,
        assemblerOptions.minReadLength,
        threadCount);

    // Initialize the KmerChecker, which has the information needed
    // to decide if a k-mer is a marker.
    createKmerChecker(assemblerOptions.k, assemblerOptions.markerDensity, threadCount);

    // Create the markers.
    createMarkers(threadCount);

    // Create MarkerKmers.
    createMarkerKmers(threadCount);

    // Create Anchors.
    createAnchors(
        assemblerOptions.minAnchorCoverage,
        assemblerOptions.maxAnchorCoverage,
        threadCount);

}





void Assembler::createKmerChecker(
    uint64_t k,
    double markerDensity,
    uint64_t threadCount)
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    assemblerInfo->k = k;

    kmerChecker = KmerCheckerFactory::createNew(
        k,
        markerDensity,
        *this);
}



void Assembler::accessKmerChecker()
{
    kmerChecker = KmerCheckerFactory::createFromBinaryData(*this);
}



void Assembler::createAnchors(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    anchorsPointer = make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        reads(),
        assemblerInfo->k,
        markers(),
        markerKmers,
        minAnchorCoverage,
        maxAnchorCoverage,
        threadCount);

    anchorsPointer->computeJourneys(threadCount);
}



void Assembler::accessAnchors()
{
     anchorsPointer = make_shared<mode3::Anchors>(
         MappedMemoryOwner(*this), reads(), assemblerInfo->k, markers());
}

