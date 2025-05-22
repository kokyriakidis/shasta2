#include "Assembler.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "AssemblerOptions.hpp"
#include "AssemblyGraph.hpp"
#include "AssemblyGraph3.hpp"
#include "Journeys.hpp"
#include "KmerCheckerFactory.hpp"
#include "MurmurHash2.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "ReadLengthDistribution.hpp"
#include "TransitionGraph.hpp"
using namespace shasta;

#include "MultithreadedObject.tpp"
template class MultithreadedObject<Assembler>;


// Construct a new Assembler.
Assembler::Assembler(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize) :
    MultithreadedObject(*this),
    MappedMemoryOwner(largeDataFileNamePrefix, largeDataPageSize)
{


    assemblerInfo.createNew(largeDataName("Info"), largeDataPageSize);
    assemblerInfo->largeDataPageSize = largeDataPageSize;

    readsPointer = make_shared<Reads>();
    readsPointer->createNew(
        largeDataName("Reads"),
        largeDataName("ReadNames"),
        largeDataName("ReadIdsSortedByName"),
        largeDataPageSize
    );
}



// Construct an Assembler from binary data. This accesses the AssemblerInfo and the Reads.
Assembler::Assembler(const string& largeDataFileNamePrefix) :
    MultithreadedObject(*this),
    MappedMemoryOwner(largeDataFileNamePrefix, 0)
{

    assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
    largeDataPageSize = assemblerInfo->largeDataPageSize;

    readsPointer = make_shared<Reads>();
    readsPointer->access(
        largeDataName("Reads"),
        largeDataName("ReadNames"),
        largeDataName("ReadIdsSortedByName")
    );

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

    addReads(
        inputFileNames,
        assemblerOptions.minReadLength,
        threadCount);
    computeReadLengthDistribution();

    createKmerChecker(assemblerOptions.k, assemblerOptions.markerDensity, threadCount);
    createMarkers(threadCount);
    createMarkerKmers(threadCount);
    createAnchors(
        assemblerOptions.minAnchorCoverage,
        assemblerOptions.maxAnchorCoverage,
        assemblerOptions.maxAnchorHomopolymerLength,
        threadCount);

    createJourneys(threadCount);
    journeys().writeAnchorGapsByRead(reads(), markers(), anchors());

    createAssemblyGraph3(assemblerOptions, threadCount);
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
    uint64_t maxHomopolymerLength,
    uint64_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    anchorsPointer = make_shared<Anchors>(
        MappedMemoryOwner(*this),
        reads(),
        assemblerInfo->k,
        markers(),
        markerKmers,
        minAnchorCoverage,
        maxAnchorCoverage,
        maxHomopolymerLength,
        threadCount);
}



void Assembler::accessAnchors(bool writeAccess)
{
     anchorsPointer = make_shared<Anchors>(
         MappedMemoryOwner(*this), reads(), assemblerInfo->k, markers(), writeAccess);
}



void Assembler::createJourneys(uint64_t threadCount)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    journeysPointer = make_shared<Journeys>(
        2 * reads().readCount(),
        anchorsPointer,
        threadCount,
        mappedMemoryOwner);

}



void Assembler::accessJourneys()
{
    journeysPointer = make_shared<Journeys>(*this);
}



void Assembler::createAssemblyGraph(
    const AssemblerOptions& options,
    uint64_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "AnchorGraph creation begins." << endl;
    const AnchorGraph anchorGraph(anchors(), journeys(),
        options.minAnchorGraphEdgeCoverageNear);
    performanceLog << timestamp << "AnchorGraph creation ends." << endl;

    performanceLog << timestamp << "AssemblyGraph creation begins." << endl;
    AssemblyGraph assemblyGraph(options, anchors(), anchorGraph, threadCount);
    performanceLog << timestamp << "AssemblyGraph creation ends." << endl;

}



void Assembler::createAssemblyGraph3(
    const AssemblerOptions& options,
    uint64_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Generate an AnchorGraph in which all edges have consistent offsets.
    anchorGraphPointer = make_shared<AnchorGraph>(
        anchors(), journeys(),
        options.minAnchorGraphEdgeCoverageNear,
        options.minAnchorGraphEdgeCoverageFar,
        options.aDrift,
        options.bDrift);
    anchorGraphPointer->save();

    // Create the AssemblyGraph3.
    AssemblyGraph3 assemblyGraph3(
        anchors(),
        *anchorGraphPointer,
        options);
    anchorGraphPointer = 0;
    assemblyGraph3.run(threadCount);

}



void Assembler::accessAnchorGraph()
{
    anchorGraphPointer = make_shared<AnchorGraph>(*this);
}


void Assembler::computeReadLengthDistribution() const
{
    ReadLengthDistribution readLengthDistribution(reads(), MappedMemoryOwner(*this));
}
