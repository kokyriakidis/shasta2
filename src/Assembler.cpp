#include "Assembler.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "Options.hpp"
#include "AssemblyGraph.hpp"
#include "Journeys.hpp"
#include "KmerCheckerFactory.hpp"
#include "MurmurHash2.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "ReadLengthDistribution.hpp"
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
    const Options& options,
    vector<string> inputFileNames)
{
    cout << "Number of threads: " << options.threadCount << endl;

    addReads(
        inputFileNames,
        options.minReadLength,
        options.threadCount);
    createReadLengthDistribution();

    createKmerChecker(options.k, options.markerDensity, options.threadCount);
    createMarkers(options.threadCount);
    createMarkerKmers(options.maxMarkerErrorRate, options.threadCount);
    createAnchors(
        options.minAnchorCoverage,
        options.maxAnchorCoverage,
        options.maxAnchorHomopolymerLength,
        options.threadCount);

    createJourneys(options.threadCount);
    journeys().writeAnchorGapsByRead(reads(), markers(), anchors());

    createAnchorGraph(options);

    createAssemblyGraph(options);
}





void Assembler::createKmerChecker(
    uint64_t k,
    double markerDensity,
    uint64_t /* threadCount */) // To permit future multithreaded KmerChecker constructors.
{

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


AnchorId Assembler::readFollowing(
    AnchorId anchorId,
    uint64_t direction,
    uint64_t minCommonCount,
    double aDrift,
    double bDrift
    ) const
{
    return anchors().readFollowing(journeys(), anchorId, direction, minCommonCount, aDrift, bDrift);
}



void Assembler::createAnchorGraph(const Options& options)
{

    anchorGraphPointer = make_shared<AnchorGraph>(
        anchors(), journeys(),
        options.minAnchorGraphEdgeCoverage,
        options.minAnchorGraphContinueReadFollowingCount,
        options.aDrift,
        options.bDrift,
        options.threadCount);
    anchorGraphPointer->save("AnchorGraph");

}



void Assembler::createSimpleAnchorGraph()
{
    simpleAnchorGraphPointer = make_shared<AnchorGraph>(anchors(), journeys());
    simpleAnchorGraphPointer->save("SimpleAnchorGraph");
}



void Assembler::createAssemblyGraph(const Options& options)
{
    AssemblyGraph assemblyGraph(
        anchors(),
        journeys(),
        *anchorGraphPointer,
        options);
    assemblyGraph.run();
}



void Assembler::accessAnchorGraph()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;
    anchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner, "AnchorGraph");
}



void Assembler::accessSimpleAnchorGraph()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;
    simpleAnchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner, "SimpleAnchorGraph");
}
