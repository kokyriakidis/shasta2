// Shasta
#include "Assembler.hpp"
#include "mode3-LocalAssembly.hpp"
#include "Mode3Assembler.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>



void Assembler::mode3Assembly(
    uint64_t threadCount,
    shared_ptr<mode3::Anchors> anchorsPointer,
    const Mode3AssemblyOptions& options,
    bool debug
    )
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    mode3Assembler = make_shared<Mode3Assembler>(
        mappedMemoryOwner,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, threadCount, options, debug);
}



// Same, but use existing Anchors. Python callable.
void Assembler::mode3Reassembly(
    uint64_t threadCount,
    const Mode3AssemblyOptions& options,
    bool debug
    )
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    // Create the Anchors from binary data.
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(mappedMemoryOwner, getReads(), assemblerInfo->k, markers);

    // Run the Mode 3 assembly.
    mode3Assembler = make_shared<Mode3Assembler>(
        mappedMemoryOwner,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, threadCount, options, debug);
}

void Assembler::accessMode3Assembler()
{
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(MappedMemoryOwner(*this), getReads(), assemblerInfo->k, markers);
    mode3Assembler = make_shared<Mode3Assembler>(*this,
        assemblerInfo->k, getReads(), markers,
        anchorsPointer, httpServerData.assemblerOptions->assemblyOptions.mode3Options);
}



void Assembler::exploreAnchor(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreAnchor(request, html);
}



void Assembler::exploreAnchorPair(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreAnchorPair(request, html);
}



void Assembler::exploreJourney(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);
    mode3Assembler->exploreJourney(request, html);
}



void Assembler::exploreReadFollowing(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->assemblyMode == 3);
    mode3Assembler->exploreReadFollowing(request, html);
}



void Assembler::exploreLocalAssembly(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreLocalAssembly(request, html);
}



void Assembler::exploreLocalAnchorGraph(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreLocalAnchorGraph(request, html);
}



// Alignment-free version of mode 3 assembly.
void Assembler::alignmentFreeAssembly(
    const Mode3AssemblyOptions& mode3Options,
    uint64_t threadCount)
{
    createMarkerKmers(threadCount);

    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(
            MappedMemoryOwner(*this),
            getReads(),
            assemblerInfo->k,
            markers,
            markerKmers,
            mode3Options.minAnchorCoverage,
            mode3Options.maxAnchorCoverage,
            threadCount);

    // Compute oriented read journeys.
    anchorsPointer->computeJourneys(threadCount);

    // Run Mode 3 assembly.
    mode3Assembly(threadCount, anchorsPointer, mode3Options, false);
}
