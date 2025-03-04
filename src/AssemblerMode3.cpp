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



void Assembler::accessMode3Assembler()
{
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(MappedMemoryOwner(*this), getReads(), assemblerInfo->k, markers());
    mode3Assembler = make_shared<Mode3Assembler>(*this,
        assemblerInfo->k, getReads(), markers(),
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
    mode3Assembler->exploreJourney(request, html);
}



void Assembler::exploreReadFollowing(const vector<string>& request, ostream& html)
{
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
