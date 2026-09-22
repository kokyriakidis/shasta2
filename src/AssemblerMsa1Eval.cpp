// Support for the msa1 hard-region evaluation harness
// (scripts/FindMsa1HardRegions.py, scripts/EvaluateMsa1AgainstTruth.py).
// See the declarations in Assembler.hpp for what these do and why.

// Shasta.
#include "Assembler.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
#include "LocalAssembly7.hpp"
#include "msa1.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
using namespace shasta2;



string Assembler::getOrientedReadSequenceString(const string& orientedReadIdString) const
{
    const OrientedReadId orientedReadId(orientedReadIdString);
    const vector<Base> sequence = anchors().reads.getOrientedReadSequence(orientedReadId);
    return toString(sequence);
}



bool Assembler::anchorContainsOrientedRead(
    AnchorId anchorId,
    const string& orientedReadIdString) const
{
    const OrientedReadId orientedReadId(orientedReadIdString);
    return anchors().anchorContains(anchorId, orientedReadId);
}



uint32_t Assembler::getAnchorPositionInOrientedRead(
    AnchorId anchorId,
    const string& orientedReadIdString) const
{
    const OrientedReadId orientedReadId(orientedReadIdString);
    return anchors().getPosition(anchorId, orientedReadId);
}



std::tuple<bool, string, bool, string> Assembler::runLocalAssemblyWithAndWithoutMsa1Repair(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    const vector<string>& orientedReadIdStrings) const
{
    vector<OrientedReadId> orientedReadIds;
    orientedReadIds.reserve(orientedReadIdStrings.size());
    for(const string& s: orientedReadIdStrings) {
        orientedReadIds.push_back(OrientedReadId(s));
    }

    ostream html(0);

    // One run, with the repair on, gives both: LocalAssembly7::sequenceBeforeRepair
    // is the consensus before the repair, sequence is after. Running twice - once
    // with useMsa1 false, once true - would recompute the same alignment twice,
    // which is the expensive part of a run; the repair itself is a small fraction
    // of the cost.
    LocalAssembly7::Options options;
    options.useMsa1 = true;
    const LocalAssembly7 localAssembly(
        options, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    return std::make_tuple(
        localAssembly.success, toString(localAssembly.sequenceBeforeRepair),
        localAssembly.success, toString(localAssembly.sequence));
}
