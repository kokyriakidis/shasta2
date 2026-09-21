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

    LocalAssembly7::Options optionsNoRepair;
    const LocalAssembly7 localAssemblyNoRepair(
        optionsNoRepair, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    LocalAssembly7::Options optionsWithRepair;
    optionsWithRepair.useMsa1 = true;
    const LocalAssembly7 localAssemblyWithRepair(
        optionsWithRepair, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    return std::make_tuple(
        localAssemblyNoRepair.success, toString(localAssemblyNoRepair.sequence),
        localAssemblyWithRepair.success, toString(localAssemblyWithRepair.sequence));
}
