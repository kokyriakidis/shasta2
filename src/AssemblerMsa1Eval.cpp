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

// Standard library.
#include "iostream.hpp"



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



std::tuple<bool, string, bool, string> Assembler::runLocalAssemblyAdaptiveAndMsa1(
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

    LocalAssembly7::Options optionsAdaptive;
    optionsAdaptive.method = LocalAssembly7::Method::Adaptive;
    const LocalAssembly7 localAssemblyAdaptive(
        optionsAdaptive, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    LocalAssembly7::Options optionsMsa1;
    optionsMsa1.method = LocalAssembly7::Method::Msa1;
    const LocalAssembly7 localAssemblyMsa1(
        optionsMsa1, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    return std::make_tuple(
        localAssemblyAdaptive.success, toString(localAssemblyAdaptive.sequence),
        localAssemblyMsa1.success, toString(localAssemblyMsa1.sequence));
}



std::tuple<
    bool, string,
    vector< std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> >
    > Assembler::runLocalAssemblyMsa1WithDiagnostics(
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

    vector<Msa1ColumnDiagnostic> diagnostics;
    msa1ColumnDiagnostics = &diagnostics;

    LocalAssembly7::Options optionsMsa1;
    optionsMsa1.method = LocalAssembly7::Method::Msa1;
    const LocalAssembly7 localAssemblyMsa1(
        optionsMsa1, anchors(), anchorIdA, anchorIdB, html, orientedReadIds);

    msa1ColumnDiagnostics = nullptr;

    vector< std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> > rows;
    rows.reserve(diagnostics.size());
    for(const Msa1ColumnDiagnostic& d: diagnostics) {
        rows.push_back(std::make_tuple(
            d.totalWeight, d.maxObserved, d.medianLength, d.cumulativeAtMedian,
            d.weightAtMedian, d.weightAtMedianPlusOne, d.chosenLength));
    }

    return std::make_tuple(
        localAssemblyMsa1.success, toString(localAssemblyMsa1.sequence), rows);
}
