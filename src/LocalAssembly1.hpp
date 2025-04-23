#pragma once

#include "Anchor.hpp"

namespace shasta {
    class LocalAssembly1;
}


class shasta::LocalAssembly1 {
public:
    LocalAssembly1(
        const Anchors&,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        bool computeAlignment,
        ostream& html);

private:
    const Anchors& anchors;
    ostream& html;

    // Information for each oriented read used in this local assembly.
    // We use all common oriented reads with positive ordinal offset
    // between anchorIdA and anchorIdB.
    // Stored sorted by OrientedReadId.
    class OrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t positionInJourneyA;
        uint32_t positionInJourneyB;
        uint32_t ordinalA;
        uint32_t ordinalB;

        // Base positions of the mid point of anchorIdA and anchorIdB.
        uint32_t basePositionA;
        uint32_t basePositionB;
        uint32_t sequenceLength() const
        {
            return basePositionB - basePositionA;
        }

        // The sequence of this oriented read in the region covered by this LocalAssembly1.
        vector<Base> sequence;

        // Sort by decrease sequence length.
        bool operator<(const OrientedRead& that) const
        {
            return sequenceLength() > that.sequenceLength();
        }
    };
    vector<OrientedRead> orientedReads;
    void gatherOrientedReads(
        AnchorId anchorIdA,
        AnchorId anchorIdB);
    void writeOrientedReads() const;

    // Consensus is always computed.
    // Alignment and alignedConsensus are only computed if computeAlignment is set to true.
    void runAbpoa(bool computeAlignment);
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    void writeConsensus() const;
    void writeAlignment() const;

};
