#pragma once

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"



namespace shasta2 {
    class LocalAssembly4;

    class Anchors;
    class AnchorPair;
}




class shasta2::LocalAssembly4 {
public:

    // This assembles between anchorIdA and anchorIdB
    // of the given AnchorPair. It uses all the OrientedReadIds
    // stored in the AnchorPair, which all appear in both
    // anchorIdA and anchorIdB. In addition, it uses OrientedReadIds
    // stored in additionalOrientedReadIds that:
    // - Are not also in the AnchorPair.
    // - Appear in at least one of anchorIdA and anchorIdB.
    // The additionalOrientedReadIds must be sorted.
    // The additionalOrientedReadIds are allowed to contain OrientedReadIds
    // that are also in the AnchorPair.
    LocalAssembly4(
        const Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

    // Assembled sequence and its coverage.
    vector<Base> sequence;
    vector<uint64_t> coverage;

private:

    const Anchors& anchors;
    ostream& html;

    // The two anchors of the AnchorPair used for this assembly.
    AnchorId leftAnchorId;
    AnchorId rightAnchorId;

    // The union of the OrientedReadIds in the AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> allOrientedReadIds;
    void gatherAllOrientedReads(
        const AnchorPair&,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // Information on OrientedReadIds that are in
    // the usableOrientedReadIds and also appear on both
    // the left and right anchors.
    class CommonOrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The ordinals of this OrientedReadId in the left/right anchor.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right anchor.
        // These are mid positions of the marker in the oriented read sequence.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            return rightPosition - leftPosition;
        }

    };
    vector<CommonOrientedReadInfo> commonOrientedReadInfos;
    void gatherCommonOrientedReads();

    // Get the sequence a CommonOrientedReadInfo contributes to the assembly.
    void getSequence(const CommonOrientedReadInfo&, vector<Base>&) const;

    // This does the multiple sequence alignment and stores sequence and coverage.
    void assemble();

    // Html output.
    void writeInput(
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeAllOrientedReadIds() const;
    void writeCommonOrientedReads() const;
    void writeAssembledSequence() const;
};
