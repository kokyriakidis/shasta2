#pragma once

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"

// Standard library.
#include <iosfwd.hpp>
#include <vector.hpp>



namespace shasta {
    class LocalAssembly3;

    class Anchors;
    class AnchorPair;
}



class shasta::LocalAssembly3 {
public:

    // This assembles between anchorIdA and anchorIdB
    // of the given AnchorPair. It uses all the OrientedReadIds
    // stored in the AnchorPair, and which all appear in both
    // anchorIdA and anchorIdB. In addition, it uses OrientedReadIds
    // stored in additionalOrientedReadIds that:
    // - Are not also in the AnchorPair.
    // - Appear in at least one of anchorIdA and anchorIdB.
    // The additionalOrientedReadIds must be sorted.
    LocalAssembly3(
        const Anchors&,
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

public:

    // The two anchors that bound the region assembled by this local assembly.
    AnchorId anchorIdA;
    AnchorId anchorIdB;



    // Information for each oriented read used in this local assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        bool isOnLeftAnchor;
        bool isOnRightAnchor;
        bool isOnBothAnchors() const
        {
            return isOnLeftAnchor and isOnRightAnchor;
        }

        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnLeftAnchor);
            SHASTA_ASSERT(isOnRightAnchor);
            return rightOrdinal - leftOrdinal;
        }

        // The left/right positions of this marker.
        // These are positions in the oriented read sequence
        // of the leftmost base of the marker.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            SHASTA_ASSERT(isOnLeftAnchor);
            SHASTA_ASSERT(isOnRightAnchor);
            return rightPosition - leftPosition;
        }
    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads(
        const Anchors&,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // Html output.
    void writeInput(
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeOrientedReads(
        const Anchors&,
        ostream& html) const;
};
