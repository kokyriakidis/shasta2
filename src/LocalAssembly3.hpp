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
    class Markers;
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
    AnchorId leftAnchorId;
    AnchorId rightAnchorId;



    // Information for each oriented read used in this local assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // Whether this OrientedReadId appear in the left and right anchors.
        bool isOnLeftAnchor;
        bool isOnRightAnchor;
        bool isOnBothAnchors() const
        {
            return isOnLeftAnchor and isOnRightAnchor;
        }

        // The ordinals of this OrientedReadId in the left/right anchor, if any.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnLeftAnchor);
            SHASTA_ASSERT(isOnRightAnchor);
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right anchor, if any.
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

        // All data members up to here are filled by gatherOrientedReads.

        // The first and last ordinal of the portion of this OrientedReadId
        // sequence that will be used in this local assembly.
        // These are filled by fillFirstLastOrdinalForAssembly.
        uint32_t firstOrdinalForAssembly;
        uint32_t lastOrdinalForAssembly;
        void fillFirstLastOrdinalForAssembly(const Markers&, uint32_t length);
        uint32_t firstPositionForAssembly(const Markers&) const;
        uint32_t lastPositionForAssembly(const Markers&) const;
    };

    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads(
        const Anchors&,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // The estimated base offset between the left and right anchor.
    // It is estimated using the oriented reads that appear in both anchors.
    uint32_t offset;
    void estimateOffset();

    // Use the estimated offset to fill in the firstOrdinal and lastOrdinal
    // of each oriented read. These define the portion of this OrientedReadId
    // sequence that will be used in this local assembly.
    void fillFirstLastOrdinalForAssembly(const Markers&, double drift);


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
