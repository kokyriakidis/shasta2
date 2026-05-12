#pragma once

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "iosfwd.hpp"
#include "memory.hpp"
#include "vector.hpp"

namespace shasta2 {
    class LocalAssembly6;

    class Anchors;
    class Base;
}



class shasta2::LocalAssembly6 {
public:

    // This assembles sequence between two anchors using an input vector
    // of oriented reads, which must be sorted.
    // * Of the oriented reads given on input, only the ones that appear
    //   in at least one of the two anchors can be used for assembly.
    // * At least one of the input oriented reads must appear on both anchors.
    // When this completes successfully, the assembled sequence is stored
    // in the sequence vector below. If an error occurs,
    // this throws a std::runtime_error.
    LocalAssembly6(
        const Anchors&,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        ostream& html,
        bool debug,
        const vector<OrientedReadId>& orientedReadIds);

    // Assembled sequence and its coverage.
    vector<Base> sequence;

private:

    const Anchors& anchors;
    AnchorId anchorIdA;     // Left Anchor.
    AnchorId anchorIdB;     // Right Anchor.
    ostream& html;          // Pass ostream(0) ti suppress html output.
    bool debug;



    // EXPOSE WHEN CODE STABILIZES.

    // For reads fixed on one side only, we use a sequence length
    // equal to aDrift * offset + bDrift.
    const double aDrift = 0.1;
    const double bDrift = 30.;



    // The oriented reads used in this local assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The ordinals of this OrientedReadId in the two Anchors, if any.
        uint32_t ordinalA = invalid<uint32_t>;
        uint32_t ordinalB = invalid<uint32_t>;
        uint32_t ordinalOffsetAB() const
        {
            SHASTA2_ASSERT(isOnBothAnchors());
            SHASTA2_ASSERT(ordinalB > ordinalA);
            return ordinalB - ordinalA;
        }

        // The base positions of this OrientedReadId in the two Anchors, if any.
        // These are positions in the oriented read sequence
        // of the mid point of the marker.
        uint32_t positionA = invalid<uint32_t>;
        uint32_t positionB = invalid<uint32_t>;
        uint32_t positionOffsetAB() const
        {
            SHASTA2_ASSERT(isOnBothAnchors());
            SHASTA2_ASSERT(positionB > positionA);
            return positionB - positionA;
        }

        // Whether this OrientedReadId appears in the two Anchors.
        bool isOnAnchorA() const {return isValid(ordinalA);}
        bool isOnAnchorB() const {return isValid(ordinalB);}
        bool isOnBothAnchors() const
        {
            return isOnAnchorA() and isOnAnchorB();
        }

        // The above portion is filled in by gatherOrientedReads().
        // The rest is filled in by gatherOrientedReadSequences().

        // The region of this oriented read that will be used in this local assembly.
        // The positions are positions of the marker midpoints.
        uint32_t positionBegin = invalid<uint32_t>;
        uint32_t positionEnd = invalid<uint32_t>;

        uint32_t positionOffsetForAssembly() const
        {
            return positionEnd - positionBegin;
        }

        uint64_t sequenceId = invalid<uint64_t>;

    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads(const vector<OrientedReadId>&);
    void gatherOrientedReadsSequences();
    void writeOrientedReads() const;

    uint32_t offset;
    void estimateOffset();



    // The distinct sequences to be used for assembly.
    // They get sorted in order of decreasing coverage.
    class SequenceInfo {
    public:
        shared_ptr< vector<Base> > sequencePointer;
        vector<OrientedReadId> orientedReadIds;
        uint64_t id = invalid<uint64_t>;
        bool operator< (const SequenceInfo& that) const
        {
            return orientedReadIds.size() > that.orientedReadIds.size();
        }
        SequenceInfo(OrientedReadId orientedReadId, const vector<Base>& sequence) :
            sequencePointer(make_shared< vector<Base> >(sequence)),
            orientedReadIds(1, orientedReadId)
        {}
    };
    vector<SequenceInfo> fixedSequencesTable;
    vector<SequenceInfo> leftFixedSequencesTable;
    vector<SequenceInfo> rightFixedSequencesTable;
    void writeSequenceTable(
        const vector<SequenceInfo>& sequenceTable,
        bool fixedOnA,
        bool fixedOnB) const;

    void assemble();

};
