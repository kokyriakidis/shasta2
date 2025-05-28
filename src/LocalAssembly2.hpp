#pragma once

#include "Anchor.hpp"

namespace shasta {
    class LocalAssembly2;
    class AnchorPair;
}


class shasta::LocalAssembly2 {
public:

    // This constructor uses an AnchorPair between the two specified anchors.
    LocalAssembly2(
        const Anchors&,
        ostream& html,
        bool debug,
        double aDrift,
        double bDrift,
        AnchorId anchorIdA,
        AnchorId anchorIdB);

    // This constructor uses the AnchorPair passed in.
    LocalAssembly2(
        const Anchors&,
        ostream& html,
        bool debug,
        double aDrift,
        double bDrift,
        const AnchorPair&);

    // After calling one of the above constructors, call run.
    void run(
        bool computeAlignment,
        uint64_t maxAbpoaLength);
    void getSequence(vector<Base>&) const;

private:
    // This checks for shortcuts (e. g. all sequences are identical)
    // that can be used to run faster. It returns true if successful.
    bool runFast(bool computeAlignment);
    void writeOrientedReadsFast() const;

    const Anchors& anchors;
    ostream& html;
    bool debug;

    // A table of the distinct marker k-mers involved in this LocalAssembly2.
    // This includes the marker k-mers between ordinalA (included) and ordinalB (included)
    // for all oriented reads used in this LocalAssembly2.
    // The index in this vector provides a perfect hash function for these k-mers.
    vector<Kmer> kmers;

    class MarkerInfo {
    public:
        uint32_t ordinal;
        uint32_t position;          // The position in the oriented read of the first base of this marker.
        uint64_t kmerId;            // The kmerId (index in the kmers vector) for this marker k-mer.
        Kmer kmer;
    };



    // Information for each oriented read used in this local assembly.
    // We use all common oriented reads with positive ordinal offset
    // between anchorIdA and anchorIdB.
    // Stored sorted by OrientedReadId.
    class OrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinalA;
        uint32_t ordinalB;
        vector<MarkerInfo> markerInfos;
        OrientedRead(OrientedReadId orientedReadId, uint32_t ordinalA, uint32_t ordinalB) :
            orientedReadId(orientedReadId), ordinalA(ordinalA), ordinalB(ordinalB) {}

        uint32_t ordinalOffset() const
        {
            return ordinalB - ordinalA;
        }
        uint32_t positionA() const
        {
            return markerInfos.front().position;
        }
        uint32_t positionB() const
        {
            return markerInfos.back().position;
        }
        uint32_t sequenceLength() const
        {
            return positionB() - positionA();
        }

        // Get the MarkerInfo for a given ordinal.
        const MarkerInfo& getMarkerInfo(uint32_t ordinal) const
        {
            SHASTA_ASSERT(ordinal >= ordinalA);
            SHASTA_ASSERT(ordinal <= ordinalB);
            const MarkerInfo& markerInfo = markerInfos[ordinal - ordinalA];
            SHASTA_ASSERT(markerInfo.ordinal == ordinal);
            return markerInfo;
        }


        // Work vector used by split.
        // It is a vector of pairs (kmerId, ordinal) for the markers
        // internal to the interval being split.
        vector< pair<uint64_t, uint32_t> > internalMarkers;

        // Same as above, but only including marker k-mers
        // that are appear exactly once in each oriented read.
        vector< pair<uint64_t, uint32_t> > commonUniqueInternalMarkers;

    };
    vector<OrientedRead> orientedReads;
    // This does not fill in the markerInfos.
    void gatherOrientedReads(
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        double aDrift,
        double bDrift);
    void gatherOrientedReads(
        const AnchorPair&,
        double aDrift,
        double bDrift);

    // This assumes that gatherKmers has already been called.
    void writeOrientedReads();

    // This gathers the marker k-mers of all reads and fills in
    // the kmers vector and the markerInfos of each OrientedRead.
    void gatherKmers();


    class AlignedMarkers {
    public:
        vector<uint32_t> ordinals;
    };
    std::list<AlignedMarkers> allAlignedMarkers;
    vector<AlignedMarkers> allAlignedMarkersVector; // For convenience during assembly.
    void alignMarkers();
    void split(const AlignedMarkers&, const AlignedMarkers&, vector<AlignedMarkers>&);

    // The call to alignMarkers can throw Failure, which contains a bit vector
    // that says which OrientedReads we should keep.
    // The remaining OrientedReads are removed and the process is restarted.
    class Failure {
    public:
        vector<bool> keep;
    };


    // Sequence assembly.
    // Each assembly step does a MSA between a pair of adjacent AlignedMarkers.
    void assemble(bool computeAlignment, uint64_t maxAbpoaLength);                   // Assemble all steps.
    void assemble(bool computeAlignment, uint64_t maxAbpoaLength, uint64_t step);    // Assemble one step.
    void writeConsensus(const vector< pair<Base, uint64_t> >& stepConsensus) const;
    void writeAlignment(
        const vector< vector<Base> >& inputSequences,
        const vector< pair<Base, uint64_t> >& consensus,
        const vector< vector<AlignedBase> >& alignment,
        const vector<AlignedBase>& alignedConsensus) const;
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    void writeConsensus() const;
    void writeAlignment();
};
