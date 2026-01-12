#pragma once

// Shasta.
#include "AnchorPair.hpp"
#include "Kmer.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MarkerInfo.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include "memory.hpp"
#include "span.hpp"



namespace shasta2 {

    class Base;
    class Marker;
    class MarkerKmers;
    class Markers;
    class MarkerInfo;
    class Reads;

    using AnchorId = uint64_t;
    class Anchor;
    class AnchorMarkerInfo;
    class Anchors;
    class AnchorInfo;
    class AnchorPairInfo;

    using AnchorBaseClass = span<const AnchorMarkerInfo>;

    class Journeys;

    string anchorIdToString(AnchorId);
    AnchorId anchorIdFromString(const string&);
}



// An Anchor is a set of AnchorMarkerInfos.
class shasta2::AnchorMarkerInfo : public MarkerInfo {
public:
    uint32_t positionInJourney = invalid<uint32_t>;

    // Default constructor.
    AnchorMarkerInfo() {}

    // Constructor from a MarkerInfo.
    AnchorMarkerInfo(
        const MarkerInfo& markerInfo) :
        MarkerInfo(markerInfo)
    {}

    // Constructor from OrientedReadId and ordinal.
    AnchorMarkerInfo(
        OrientedReadId orientedReadId,
        uint32_t ordinal) :
        MarkerInfo(orientedReadId, ordinal)
    {}

    bool operator<(const AnchorMarkerInfo& that) const
    {
        return orientedReadId < that.orientedReadId;
    }
};



class shasta2::AnchorInfo {
public:
    // The k-mer index in the MarkerKmers for the k-mer
    // that generated this anchor and its reverse complement.
    // This is only used in constructThreadFunctionPass2.
    // When using ExternalAnchors, it is not filled in.
    uint64_t kmerIndex = invalid<uint64_t>;

    AnchorInfo(uint64_t kmerIndex = invalid<uint64_t>) : kmerIndex(kmerIndex) {}
};



// An Anchor is a set of AnchorMarkerInfos.
class shasta2::Anchor : public AnchorBaseClass {
public:

    Anchor(const AnchorBaseClass& s) : AnchorBaseClass(s) {}

    void check() const;

    uint64_t coverage() const
    {
        return size();
    }
};



class shasta2::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    // Constructor to create Anchors from MarkerKmers.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        const MarkerKmers&,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        const vector<uint64_t>& maxAnchorRepeatLength,
        uint64_t threadCount);

    // Constructor to read Anchors from ExternalAnchors.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        const MarkerKmers&,
        const string& externalAnchorsName);

    // Constructor to create an empty Anchors object.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        const MarkerKmers&);

    // Constructor that makes a copy of of a source Anchors object, but removing
    // a specified set of AnchorMarkerInfos.
    // The keep vector specifies which AnchorMarkerInfos should be kept.
    // It must be of size that.anchorMarkerInfos.totalSize() and
    // is indexed by the global position of the AnchorMarkerInfo
    // in that.anchorMarkerInfos, that is, &anchorMarkerInfo-that.anchorMarkerInfos.begin().
    Anchors(
        const Anchors&,
        const string& baseName,
        const vector<bool>& keep);

    // Constructor to accesses existing Anchors.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        const MarkerKmers&,
        bool writeAccess);

    void remove();

    Anchor operator[](AnchorId) const;
    uint64_t size() const;

    // This returns the sequence of the marker k-mer
    // that this anchor was created from.
    vector<Base> anchorKmerSequence(AnchorId) const;
    Kmer anchorKmer(AnchorId) const;

    // Return the number of common oriented reads between two Anchors,
    // counting only oriented reads that have a greater ordinal on anchorId1
    // than they have on anchorId0.
    uint64_t countCommon(AnchorId anchorId0, AnchorId anchorId1) const;

    // Same as above, but also compute the average offset in bases.
    uint64_t countCommon(AnchorId anchorId0, AnchorId anchorId1, uint64_t& baseOffset) const;

    // Analyze the oriented read composition of two anchors.
    void analyzeAnchorPair(AnchorId, AnchorId, AnchorPairInfo&) const;
    void writeHtml(AnchorId, AnchorId, AnchorPairInfo&, const Journeys&
        , ostream&) const;

    // Cluster oriented reads in an anchor pair using their journey
    // portions between AnchorIdA and AnchorIdB.
    // Output to html if it is open.
    // Returns clusters in order of decreasing length.
    // Each cluster contains indices in AnchoirPair::orientedReadIds
    // of the OrientedReadIds that belong to that cluster.
    void clusterAnchorPairOrientedReads(
        const AnchorPair&,
        const Journeys&,
        double clusteringMinJaccard,
        vector< vector<uint64_t> >& clusters,
        ostream& html) const;


    void writeCoverageHistogram() const;

    MemoryMapped::VectorOfVectors<AnchorMarkerInfo, uint64_t> anchorMarkerInfos;

    // A MemoryMapped::Vector that gives, for each k-mer in the Marker K-mers,
    // the AnchordId of the first of the two anchors generated by that k-mer,
    // or invalid<AnchorId> if that k-mer did not generate any anchors.
    // This is only used in the http server.
    // When using ExternalAnchors, it is filled with invalid<AnchorId>.
    MemoryMapped::Vector<AnchorId> kmerToAnchorTable;

    // Use the above table to get the AnchorId corresponding to a given Kmer.
    // When using ExternalAnchors, this always returns invalid<AnchorId>.
    AnchorId getAnchorIdFromKmer(const Kmer&) const;

    // Get the ordinal for the AnchorMarkerInfo corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    uint32_t getOrdinal(AnchorId, OrientedReadId) const;

    // Get the positioInJourney for the AnchorMarkerInfo corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    uint32_t getPositionInJourney(AnchorId, OrientedReadId) const;

    // Get the AnchorMarkerInfo corresponding to a given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    const AnchorMarkerInfo& getAnchorMarkerInfo(AnchorId, OrientedReadId) const;


    // Find out if the given AnchorId contains the specified OrientedReadId.
    bool anchorContains(AnchorId, OrientedReadId) const;

    const string baseName;
    const Reads& reads;
    const uint64_t k;
    const uint64_t kHalf;
    const Markers& markers;
    const MarkerKmers& markerKmers;

private:

    void check() const;

public:

    // For a given AnchorId, follow the read journeys forward/backward by one step.
    // Return a vector of the AnchorIds reached in this way.
    // The count vector is the number of oriented reads each of the AnchorIds.
    void findChildren(
        const Journeys&,
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;
    void findParents(
        const Journeys&,
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;


    // In addition to the marker intervals, we also store an AnchorInfo for each Anchor.
    MemoryMapped::Vector<AnchorInfo> anchorInfos;
    void storeAnchorInfo(AnchorId)
    {
        // For now there is nothing to store.
        // AnchorInfo& anchorInfo = anchorInfos[anchorId];
    }

    // Read following.
    void followOrientedReads(
        const Journeys& journeys,
        AnchorId,
        uint64_t direction,                         // 0 = forward, 1 = backward
        uint64_t minCommonCount,
        double minJaccard,
        double minCorrectedJaccard,
        vector< pair<AnchorId, AnchorPairInfo> >&
        ) const;

    // This is fast as it uses a priority queue.
    AnchorId readFollowing(
        const Journeys& journeys,
        AnchorId,
        uint64_t direction,                         // 0 = forward, 1 = backward
        uint64_t minCommonCount,
        double aDrift,
        double bDrift
        ) const;



    // A more sophisticated version of read following with a minimum offset guarantee.
    // It finds the AnchorPair with minimum offset among all AnchorPairs that
    // satisfy the following:
    // - If direction is 0, they start at anchorId0. If direction is 1, they end at anchorId0.
    // - They have consistent offsets using the given values of aDrift, bDrift.
    // - They have at least minCommonCount oriented reads.
    // This returns false if no solution is found.
    bool readFollowing(
        const Journeys& journeys,
        AnchorId anchorId0,
        uint64_t direction,                         // 0 = forward, 1 = backward
        uint64_t minCommonCount,
        double aDrift,
        double bDrift,
        AnchorPair& anchorPair, // Filled in only when returning true.
        uint32_t& offset        // Filled in only when returning true.
        ) const;

    // Same, but can find more than one AnchorPair.
    void readFollowing(
        const Journeys& journeys,
        AnchorId anchorId0,
        uint64_t direction,                         // 0 = forward, 1 = backward
        uint64_t minCommonCount,
        uint64_t minContinueReadFollowingCount,
        double aDrift,
        double bDrift,
        vector< pair<AnchorPair, uint32_t> >&   // AnchorPairs and offsets.
        ) const;



private:

    // Data and functions used when constructing the Anchors.
    class ConstructData {
    public:
        uint64_t minAnchorCoverage;
        uint64_t maxAnchorCoverage;
        vector<uint64_t> maxAnchorRepeatLength;

        // During multithreaded pass 1 we loop over all marker k-mers
        // and for each one we find out if it can be used to generate
        // a pair of anchors or not. If it can be used,
        // we also fill in the coverage - that is,
        // the number of usable MarkerInfos that will go in each of the
        // two anchors.
        MemoryMapped::Vector<uint64_t> coverage;
    };
    ConstructData constructData;
    void constructThreadFunctionPass1(uint64_t threadId);
    void constructThreadFunctionPass2(uint64_t threadId);

};



// Information about the read composition similarity of two anchors A and B.
class shasta2::AnchorPairInfo {
public:

    // The total number of OrientedReadIds in each of the anchors A and B.
    uint64_t totalA = 0;
    uint64_t totalB = 0;

    // The number of common oriented reads.
    uint64_t common = 0;

    // The number of oriented reads present in A but not in B.
    uint64_t onlyA = 0;

    // The number of oriented reads present in B but not in A.
    uint64_t onlyB = 0;

    // The rest of the statistics are only valid if the number
    // of common oriented reads is not 0.

    // The estimated offset between the two Anchors.
    // The estimate is done using the common oriented reads.
    int64_t offsetInMarkers = invalid<int64_t>;
    int64_t offsetInBases = invalid<int64_t>;

    // The number of onlyA reads which are too short to be on edge B,
    // based on the above estimated offset.
    uint64_t onlyAShort = invalid<uint64_t>;

    // The number of onlyB reads which are too short to be on edge A,
    // based on the above estimated offset.
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCount() const
    {
        return common;
    }
    uint64_t unionCount() const {
        return totalA + totalB - common;
    }
    uint64_t correctedUnionCount() const
    {
        return unionCount() - onlyAShort - onlyBShort;
    }
    double jaccard() const
    {
        return double(intersectionCount()) / double(unionCount());
    }
    double correctedJaccard() const
    {
        return double(intersectionCount()) / double(correctedUnionCount());
    }

    void reverse()
    {
        swap(totalA, totalB);
        swap(onlyA, onlyB);
        swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }

};
