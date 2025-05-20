// Shasta2
#include "LocalAssembly2.hpp"
#include "abpoaWrapper.hpp"
#include "AnchorPair.hpp"
#include "deduplicate.hpp"
#include "extractKmer128.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "poastaWrapper.hpp"
#include "Reads.hpp"
using namespace shasta;

// Standard library.
#include <chrono.hpp>
#include <fstream.hpp>
#include <queue>



LocalAssembly2::LocalAssembly2(
    const Anchors& anchors,
    ostream& html,
    bool debugArgument,
    double aDrift,
    double bDrift,
    AnchorId anchorIdA,
    AnchorId anchorIdB) :
    anchors(anchors),
    html(html),
    debug(debugArgument)
{
    gatherOrientedReads(anchorIdA, anchorIdB, aDrift, bDrift);
}



LocalAssembly2::LocalAssembly2(
    const Anchors& anchors,
    ostream& html,
    bool debugArgument,
    double aDrift,
    double bDrift,
    const AnchorPair& anchorPair) :
    anchors(anchors),
    html(html),
    debug(debugArgument)
{
    gatherOrientedReads(anchorPair, aDrift, bDrift);
}



void LocalAssembly2::run(
    bool computeAlignment,
    uint64_t maxAbpoaLength) {

    // Iterate until alignMarkers is successful.
    // Each failed iteration removes one or more OrientedReads.
    while(true) {
        gatherKmers();

        if(html) {
            writeOrientedReads();
        }

        try {
            alignMarkers();
        } catch(const Failure& failure) {

            // We must remove some OrientedReads.
            vector<OrientedRead> newOrientedReads;
            uint64_t discardedCount = 0;
            for(uint64_t i=0; i<orientedReads.size(); i++) {
                const OrientedRead& orientedRead = orientedReads[i];
                if(failure.keep[i]) {
                    newOrientedReads.emplace_back(orientedRead.orientedReadId, orientedRead.ordinalA, orientedRead.ordinalB);
                } else {
                    if(html) {
                        ++discardedCount;
                        html << "<br>Discarding " << orientedRead.orientedReadId;
                    }
                }
            }
            orientedReads.swap(newOrientedReads);
            allAlignedMarkers.clear();

            // Try again.
            if(html) {
                html << "<br>Restarting after discarding " << discardedCount << " oriented reads.";
            }
            continue;
        }

        // Success.
        break;
    }

    // Assemble sequence using the AlignedMarkers we found.
    assemble(computeAlignment, maxAbpoaLength);

}



// This finds the OrientedReads to be used in this LocalAssembly2.
// We use all common oriented reads with positive ordinal offset
// between anchorIdA and anchorIdB.
// This does not fill in the markerInfos.
void LocalAssembly2::gatherOrientedReads(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    double aDrift,
    double bDrift)
{
    // Create an AnchorPair that contains the OrientedReadIds we want.
    // The last argument means that we just require the oriented reads
    // to visit anchorIdB after anchorIdA, not immediately after.
    const AnchorPair anchorPairBeforeSplitting(anchors, anchorIdA, anchorIdB, false);

    gatherOrientedReads(anchorPairBeforeSplitting, aDrift, bDrift);
}



void LocalAssembly2::gatherOrientedReads(
    const AnchorPair& anchorPairBeforeSplitting,
    double aDrift,
    double bDrift)
{


    // Split this AnchorPair into AnchorPairs with consistent offsets
    // and only use the larger one.
    vector<AnchorPair> splitAnchorPairs;
    anchorPairBeforeSplitting.split(anchors, aDrift, bDrift, splitAnchorPairs);
    const AnchorPair& anchorPair = splitAnchorPairs.front();

    const uint64_t n = anchorPair.orientedReadIds.size();

    // Also get the ordinals of each oriented read in anchorIdA and anchorIdB.
    vector<pair<uint32_t, uint32_t> > ordinals;
    anchorPair.getOrdinals(anchors, ordinals);
    SHASTA_ASSERT(ordinals.size() == n);

    // Store the OrientedReads. This does not fill in the markerInfos.
    orientedReads.clear();
    for(uint64_t i=0; i<n; i++) {
        const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
        const auto& p = ordinals[i];
        const uint32_t ordinalA = p.first;
        const uint32_t ordinalB = p.second;
        orientedReads.emplace_back(orientedReadId, ordinalA, ordinalB);
    }
}



// This gathers the marker k-mers of all reads and fills in
// the kmers vector and the markerInfos of each OrientedRead.
void LocalAssembly2::gatherKmers()
{
    const uint64_t k = anchors.k;

    // Loop over the OrientedReads to fill in their MarkerInfos
    // and also to fill the kmers vector.
    for(OrientedRead& orientedRead: orientedReads) {

        // Access the information we need for this OrientedRead.
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const auto read = anchors.reads.getRead(readId);
        const uint64_t readLength = read.baseCount;
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

        // Loop over all ordinals in [ordinalA, ordinalB].
        orientedRead.markerInfos.resize(orientedRead.ordinalB + 1 - orientedRead.ordinalA);
        for(uint32_t ordinal=orientedRead.ordinalA; ordinal<=orientedRead.ordinalB; ordinal++) {
            MarkerInfo& markerInfo = orientedRead.markerInfos[ordinal - orientedRead.ordinalA];
            markerInfo.ordinal = ordinal;
            const Marker& marker = orientedReadMarkers[ordinal];
            markerInfo.position = marker.position;

            if(strand == 0) {
                extractKmer128(read, markerInfo.position, k, markerInfo.kmer);
            } else {
                extractKmer128(read, readLength - k - markerInfo.position, k, markerInfo.kmer);
                markerInfo.kmer = markerInfo.kmer.reverseComplement(k);
            }

            kmers.push_back(markerInfo.kmer);
        }
    }
    deduplicate(kmers);

    // Now we can fill in the kmerId fields in the MarkerInfos.
    for(OrientedRead& orientedRead: orientedReads) {
        for(MarkerInfo& markerInfo: orientedRead.markerInfos) {
            auto it = std::lower_bound(kmers.begin(), kmers.end(), markerInfo.kmer);
            SHASTA_ASSERT(it != kmers.end());
            SHASTA_ASSERT(*it == markerInfo.kmer);
            markerInfo.kmerId = it - kmers.begin();
        }
    }
}



// This assumes that gatherKmers has already been called.
void LocalAssembly2::writeOrientedReads()
{
    html <<
        "<h3>Oriented reads portions used in this local assembly</h3>"
        "<p>This local assembly uses " << orientedReads.size() << " oriented reads."
        "<table><tr>"
        "<th>OrientedReadId"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Ordinal<br>offset"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Sequence<br>length";

    for(const OrientedRead& orientedRead: orientedReads) {
        html <<
            "<tr>"
            "<td class=centered>" << orientedRead.orientedReadId <<
            "<td class=centered>" << orientedRead.ordinalA <<
            "<td class=centered>" << orientedRead.ordinalB <<
            "<td class=centered>" << orientedRead.ordinalOffset() <<
            "<td class=centered>" << orientedRead.positionA() <<
            "<td class=centered>" << orientedRead.positionB() <<
            "<td class=centered>" << orientedRead.sequenceLength();
    }
    html << "</table>";
}


#if 0
void LocalAssembly2::gatherOrientedReadsKmers()
{
    const uint64_t k = anchors.k;

    uint64_t totalKmerCount = 0;
    for(OrientedRead& orientedRead: orientedReads) {

        // Access the information we need for this OrientedRead.
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const auto read = anchors.reads.getRead(readId);
        const uint64_t readLength = read.baseCount;
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

        // Loop over marker k-mers internal to (ordinalA, ordinalB).
        orientedRead.kmers.clear();
        for(uint32_t ordinal=orientedRead.ordinalA + 1; ordinal!= orientedRead.ordinalB; ++ordinal) {
            const Marker& marker = orientedReadMarkers[ordinal];
            const uint64_t position = marker.position;

            Kmer kmer;
            if(strand == 0) {
                extractKmer128(read, position, k, kmer);
            } else {
                extractKmer128(read, readLength - k - position, k, kmer);
                kmer = kmer.reverseComplement(k);
            }

            orientedRead.kmers.push_back(make_pair(kmer, ordinal));
        }

        SHASTA_ASSERT(orientedRead.kmers.size() == orientedRead.ordinalB - orientedRead.ordinalA - 1);
        totalKmerCount += orientedRead.kmers.size();
    }

    if(html) {
        html << "<p>Total number of internal marker k-mers on all reads is " << totalKmerCount;
    }

}



void LocalAssembly2::writeOrientedReads() const
{
    html <<
        "<h3>Oriented read sequences</h3>"
        "<p>This local assembly will use the following "<< orientedReads.size() << " oriented reads, "
        "ordered by decreasing sequence length."
        "<table>"
        "<tr>"
        "<th>OrientedReadId"
        "<th>PositionA<br>in journey<th>PositionB<br>in journey<th>Journey<br>offset"
        "<th>OrdinalA<th>OrdinalB<th>Ordinal<br>offset"
        "<th>PositionA<th>PositionB<th>Sequence<br>length";

    for(const OrientedRead& orientedRead: orientedReads) {
        html <<
            "<tr>"
            "<td class=centered>" << orientedRead.orientedReadId <<
            "<td class=centered>" << orientedRead.positionInJourneyA <<
            "<td class=centered>" << orientedRead.positionInJourneyB <<
            "<td class=centered>" << orientedRead.positionInJourneyB - orientedRead.positionInJourneyA <<
            "<td class=centered>" << orientedRead.ordinalA <<
            "<td class=centered>" << orientedRead.ordinalB <<
            "<td class=centered>" << orientedRead.ordinalB - orientedRead.ordinalA <<
            "<td class=centered>" << orientedRead.basePositionA <<
            "<td class=centered>" << orientedRead.basePositionB <<
            "<td class=centered>" << orientedRead.basePositionB - orientedRead.basePositionA;
    }

    html << "</table>";
}



void LocalAssembly2::gatherKmers()
{
    kmers.clear();
    for(const OrientedRead& orientedRead: orientedReads) {
        for(const auto& p: orientedRead.kmers) {
            kmers.push_back(p.first);
        }
    }
    deduplicate(kmers);

    if(html) {
        if(html) {
            html << "<p>Total number of distinct internal marker k-mers on all reads is " << kmers.size();
        }
    }

    // Fill in the k-mer id in the OrientedReads.
    for(OrientedRead& orientedRead: orientedReads) {
        for(auto& p: orientedRead.kmers) {
            const Kmer& kmer = p.first;
            auto it = std::lower_bound(kmers.begin(), kmers.end(), kmer);
            SHASTA_ASSERT(it != kmers.end());
            SHASTA_ASSERT(*it == kmer);
            p.second = it - kmers.begin();
        }

        orientedRead.kmerFrequency.clear();
        orientedRead.kmerFrequency.resize(kmers.size(), 0);
        for(const auto& p: orientedRead.kmers) {
            const uint64_t kmerId = p.second;
            ++orientedRead.kmerFrequency[kmerId];
        }
    }

    if(html) {
        html <<
            "<h3>Marker k-mers</h3>"
            "<table><tr><th<OrientedReadId<th class=left>Marker k-mers";
        for(OrientedRead& orientedRead: orientedReads) {
            html << "<tr><td class=centered>" << orientedRead.orientedReadId << "<td class=left>";
            for(auto& p: orientedRead.kmers) {
                html << p.second << " ";
            }
        }
        html << "</table>";

        html << "<h3>Marker k-mer frequency</h3>"
            "<table><tr><th>K-mer id";
        for(OrientedRead& orientedRead: orientedReads) {
            html << "<th>" << orientedRead.orientedReadId;
        }
        for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
            html << "<tr><th>" << kmerId;
            for(OrientedRead& orientedRead: orientedReads) {
                html << "<td class=centered>" << orientedRead.kmerFrequency[kmerId];
            }
        }
        html << "</table>";
    }

}
#endif




void LocalAssembly2::alignMarkers()
{

    // Start with two AlignedMarkers at A and B.
    AlignedMarkers alignedMarkersA;
    for(const OrientedRead& orientedRead: orientedReads) {
        alignedMarkersA.ordinals.push_back(orientedRead.ordinalA);
    }
    AlignedMarkers alignedMarkersB;
    for(const OrientedRead& orientedRead: orientedReads) {
        alignedMarkersB.ordinals.push_back(orientedRead.ordinalB);
    }
    allAlignedMarkers.push_back(alignedMarkersA);
    allAlignedMarkers.push_back(alignedMarkersB);



    // Recursively split the interval between adjacent AlignedMarkers.
    using ListIterator = std::list<AlignedMarkers>::iterator;
    std::queue<ListIterator> q;
    q.push(allAlignedMarkers.begin());
    while(not q.empty()) {
        const ListIterator it0 = q.front();
        q.pop();
        ListIterator it1 = it0;
        ++it1;

        // Try and split between alignedMarkers0 and alignedMarkers1.
        const AlignedMarkers& alignedMarkers0 = *it0;
        const AlignedMarkers& alignedMarkers1 = *it1;
        vector<AlignedMarkers> allNewAlignedMarkers;
        split(alignedMarkers0, alignedMarkers1, allNewAlignedMarkers);
        for(const AlignedMarkers& alignedMarkers: allNewAlignedMarkers) {
            const ListIterator it = allAlignedMarkers.insert(it1, alignedMarkers);
            q.push(it);
        }
    }

    // Sanity check: on each AlignedMarkers, all oriented reads must have the same kmerId.
    for(const AlignedMarkers& alignedMarkers: allAlignedMarkers) {
        for(uint64_t i=0; i<orientedReads.size(); i++) {
            const OrientedRead& orientedRead = orientedReads[i];
            const uint32_t ordinal = alignedMarkers.ordinals[i];
            const MarkerInfo& markerInfo = orientedRead.getMarkerInfo(ordinal);
            const OrientedRead& firstOrientedRead = orientedReads.front();
            const uint32_t firstOrientedReadOrdinal = alignedMarkers.ordinals.front();
            const MarkerInfo& firstOrientedReadMarkerInfo = firstOrientedRead.getMarkerInfo(firstOrientedReadOrdinal);
            SHASTA_ASSERT(markerInfo.kmerId == firstOrientedReadMarkerInfo.kmerId);
        }
    }

    // For convenience also store the AlignedMarkers in a vector.
    allAlignedMarkersVector.clear();;
    copy(allAlignedMarkers.begin(), allAlignedMarkers.end(), back_inserter(allAlignedMarkersVector));
    allAlignedMarkers.clear();


    // Write out allAlignedMarkers.
    if(html and debug) {
        html <<
            "<h3>Aligned markers</h3>"
            "<p>Found " << allAlignedMarkersVector.size() << " sets of aligned markers, "
            "shown one per row in the following table."
            "<table><tr><th>Step";
        for(const OrientedRead& orientedRead: orientedReads) {
            html << "<th>" << orientedRead.orientedReadId;
        }
        for(uint64_t step=0; step<allAlignedMarkersVector.size(); step++) {
            const AlignedMarkers& alignedMarkers = allAlignedMarkersVector[step];
            html << "<tr><th>" << step;
            for(uint64_t i=0; i<orientedReads.size(); i++) {
                const uint32_t ordinal = alignedMarkers.ordinals[i];
                html << "<td class=centered>" << ordinal;

            }
        }
        html << "</table>";
    }
}



void LocalAssembly2::split(
    const AlignedMarkers& alignedMarkers0,
    const AlignedMarkers& alignedMarkers1,
    vector<AlignedMarkers>& newAlignedMarkers)
{
    newAlignedMarkers.clear();

    // If all OrientedReads have no internal markers, do nothing.
    bool internalMarkersExist = false;
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const uint32_t ordinal0 = alignedMarkers0.ordinals[i];
        const uint32_t ordinal1 = alignedMarkers1.ordinals[i];
        if(ordinal1 - ordinal0 > 1) {
            internalMarkersExist = true;
            break;
        }
    }
    if(not internalMarkersExist) {
        return;
    }


    if(debug and html) {
        html <<
            "<h3>Split operation</h3>"
            "<table><tr><th>OrientedReadId<th>Ordinal0<th>Ordinal1";
        for(uint64_t i=0; i<orientedReads.size(); i++) {
            html <<
                "<tr>"
                "<td class=centered>" << orientedReads[i].orientedReadId <<
                "<td class=centered>" << alignedMarkers0.ordinals[i] <<
                "<td class=centered>" << alignedMarkers1.ordinals[i];
        }
        html << "</table>";
    }



    // For each OrientedRead, gather marker k-mers internal to the (0, 1) interval we want to split.
    // For each OrientedRead, only keep the ones that appear once.
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        OrientedRead& orientedRead = orientedReads[i];
        orientedRead.internalMarkers.clear();
        const uint32_t ordinal0 = alignedMarkers0.ordinals[i];
        const uint32_t ordinal1 = alignedMarkers1.ordinals[i];
        for(uint32_t ordinal=ordinal0+1; ordinal<ordinal1; ordinal++) {
            const MarkerInfo& markerInfo = orientedRead.getMarkerInfo(ordinal);
            orientedRead.internalMarkers.push_back(make_pair(markerInfo.kmerId, ordinal));
        }
        if(debug and html) {
            html << "<br>" << orientedRead.orientedReadId << " has " <<
                orientedRead.internalMarkers.size() << " internal markers.";
        }
    }

    // Only keep the k-mers that appear exactly once in each read.
    for(OrientedRead& orientedRead: orientedReads) {
        deduplicateAndCountAndKeepUnique(orientedRead.internalMarkers, OrderPairsByFirstOnly<uint64_t, uint32_t>());
        if(debug and html) {
            html << "<br>" << orientedRead.orientedReadId << " has " <<
                orientedRead.internalMarkers.size() << " unique internal marker k-mers.";
        }
    }

    // Count how many oriented reads each of those k-mers appear in.
    vector<uint64_t> uniqueKmerIds;
    for(OrientedRead& orientedRead: orientedReads) {
        for(const auto& p: orientedRead.internalMarkers) {
            const uint64_t kmerId = p.first;
            uniqueKmerIds.push_back(kmerId);
        }
    }
    vector<uint64_t> frequency;
    deduplicateAndCount(uniqueKmerIds, frequency);

    // Find the marker k-mers that appear exactly once in each oriented read.
    vector<uint64_t> commonUniqueKmerIds;
    for(uint64_t i=0; i<uniqueKmerIds.size(); i++) {
        if(frequency[i] == orientedReads.size()) {
            commonUniqueKmerIds.push_back(uniqueKmerIds[i]);
        }
    }
    if(debug and html) {
        html << "<br>Found " << commonUniqueKmerIds.size() <<
            " common unique internal marker k-mers.";
    }

    // Store the commonUniqueInternalMarkers in each OrientedRead.
    for(OrientedRead& orientedRead: orientedReads) {
        orientedRead.commonUniqueInternalMarkers.clear();
        for(const auto& p: orientedRead.internalMarkers) {
            const uint64_t kmerId = p.first;
            if(binary_search(commonUniqueKmerIds.begin(), commonUniqueKmerIds.end(), kmerId)) {
                orientedRead.commonUniqueInternalMarkers.push_back(p);
            }
        }

        // Sort them by ordinal.
        sort(orientedRead.commonUniqueInternalMarkers.begin(), orientedRead.commonUniqueInternalMarkers.end(),
            OrderPairsBySecondOnly<uint64_t, uint32_t>());
    }

    // All OrientedReads must have exactly the same number of commonUniqueInternalMarkers.
    const uint64_t commonUniqueInternalMarkersCount = orientedReads.front().commonUniqueInternalMarkers.size();
    for(const OrientedRead& orientedRead: orientedReads) {
        SHASTA_ASSERT(orientedRead.commonUniqueInternalMarkers.size() == commonUniqueInternalMarkersCount);
    }
    if(debug and html) {
        html << "<br>Found " << commonUniqueInternalMarkersCount << " common internal unique marker k-mers.";

        html << "<table>";
        for(const OrientedRead& orientedRead: orientedReads) {
            html << "<tr><th>" << orientedRead.orientedReadId;
            for(const auto& p: orientedRead.commonUniqueInternalMarkers)
            {
                html << "<td class=centered>" << p.first << " " << p.second;
            }
        }
        html << "<p></table>";

        // Write the table again, without the ordinal.
        html << "<table>";
        for(const OrientedRead& orientedRead: orientedReads) {
            html << "<tr><th>" << orientedRead.orientedReadId;
            for(const auto& p: orientedRead.commonUniqueInternalMarkers)
            {
                html << "<td class=centered>" << p.first;
            }
        }
        html << "</table>";
    }



    // In the easy and most common case, the common unique marker k-mers appear
    // in the same order in all OrientedReads.
    // Check if this is the case.
    bool isEasyCase = true;
    for(const OrientedRead& orientedRead: orientedReads) {
        for(uint64_t i=0; i<commonUniqueInternalMarkersCount; i++) {
            if(orientedRead.commonUniqueInternalMarkers[i].first !=
                orientedReads.front().commonUniqueInternalMarkers[i].first) {
                isEasyCase = false;
                break;
            }
        }
        if(not isEasyCase) {
            break;
        }
    }
    if(debug and html) {
        if(isEasyCase) {
            html << "<br>All common unique internal marker k-mers appear in the same order in all oriented reads.";
        } else {
            html << "<br>Common unique internal marker k-mers don't appear in the same order in all oriented reads.";
        }
    }


    // If not in the easy case, we flag some OrientedReadsa for removal and throw a Failure.
    if(not isEasyCase) {

        // Group identical sequences of common unique internal markers.
        std::map<vector<uint64_t>, vector<uint64_t> > m;
        for(uint64_t i=0; i<orientedReads.size(); i++) {
            vector<uint64_t> kmerIds;
            for(auto& p: orientedReads[i].commonUniqueInternalMarkers) {
                const uint64_t kmerId = p.first;
                kmerIds.push_back(kmerId);
            }
            m[kmerIds].push_back(i);
        }

        // Gather groups of OrientedReads with identical sequences of common unique internal markers..
        vector< vector<uint64_t> > groups;
        for(const auto& p: m) {
            groups.push_back(p.second);
        }
        sort(groups.begin(), groups.end(), OrderVectorsByDecreasingSize<uint64_t>());

        if(debug and html) {
            html << "<p>Groups of consistent oriented reads:";
            for(const vector<uint64_t>& group: groups) {
                html << "<br>";
                for(const uint64_t i: group) {
                    html << orientedReads[i].orientedReadId << " ";
                }
            }
            html << "<p>Only the first and largest group will be kept.";
        }

        Failure failure;
        failure.keep.resize(orientedReads.size(), false);
        for(const uint64_t i: groups.front()) {
            failure.keep[i] = true;
        }
        throw failure;
    }



    // Each common unique internal marker generates a new AlignedMarkers.
    newAlignedMarkers.clear();
    for(uint64_t i=0; i<commonUniqueInternalMarkersCount; i++) {
        AlignedMarkers alignedMarkers;
        for(OrientedRead& orientedRead: orientedReads) {
            orientedRead.internalMarkers.clear();
            orientedRead.commonUniqueInternalMarkers.clear();
            alignedMarkers.ordinals.push_back(orientedRead.commonUniqueInternalMarkers[i].second);
        }
        newAlignedMarkers.push_back(alignedMarkers);
    }

    if(debug and html) {
        html << "Found " << newAlignedMarkers.size() << " new aligned markers.";
        html << "<table><tr>";
        for(const OrientedRead& orientedRead: orientedReads) {
            html << "<th>" << orientedRead.orientedReadId;
        }
        for(const AlignedMarkers& alignedMarkers: newAlignedMarkers) {
            html << "<tr>";
            for(const uint32_t ordinal: alignedMarkers.ordinals) {
                html << "<td class=centered>" << ordinal;
            }
        }
        html << "</table>";

    }
}



void LocalAssembly2::assemble(bool computeAlignment, uint64_t maxAbpoaLength)
{
    consensus.clear();
    if(computeAlignment) {
        alignment.clear();
        alignment.resize(orientedReads.size());
        alignedConsensus.clear();
    }
    for(uint64_t step=0; step<allAlignedMarkersVector.size()-1; step++) {
        assemble(computeAlignment, maxAbpoaLength, step);
    }

    if(html) {
        writeConsensus();
        if(computeAlignment) {
            writeAlignment();
        }
    }
}



void LocalAssembly2::assemble(
    bool computeAlignment,
    uint64_t maxAbpoaLength,
    uint64_t step)
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    const AlignedMarkers& alignedMarkers0 = allAlignedMarkersVector[step];
    const AlignedMarkers& alignedMarkers1 = allAlignedMarkersVector[step + 1];



    // Gather the sequences to be used for the MSA of this step.
    vector< vector<Base> > inputSequences(orientedReads.size());
    uint64_t maxSequenceLength = 0;
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        OrientedRead& orientedRead = orientedReads[i];
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;
        const uint32_t ordinal0 = alignedMarkers0.ordinals[i];
        const uint32_t ordinal1 = alignedMarkers1.ordinals[i];
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
        const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
        const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;

         for(uint32_t position=position0; position!=position1; position++) {
             inputSequences[i].push_back(anchors.reads.getOrientedReadBase(orientedReadId, position));
        }
         maxSequenceLength = max(maxSequenceLength, inputSequences[i].size());
    }

    if(debug and html) {
        html <<
            "<h3>Assembly step " << step << "</h3>"
            "<table><tr>"
            "<th>OrientedReadId"
            "<th>Ordinal0"
            "<th>Ordinal1"
            "<th>Position0"
            "<th>Position1"
            "<th>Sequence<br>length"
            "<th class=left>Sequence";

        for(uint64_t i=0; i<orientedReads.size(); i++) {
            OrientedRead& orientedRead = orientedReads[i];
            const OrientedReadId orientedReadId = orientedRead.orientedReadId;
            const uint32_t ordinal0 = alignedMarkers0.ordinals[i];
            const uint32_t ordinal1 = alignedMarkers1.ordinals[i];
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
            const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
            const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;

            html <<
                "<tr>"
                "<td class=centered>" << orientedReadId <<
                "<td class=centered>" << ordinal0 <<
                "<td class=centered>" << ordinal1 <<
                "<td class=centered>" << position0 <<
                "<td class=centered>" << position1 <<
                "<td class=centered>" << position1 - position0 <<
                "<td style='font-family:monospace'>";
            copy(inputSequences[i].begin(), inputSequences[i].end(), ostream_iterator<Base>(html));
        }
        html << "</table>";
    }

    // Do the MSA.
    const bool usePoasta = (maxSequenceLength > maxAbpoaLength);
    vector< pair<Base, uint64_t> > stepConsensus;
    vector< vector<AlignedBase> > stepAlignment;
    vector<AlignedBase> stepAlignedConsensus;
    const auto t0 = steady_clock::now();
    if(usePoasta) {
        poasta(inputSequences, stepConsensus, stepAlignment, stepAlignedConsensus);
    } else {
        abpoa(inputSequences, stepConsensus, stepAlignment, stepAlignedConsensus, computeAlignment);
    }
    const auto t1 = steady_clock::now();
    if(debug and html) {
        html << "<p>MSA with " << (usePoasta ? "poasta" : "abpoa") << " took " << seconds(t1-t0) << " s.";
        writeConsensus(stepConsensus);
        if(computeAlignment) {
            writeAlignment(inputSequences, stepConsensus, stepAlignment, stepAlignedConsensus);
        }
    }

    // Append the consensus and the alignment for this step to the global consensus and alignment.
    const uint64_t begin = consensus.size();
    copy(stepConsensus.begin(), stepConsensus.end(), back_inserter(consensus));
    const uint64_t end = consensus.size();
    if(debug and html) {
        html << "<p>This assembly step contributed positions " << begin << "-" << end <<
            " of the global consensus.";
    }
    if(computeAlignment) {
        for(uint64_t i=0; i<orientedReads.size(); i++) {
            copy(stepAlignment[i].begin(), stepAlignment[i].end(), back_inserter(alignment[i]));
        }
        copy(stepAlignedConsensus.begin(), stepAlignedConsensus.end(), back_inserter(alignedConsensus));
    }
}



void LocalAssembly2::writeConsensus(const vector< pair<Base, uint64_t> >& consensus) const
{
    html <<
        "<h4>Consensus</h4>"
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << consensus.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace'>";

    for(uint64_t position=0; position<consensus.size(); position++) {
        const Base b = consensus[position].first;
        html << "<span title='" << position << "'>" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace'>";

    std::map<char, uint64_t> coverageLegend;

    for(const auto& p: consensus) {
        const uint64_t coverage = p.second;
        const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');
        coverageLegend.insert(make_pair(c, coverage));

        if(coverage < orientedReads.size()) {
            html << "<span style='background-color:Pink'>";
        }

        html << c;

        if(coverage < orientedReads.size()) {
            html << "</span>";
        }
    }

    html << "</table>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Symbol<th>Coverage";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";

}



void LocalAssembly2::writeConsensus() const
{
    html <<
        "<h3>Consensus</h3>"
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << consensus.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace'>";

    for(uint64_t position=0; position<consensus.size(); position++) {
        const Base b = consensus[position].first;
        html << "<span title='" << position << "'>" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace'>";

    std::map<char, uint64_t> coverageLegend;

    for(const auto& p: consensus) {
        const uint64_t coverage = p.second;
        const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');
        coverageLegend.insert(make_pair(c, coverage));

        if(coverage < orientedReads.size()) {
            html << "<span style='background-color:Pink'>";
        }

        html << c;

        if(coverage < orientedReads.size()) {
            html << "</span>";
        }
    }

    html << "</table>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Symbol<th>Coverage";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";

}



// Write the alignment for one step.
void LocalAssembly2::writeAlignment(
    const vector< vector<Base> >& inputSequences,
    const vector< pair<Base, uint64_t> >& consensus,
    const vector< vector<AlignedBase> >& alignment,
    const vector<AlignedBase>& alignedConsensus) const
{
    html <<
        "<h4>Alignment</h4>"
        "<table>"
        "<tr><th class=left>OrientedReadId"
        "<th class=left>Sequence<br>length"
        "<th class=left>Aligned sequence";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];

        html << "<tr><th>" << orientedReads[i].orientedReadId <<
            "<td class=centered>" << inputSequences[i].size() <<
            "<td style='font-family:monospace;white-space: nowrap'>";

        for(uint64_t j=0; j<alignmentRow.size(); j++) {
            const AlignedBase b = alignmentRow[j];
            const bool isMatch = (b == alignedConsensus[j]);

            if(not isMatch) {
                html << "<span style='background-color:Pink'>";
            }
            html << b;
            if(not isMatch) {
                html << "</span>";
            }
        }

    }

    html << "<tr><th>Consensus<td class=centered>" << consensus.size() <<
        "<td style='font-family:monospace;background-color:LightCyan;white-space:nowrap'>";

    uint64_t position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(not b.isGap()) {
            html << "<span title='" << position << "'>";
        }

        html << b;

        if(not b.isGap()) {
            html << "</span>";
            ++position;
        }
    }

    html << "<tr><th>Consensus coverage<td>"
        "<td style='font-family:monospace;white-space:nowrap'>";

    position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(b.isGap()) {
            html << "-";
        } else {
            const uint64_t coverage = consensus[position].second;
            const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');

            if(coverage < orientedReads.size()) {
                html << "<span style='background-color:Pink'>";
            }

            html << c;

            if(coverage < orientedReads.size()) {
                html << "</span>";
            }

            ++position;
        }
    }

    html << "</table>";

}



// Write the global alignment (all steps).
void LocalAssembly2::writeAlignment()
{
    html << "<h3>Alignment</h3>";

    // Get the input sequences for all reads.
    vector< vector<Base> > inputSequences(orientedReads.size());
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const OrientedRead& orientedRead = orientedReads[i];
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;
        for(uint32_t position=orientedRead.positionA(); position!=orientedRead.positionB(); position++) {
            inputSequences[i].push_back(anchors.reads.getOrientedReadBase(orientedReadId, position));
        }
    }

    html <<
        "<table>"
        "<tr><th class=left>OrientedReadId"
        "<th class=left>Sequence<br>length"
        "<th class=left>Aligned sequence";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];

        html << "<tr><th>" << orientedReads[i].orientedReadId <<
            "<td class=centered>" << inputSequences[i].size() <<
            "<td style='font-family:monospace;white-space: nowrap'>";

        for(uint64_t j=0; j<alignmentRow.size(); j++) {
            const AlignedBase b = alignmentRow[j];
            const bool isMatch = (b == alignedConsensus[j]);

            if(not isMatch) {
                html << "<span style='background-color:Pink'>";
            }
            html << b;
            if(not isMatch) {
                html << "</span>";
            }
        }

    }

    html << "<tr><th>Consensus<td class=centered>" << consensus.size() <<
        "<td style='font-family:monospace;background-color:LightCyan;white-space:nowrap'>";

    uint64_t position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(not b.isGap()) {
            html << "<span title='" << position << "'>";
        }

        html << b;

        if(not b.isGap()) {
            html << "</span>";
            ++position;
        }
    }

    html << "<tr><th>Consensus coverage<td>"
        "<td style='font-family:monospace;white-space:nowrap'>";

    position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(b.isGap()) {
            html << "-";
        } else {
            const uint64_t coverage = consensus[position].second;
            const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');

            if(coverage < orientedReads.size()) {
                html << "<span style='background-color:Pink'>";
            }

            html << c;

            if(coverage < orientedReads.size()) {
                html << "</span>";
            }

            ++position;
        }
    }

    html << "</table>";
}



void LocalAssembly2::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(const auto& p: consensus) {
        sequence.push_back(p.first);
    }
}
