// Shasta.
#include "Assembler.hpp"
#include "KmerChecker.hpp"
#include "Markers.hpp"
#include "MarkerKmers.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>



void Assembler::exploreReadMarkers(const vector<string>& request, ostream& html)
{

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);


    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << readsPointer->readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= readsPointer->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }


    // Access the read information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const auto sequence = readsPointer->getRead(readId);
    const span<const Marker> orientedReadMarkers = markers()[orientedReadId.getValue()];




    // Page title.
    html << "<h2 title='Markers of read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Markers of oriented read " << orientedReadId << "</h2>";

    // Write a table with some summary information for the markers of this oriented read.
    const double readMarkerDensity = double(orientedReadMarkers.size()) / double(sequence.baseCount);
    html <<
        "<table>"
        "<tr><th class=left>Length in bases<td class=centered>" << sequence.baseCount <<
        "<tr><th class=left>Number of markers<td class=centered>" << orientedReadMarkers.size() <<
        "<tr><th class=left>Average marker density for this read<td class=centered>" <<
        readMarkerDensity <<
        "</table>";


    // Count k-mers in this oriented read.
    // Reverse complemented k-mers are considered equivalent.
    const uint64_t k = assemblerInfo->k;
    std::map<Kmer, uint64_t> kmerFrequencyMap;
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const Kmer kmer = markers().getKmer(orientedReadId, uint32_t(ordinal));
        const Kmer rcKmer = kmer.reverseComplement(k);
        const Kmer& canonicalKmer = (kmer <= rcKmer) ? kmer : rcKmer;
        const auto it = kmerFrequencyMap.find(canonicalKmer);
        if(it == kmerFrequencyMap.end()) {
            kmerFrequencyMap.insert({canonicalKmer, 1});
        } else {
            ++(it->second);
        }
    }



    // Begin the main table containing one row for each marker.
    const uint64_t maxPeriod = 6;
    html <<
        "<p><table>"
        "<tr>"
        "<th>Marker<br>ordinal"
        "<th>Begin<br>position"
        "<th>End<br>position"
        "<th>Kmer"
        "<th>Frequency<br>in this<br>oriented read"
        ;
    if(markerKmers and markerKmers->isOpen()) {
        html << "<th>Global<br>frequency";
    }
    for(uint64_t period=1; period<=maxPeriod; period++) {
        html << "<th>Max copy<br>number<br>for repeats<br>of period " << period;
    }


    // Write one row for each marker.
    const auto& maxAnchorRepeatLength = httpServerData.options->maxAnchorRepeatLength;
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const uint64_t position = orientedReadMarkers[ordinal].position;
        const Kmer kmer = markers().getKmer(orientedReadId, uint32_t(ordinal));
        const Kmer rcKmer = kmer.reverseComplement(k);
        const Kmer canonicalKmer = (kmer <= rcKmer) ? kmer : rcKmer;

        html <<
            "<tr>"
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position <<
            "<td class=centered>" << position + k <<
            "<td class=centered style='font-family:Courier New;'>";
        kmer.write(html, k);

        // Frequency of this Kmer in this read.
        html << "<td class=centered>" << kmerFrequencyMap[canonicalKmer];

        // Global frequency of this Kmer.
        if(markerKmers and markerKmers->isOpen()) {
            html << "<td class=centered>" << markerKmers->getFrequency(kmer);
        }

        for(uint64_t period=1; period<=maxPeriod; period++) {
            const uint64_t copyNumber = kmer.countExactRepeatCopies(period, k);
            bool isHighCopyNumber = false;
            if(period - 1 < maxAnchorRepeatLength.size()) {
                const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[period - 1];
                isHighCopyNumber = (copyNumber > maxAllowedCopyNumber);

            }
            html << "<td class=centered";
            if(isHighCopyNumber) {
                html << " style='background-color:LightPink'";
            }
            html << ">" << copyNumber;
        }
    }

    html << "</table>";

}



void Assembler::exploreMarkerKmers(const vector<string>& request, ostream& html)
{
    SHASTA_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Marker k-mer</h2>";

    // Get the request parameters.
    string kmerString;
    getParameterValue(request, "kmer", kmerString);
    boost::trim(kmerString);

    // Write the form.
    html <<
        "<p><form><table>"
        "<tr><th class=left>K-mer"
        "<td><input type=text name=kmer style='font-family:monospace' "
        "size=" << k + 10 << " "
        "value='" << kmerString << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"
        "</table><input type=submit value='Get k-mer information'></form>";

    // If the k-mer string is empty, do nothing.
    if(kmerString.empty()) {
        return;
    }

    // Check the length.
    if(kmerString.size() != k) {
        html << "This k-mer is " << kmerString.size() << " bases long. "
            "This assembly uses " << k << "-base markers. "
            "Specify a " << k << "-mer." << endl;
        return;
    }

    // Construct the k-mer.
    Kmer kmer;
    for(uint64_t i=0; i<kmerString.size(); i++) {
        const char c = kmerString[i];
        const Base b = Base::fromCharacterNoException(c);
        if(not b.isValid()) {
            throw runtime_error("Invalid base character " + string(1, c) + " at k-mer position " + to_string(i));
        }
        kmer.set(i, b);
    }

    // Check if it is a marker.
    SHASTA_ASSERT(kmerChecker);
    if(not kmerChecker->isMarker(kmer)) {
        throw runtime_error("This assembly does not use this as a marker.");
    }



    // Summary table.
    const uint64_t coverage = markerKmers->getFrequency(kmer);
    const Kmer kmerRc = kmer.reverseComplement(k);
    html <<
        "<table><tr><th class=left>K-mer<td class=centered style='font-family:monospace'>";
    kmer.write(html, k);
    html <<
        "<tr><th class=left>Reverse complement K-mer<td class=centered style='font-family:monospace'>";
    kmerRc.write(html, k);
    html << "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "</table>";



    // Details table.
    vector<MarkerInfo> markerInfos;
    markerKmers->get(kmer, markerInfos);

    html <<
        "<p>"
        "<table>"
        "<tr><th>Oriented<br>read<th>Ordinal<th>Position";
    for(const MarkerInfo& markerInfo: markerInfos) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const uint32_t ordinal = markerInfo.ordinal;
        const Marker& marker = markers()[orientedReadId.getValue()][ordinal];
        const uint32_t position = marker.position;

        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position;
    }
    html << "</table>";


}
