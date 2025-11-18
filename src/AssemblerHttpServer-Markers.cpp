// Shasta.
#include "Assembler.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "graphvizToHtml.hpp"
#include "KmerChecker.hpp"
#include "Markers.hpp"
#include "MarkerKmerPair.hpp"
#include "MarkerKmers.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>
#include <queue>



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



void Assembler::exploreMarkerKmer(const vector<string>& request, ostream& html)
{
    SHASTA2_ASSERT(markerKmers and markerKmers->isOpen());

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
        "size=" << k << " "
        "value='" << kmerString << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"
        "</table><br><input type=submit value='Get k-mer information'></form>";

    // If the k-mer string is empty, do nothing.
    if(kmerString.empty()) {
        return;
    }

    // Check the length.
    if(kmerString.size() != k) {
        html << "<br>This k-mer is " << kmerString.size() << " bases long. "
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
    SHASTA2_ASSERT(kmerChecker);
    if(not kmerChecker->isMarker(kmer)) {
        throw runtime_error("This assembly does not use this as a marker.");
    }



    // Summary table.
    const uint64_t coverage = markerKmers->getFrequency(kmer);
    const Kmer kmerRc = kmer.reverseComplement(k);
    html <<
        "<br><table><tr><th class=left>K-mer<td class=centered style='font-family:monospace'>";
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
        "<br>"
        "<table>"
        "<tr><th>Oriented<br>read<th>Ordinal<th>Position<th>Repeated<br>ReadId";
    for(uint64_t i=0; i<markerInfos.size(); i++) {
        const MarkerInfo& markerInfo = markerInfos[i];
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const uint32_t ordinal = markerInfo.ordinal;
        const Marker& marker = markers()[orientedReadId.getValue()][ordinal];
        const uint32_t position = marker.position;

        // Figure out if it is a repeated ReadId.
        bool isRepeatedReadId = false;
        if(i != 0) {
            isRepeatedReadId = (readId == markerInfos[i-1].orientedReadId.getReadId());
        }
        if(i != markerInfos.size() - 1) {
            isRepeatedReadId = isRepeatedReadId or
                (readId == markerInfos[i+1].orientedReadId.getReadId());
        }

        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position <<
            "<td class=centered>";
        if(isRepeatedReadId) {
            html << "&#10003;";
        }
    }
    html << "</table>";


}


void Assembler::exploreMarkerKmerAnalysisWithMarkerOffset(const vector<string>& request, ostream& html)
{
    SHASTA2_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Marker k-mer</h2>";

    // Get the request parameters.
    string kmerString;
    getParameterValue(request, "kmer", kmerString);
    boost::trim(kmerString);

    int64_t ordinalOffset = 3;
    getParameterValue(request, "ordinalOffset", ordinalOffset);



    // Write the form.
    html <<
        "<p><form><table>"

        "<tr><th class=left>K-mer"
        "<td><input type=text name=kmer style='font-family:monospace' "
        "size=" << k << " "
        "value='" << kmerString << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"

        "<tr><th class=left>Ordinal offset"
        "<td><input type=text name=ordinalOffset size=6 value=" << ordinalOffset << " style='text-align:center'>"

        "</table><br><input type=submit value='Get k-mer information'></form>";



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
    SHASTA2_ASSERT(kmerChecker);
    if(not kmerChecker->isMarker(kmer)) {
        throw runtime_error("This assembly does not use this as a marker.");
    }



    // Summary table.
    const uint64_t coverage = markerKmers->getFrequency(kmer);
    const Kmer kmerRc = kmer.reverseComplement(k);
    html <<
        "<br><table><tr><th class=left>K-mer<td class=centered style='font-family:monospace'>";
    kmer.write(html, k);
    html <<
        "<tr><th class=left>Reverse complement K-mer<td class=centered style='font-family:monospace'>";
    kmerRc.write(html, k);
    html << "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "</table>";

    // Get the MarkerInfos for this K-mer.
    vector<MarkerInfo> markerInfos;
    markerKmers->get(kmer, markerInfos);



    // For each of these MarkerInfos, gather the k-mers by moving
    // forward or backward by up to ordinalOffset markers.
    class KmerInfo {
    public:
        bool isAvailable = false;
        Kmer kmer;
        uint64_t id = invalid<uint64_t>;    // Index in the kmers vector below.
        KmerInfo() {}
        KmerInfo(const Kmer& kmer) : isAvailable(true), kmer(kmer) {}
    };
    vector< vector<KmerInfo> > kmerInfos(markerInfos.size());
    vector<Kmer> kmers;
    for(uint64_t i=0; i<markerInfos.size(); i++) {
        const MarkerInfo& markerInfo = markerInfos[i];
        kmerInfos[i].reserve(2 * ordinalOffset + 1);
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto orientedReadMarkers = markers()[orientedReadId.getValue()];

        // Use signed arithmethic to simplify checks.
        const int64_t markerCount = int64_t(orientedReadMarkers.size());
        const int64_t ordinal0 = markerInfo.ordinal;
        for(int64_t offset=-ordinalOffset; offset<=ordinalOffset; offset++) {
            const int64_t ordinal1 = ordinal0 + offset;
            if((ordinal1 < 0) or (ordinal1 >= markerCount)) {
                kmerInfos[i].push_back(KmerInfo());
            } else {
                const Kmer kmer = markers().getKmer(orientedReadId, uint32_t(ordinal1));
                kmerInfos[i].push_back(KmerInfo(kmer));
                kmers.push_back(kmer);
            }
        }
        SHASTA2_ASSERT(kmerInfos[i].size() == 2 * uint64_t(ordinalOffset) + 1);
    }
    vector<uint64_t> count;
    deduplicateAndCount(kmers, count);

    // Locate our start k-mer in the table.
    const auto it = lower_bound(kmers.begin(), kmers.end(), kmer);
    SHASTA2_ASSERT(it != kmers.end());
    SHASTA2_ASSERT(*it == kmer);
    const uint64_t iStart = it - kmers.begin();

    html <<
        "<br>Found " << kmers.size() << " distinct k-mers."
        "<br><table><tr><th>Id<th>Number of<br>occurrences<th>K-mer";
    for(uint64_t i=0; i<kmers.size(); i++) {
        html << "<tr";
        if(i == iStart) {
             html << " style='background-color:Pink'";
        }
        html <<
            "><td class=centered>" << i <<
            "<td class=centered>" << count[i] <<
            "<td class=centered style='font-family:monospace'>";
        kmers[i].write(html, assemblerInfo->k);
    }
    html << "</table>";


    html <<
        "<br>"
        "<table>"
        "<tr><th>Oriented<br>read<th>Ordinal<th>Position<th>Repeated<br>ReadId";
    for(int64_t offset=-ordinalOffset; offset<=ordinalOffset; offset++) {
        html << "<th>" << offset;
    }
    for(uint64_t i=0; i<markerInfos.size(); i++) {
        const MarkerInfo& markerInfo = markerInfos[i];
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const uint32_t ordinal = markerInfo.ordinal;
        const Marker& marker = markers()[orientedReadId.getValue()][ordinal];
        const uint32_t position = marker.position;

        // Figure out if it is a repeated ReadId.
        bool isRepeatedReadId = false;
        if(i != 0) {
            isRepeatedReadId = (readId == markerInfos[i-1].orientedReadId.getReadId());
        }
        if(i != markerInfos.size() - 1) {
            isRepeatedReadId = isRepeatedReadId or
                (readId == markerInfos[i+1].orientedReadId.getReadId());
        }

        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position <<
            "<td class=centered>";
        if(isRepeatedReadId) {
            html << "&#10003;";
        }
        for(int64_t offset=-ordinalOffset; offset<=ordinalOffset; offset++) {
            KmerInfo& kmerInfo = kmerInfos[i][offset + ordinalOffset];
            html << "<td class=centered>";
            if(kmerInfo.isAvailable) {
                auto it = std::lower_bound(kmers.begin(), kmers.end(), kmerInfo.kmer);
                kmerInfo.id = it - kmers.begin();
                html << kmerInfo.id;
            }
        }
    }
    html << "</table>";



    // Create a graph with one vertex for each of these k-mers.
    // Generate edges by following the reads.
    // Each vertex and edge stores coverage.
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        uint64_t,
        uint64_t>;
    Graph graph(kmers.size());
    for(uint64_t i=0; i<kmers.size(); i++) {
        graph[i] = count[i];
    }
    for(uint64_t i=0; i<kmerInfos.size(); i++) {
        for(uint64_t j1=1; j1<kmerInfos[i].size(); j1++) {
            const uint64_t j0 = j1 - 1;
            const KmerInfo& kmerInfo0 = kmerInfos[i][j0];
            const KmerInfo& kmerInfo1 = kmerInfos[i][j1];
            if(not (kmerInfo0.isAvailable and kmerInfo1.isAvailable)) {
                continue;
            }
            const uint64_t id0 = kmerInfo0.id;
            const uint64_t id1 = kmerInfo1.id;
            Graph::edge_descriptor e;
            bool edgeWasFound = false;
            tie(e, edgeWasFound) = edge(id0, id1, graph);
            if(edgeWasFound) {
                ++graph[e];
            } else {
                tie(e, edgeWasFound) = add_edge(id0, id1, 1, graph);
            }
        }
    }



    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dot(dotFileName);
    dot << "digraph KmerGraph {\n";
    BGL_FORALL_VERTICES(v, graph, Graph) {
        dot << v;
        dot << " [";
        dot << "label=\"" << v << "\\n" << graph[v] << "\"";
        if(v == iStart) {
            dot << " style=filled fillcolor=Pink";
        }
        dot << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);
        dot << v0 << "->" << v1 <<
            " [label=\"" << graph[e] << "\""
            " penwidth=" << std::setprecision(2) << 0.2 * double(graph[e]) <<
            "];\n";
    }
    dot << "}\n";
    dot.close();

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle";
    html << "<p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);

}



void Assembler::exploreMarkerKmerAnalysisWithBaseOffset(const vector<string>& request, ostream& html)
{
    SHASTA2_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Marker k-mer</h2>";

    // Get the request parameters.
    string kmerString;
    getParameterValue(request, "kmer", kmerString);
    boost::trim(kmerString);

    int64_t baseOffset = 1000;
    getParameterValue(request, "baseOffset", baseOffset);



    // Write the form.
    html <<
        "<p><form><table>"

        "<tr><th class=left>K-mer"
        "<td><input type=text name=kmer style='font-family:monospace' "
        "size=" << k << " "
        "value='" << kmerString << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"

        "<tr><th class=left>Base offset"
        "<td><input type=text name=baseOffset size=6 value=" << baseOffset << " style='text-align:center'>"

        "</table><br><input type=submit value='Get k-mer information'></form>";



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
    Kmer kmer0;
    for(uint64_t i=0; i<kmerString.size(); i++) {
        const char c = kmerString[i];
        const Base b = Base::fromCharacterNoException(c);
        if(not b.isValid()) {
            throw runtime_error("Invalid base character " + string(1, c) + " at k-mer position " + to_string(i));
        }
        kmer0.set(i, b);
    }

    // Check if it is a marker.
    SHASTA2_ASSERT(kmerChecker);
    if(not kmerChecker->isMarker(kmer0)) {
        throw runtime_error("This assembly does not use this as a marker.");
    }



    // Summary table.
    const uint64_t coverage = markerKmers->getFrequency(kmer0);
    const Kmer kmer0Rc = kmer0.reverseComplement(k);
    html <<
        "<br><table><tr><th class=left>K-mer<td class=centered style='font-family:monospace'>";
    kmer0.write(html, k);
    html <<
        "<tr><th class=left>Reverse complement K-mer<td class=centered style='font-family:monospace'>";
    kmer0Rc.write(html, k);
    html << "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "</table>";

    // Get the MarkerInfos for this K-mer.
    vector<MarkerInfo> markerInfos;
    markerKmers->get(kmer0, markerInfos);



    // For each of these MarkerInfos, gather k-mers by moving
    // forward or backward by up to baseOffset bases.
    class KmerInfo {
    public:
        int32_t ordinalOffset;
        int32_t positionOffset;
        Kmer kmer;
        uint64_t id = invalid<uint64_t>;    // Index in the kmers vector below.
        KmerInfo(
            int32_t ordinalOffset,
            int32_t positionOffset,
            const Kmer& kmer) :
            ordinalOffset(ordinalOffset),
            positionOffset(positionOffset),
            kmer(kmer)
            {}
    };

    vector< vector<KmerInfo> > kmerInfos(markerInfos.size());
    vector<Kmer> kmers;

    vector<KmerInfo> forwardKmerInfos;
    vector<KmerInfo> backwardKmerInfos;
    for(uint64_t i=0; i<markerInfos.size(); i++) {
        const MarkerInfo& markerInfo = markerInfos[i];
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto orientedReadMarkers = markers()[orientedReadId.getValue()];
        const int32_t markerCount = int32_t(orientedReadMarkers.size());
        const int32_t ordinal0 = int32_t(markerInfo.ordinal);
        const int32_t position0 = int32_t(orientedReadMarkers[ordinal0].position);

        // Move forward.
        forwardKmerInfos.clear();
        for(int32_t ordinalOffset=0; ; ordinalOffset++) {
            const int64_t ordinal1 = ordinal0 + ordinalOffset;
            if(ordinal1 >= markerCount) {
                break;
            }
            const int32_t position1 = int32_t(orientedReadMarkers[ordinal1].position);
            const int32_t positionOffset = position1 - position0;
            if(positionOffset > baseOffset) {
                break;
            }
            const Kmer kmer1 = markers().getKmer(orientedReadId, uint32_t(ordinal1));
            forwardKmerInfos.push_back(KmerInfo(ordinalOffset, positionOffset, kmer1));
            kmers.push_back(kmer1);
        }

        // Move backward.
        backwardKmerInfos.clear();
        for(int32_t ordinalOffset=1; ; ordinalOffset++) {
            const int64_t ordinal1 = ordinal0 - ordinalOffset;
            if(ordinal1 < 0) {
                break;
            }
            const int32_t position1 = int32_t(orientedReadMarkers[ordinal1].position);
            const int32_t positionOffset = position0 - position1;
            if(positionOffset > baseOffset) {
                break;
            }
            const Kmer kmer1 = markers().getKmer(orientedReadId, uint32_t(ordinal1));
            backwardKmerInfos.push_back(KmerInfo(ordinalOffset, positionOffset, kmer1));
            kmers.push_back(kmer1);
        }

        copy(backwardKmerInfos.rbegin(), backwardKmerInfos.rend(), back_inserter(kmerInfos[i]));
        copy(forwardKmerInfos.begin(), forwardKmerInfos.end(), back_inserter(kmerInfos[i]));
    }

    vector<uint64_t> count;
    deduplicateAndCount(kmers, count);



    // Fill in the id fields in all the KmerInfos.
    for(vector<KmerInfo>& v: kmerInfos) {
        for(KmerInfo& kmerInfo: v) {
            const auto it = std::lower_bound(kmers.begin(), kmers.end(), kmerInfo.kmer);
            SHASTA2_ASSERT(it != kmers.end());
            SHASTA2_ASSERT(*it == kmerInfo.kmer);
            kmerInfo.id = it - kmers.begin();
        }
    }



    // Locate our start k-mer in the table.
    const auto it0 = lower_bound(kmers.begin(), kmers.end(), kmer0);
    SHASTA2_ASSERT(it0 != kmers.end());
    SHASTA2_ASSERT(*it0 == kmer0);
    const uint64_t i0 = it0 - kmers.begin();

    // Write the k-mer table.
    html <<
        "<br>Found " << kmers.size() << " distinct k-mers."
        "<br><br><table><tr><th>Id<th>Number of<br>occurrences<th>K-mer";
    for(uint64_t i=0; i<kmers.size(); i++) {
        html << "<tr";
        if(i == i0) {
             html << " style='background-color:Pink'";
        }
        html <<
            "><td class=centered>" << i <<
            "<td class=centered>" << count[i] <<
            "<td class=centered style='font-family:monospace'>";
        kmers[i].write(html, assemblerInfo->k);
    }
    html << "</table>";



    // Create a graph with one vertex for each of these k-mers.
    // Generate edges by following the reads.
    // Each vertex and edge stores coverage.
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        uint64_t,
        uint64_t>;
    Graph graph(kmers.size());
    for(uint64_t i=0; i<kmers.size(); i++) {
        graph[i] = count[i];
    }
    for(uint64_t i=0; i<kmerInfos.size(); i++) {
        for(uint64_t j1=1; j1<kmerInfos[i].size(); j1++) {
            const uint64_t j0 = j1 - 1;
            const KmerInfo& kmerInfo0 = kmerInfos[i][j0];
            const KmerInfo& kmerInfo1 = kmerInfos[i][j1];
            const uint64_t id0 = kmerInfo0.id;
            const uint64_t id1 = kmerInfo1.id;
            Graph::edge_descriptor e;
            bool edgeWasFound = false;
            tie(e, edgeWasFound) = edge(id0, id1, graph);
            if(edgeWasFound) {
                ++graph[e];
            } else {
                tie(e, edgeWasFound) = add_edge(id0, id1, 1, graph);
            }
        }
    }



    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dot(dotFileName);
    dot << "digraph KmerGraph {\n";
    BGL_FORALL_VERTICES(v, graph, Graph) {
        dot << v;
        dot << " [";
        dot << "label=\"" << v << "\\n" << graph[v] << "\"";
        if(v == i0) {
            dot << " style=filled fillcolor=Pink";
        }
        dot << " tooltip=\"";
        kmers[v].write(dot, k);
        dot << "\"";
        dot << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);
        dot << v0 << "->" << v1 <<
            " [label=\"" << graph[e] << "\""
            " penwidth=" << std::setprecision(2) << 0.2 * double(graph[e]) <<
            "];\n";
    }
    dot << "}\n";
    dot.close();

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle";
    html << "<p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);

}



void Assembler::exploreLocalKmerGraph(const vector<string>& request, ostream& html)
{
    SHASTA2_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Local k-mer graph</h2>";

    // Get the request parameters.
    string kmersString = getKeywordValue(request, "kmer");
    boost::trim(kmersString);

    uint64_t distance = 10;
    getParameterValue(request, "distance", distance);

    uint64_t minVertexCoverage = 2;
    getParameterValue(request, "minVertexCoverage", minVertexCoverage);

    uint64_t minEdgeCoverage = 2;
    getParameterValue(request, "minEdgeCoverage", minEdgeCoverage);



    // Write the form.
    html <<
        "<p><form><table>"

        "<tr><th class=left>Starting k-mers"
        "<td><textarea name=kmer style='font-family:monospace' "
        "rows=\"4\" cols=\"" << k << "\">" <<
        kmersString << "</textarea>"

        "<tr><th class=left>Distance"
        "<td class=centered><input type=text name=distance size=6 value=" << distance << " style='text-align:center'>"

        "<tr><th class=left>Minimum vertex coverage"
        "<td class=centered><input type=text name=minVertexCoverage size=6 value=" << minVertexCoverage << " style='text-align:center'>"

        "<tr><th class=left>Minimum edge coverage"
        "<td class=centered><input type=text name=minEdgeCoverage size=6 value=" << minEdgeCoverage << " style='text-align:center'>"

        "</table><br><input type=submit value='Get k-mer information'></form>";

    // If the k-mer string is empty, do nothing.
    if(kmersString.empty()) {
        return;
    }

    // Parse the kmersString to get the starting Kmers.
    std::istringstream s(kmersString);
    vector<Kmer> startingKmers;
    while(s) {
        string kmerString;
        getline(s, kmerString);
        boost::trim(kmerString);
        if(kmerString.empty()) {
            continue;
        }

        // Check the length.
        if(kmerString.size() != k) {
            html << "<br>This k-mer is " << kmerString.size() << " bases long. "
                "This assembly uses " << k << "-base markers. "
                "Specify a " << k << "-mer.<br>" << kmerString;
            html << "<br>";
            return;
        }

        // Construct the k-mer.
        Kmer kmer;
        for(uint64_t i=0; i<kmerString.size(); i++) {
            const char c = kmerString[i];
            const Base b = Base::fromCharacterNoException(c);
            if(not b.isValid()) {
                html << "<br>This k-mer contains an invalid base character at position " << i <<
                    "<br>" << kmerString;
                return;
            }
            kmer.set(i, b);
        }

        // Check if it is a marker.
        SHASTA2_ASSERT(kmerChecker);
        if(not kmerChecker->isMarker(kmer)) {
            html << "<br>This assembly does not use this k-mer as a marker:"
                "<br>" << kmerString;
            return;
        }

        startingKmers.push_back(kmer);
    }



    // The local k-mer graph.
    class Vertex {
    public:
        uint64_t id;
        Kmer kmer;
        uint64_t coverage;
        uint64_t distance;
        Vertex(
            uint64_t id,
            Kmer kmer,
            uint64_t coverage,
            uint64_t distance) :
            id(id), kmer(kmer), coverage(coverage), distance(distance)
            {}
    };
    class Edge {
    public:
        uint64_t coverage;
        Edge(uint64_t coverage) : coverage(coverage) {}
    };
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge>;
    using vertex_descriptor = Graph::vertex_descriptor;
    // using edge_descriptor = Graph::edge_descriptor;
    Graph graph;

    // Construct graph vertices (the k-mers).
    uint64_t nextVertexId = 0;
    std::map<Kmer, vertex_descriptor> vertexMap;
    std::queue<vertex_descriptor> q;
    for(const Kmer& kmer: startingKmers) {
        const uint64_t coverage = markerKmers->getFrequency(kmer);
        const vertex_descriptor v = add_vertex(Vertex(nextVertexId++, kmer, coverage, 0), graph);
        vertexMap.insert(make_pair(kmer, v));
        q.push(v);
    }

    vector<pair<Kmer, uint64_t> > nextKmers;
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        const Vertex& vertex0 = graph[v0];
        const Kmer& kmer0 = vertex0.kmer;
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        markerKmers->getNext(kmer0, nextKmers);
        for(const auto& p: nextKmers) {
            if(p.second < minEdgeCoverage) {
                continue;
            }
            const Kmer& kmer1 = p.first;

            const uint64_t coverage1 = markerKmers->getFrequency(kmer1);
            if(coverage1 >= minVertexCoverage) {
                auto it1 = vertexMap.find(kmer1);
                if(it1 == vertexMap.end()) {
                    const vertex_descriptor v1 = add_vertex(Vertex(nextVertexId++, kmer1, coverage1, distance1), graph);
                    vertexMap.insert(make_pair(kmer1, v1));
                    if(distance1 < distance) {
                        q.push(v1);
                    }
                }
            }
        }

        markerKmers->getPrevious(kmer0, nextKmers);
        for(const auto& p: nextKmers) {
            if(p.second < minEdgeCoverage) {
                continue;
            }
            const Kmer& kmer1 = p.first;

            const uint64_t coverage1 = markerKmers->getFrequency(kmer1);
            if(coverage1 >= minVertexCoverage) {
                auto it1 = vertexMap.find(kmer1);
                if(it1 == vertexMap.end()) {
                    const vertex_descriptor v1 = add_vertex(Vertex(nextVertexId++, kmer1, coverage1, distance1), graph);
                    vertexMap.insert(make_pair(kmer1, v1));
                    if(distance1 < distance) {
                        q.push(v1);
                    }
                }
            }
        }
    }



    // Now add edges.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        const Kmer& kmer0 = graph[v0].kmer;
        markerKmers->getNext(kmer0, nextKmers);
        for(const auto& p: nextKmers) {
            const Kmer& kmer1 = p.first;
            const uint64_t coverage1 = p.second;
            if(coverage1 >= minEdgeCoverage) {
                const auto it1 = vertexMap.find(kmer1);
                if(it1 != vertexMap.end()) {
                    const vertex_descriptor v1 = it1->second;
                    add_edge(v0, v1, Edge(coverage1), graph);
                }
            }
        }
    }
    html << "<br>The local k-mer graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges.";



    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dot(dotFileName);
    dot << "digraph KmerGraph {\n";
    dot << "node [fontname = \"Courier New\" fontsize=10]\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        dot << vertex.id;
        dot << " [";

        dot << "label=\"" << /* vertex.id << "\\n" << */ graph[v].coverage;
        for(uint64_t i=0; i<k; i++) {
            if((i % 20) == 0) {
                dot << "\\n";
            }
            dot << vertex.kmer[i];
        }
        bool hasAnchor = false;
        if(anchorsPointer) {
            const AnchorId anchorId = anchors().getAnchorIdFromKmer(vertex.kmer);
            if(anchorId != invalid<AnchorId>) {
                dot << "\\n" << anchorIdToString(anchorId);
                hasAnchor = true;
            }
        }
        dot << "\"";

        if(vertex.distance == 0) {
            dot << " style=filled fillcolor=Pink";
        } else if(hasAnchor) {
            dot << " style=filled fillcolor=LightCyan";
        }

        dot << " tooltip=\"";
        graph[v].kmer.write(dot, k);
        dot << "\"";
        dot << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);
        const uint64_t coverage = graph[e].coverage;
        dot << graph[v0].id << "->" << graph[v1].id <<
            " [label=\"" << coverage << "\""
            " penwidth=" << std::setprecision(2) << 0.2 * double(coverage) <<
            "];\n";
    }
    dot << "}\n";
    dot.close();

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle";
    html << "<br><br>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);
}



void Assembler::exploreMarkerKmerPair(
    const vector<string>& request,
    ostream& html)
{
    SHASTA2_ASSERT(markerKmers and markerKmers->isOpen());

    const uint64_t k = assemblerInfo->k;

    html << "<h2>Marker k-mer pair</h2>";

    // Get the request parameters.
    string kmerString0;
    getParameterValue(request, "kmer0", kmerString0);
    boost::trim(kmerString0);

    string kmerString1;
    getParameterValue(request, "kmer1", kmerString1);
    boost::trim(kmerString1);



    // Write the form.
    html <<
        "<p><form><table>"

        "<tr><th class=left>Left k-mer (kmer0)"
        "<td><input type=text name=kmer0 style='font-family:monospace' "
        "size=" << k << " "
        "value='" << kmerString0 << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"

        "<tr><th class=left>Right k-mer (kmer1)"
        "<td><input type=text name=kmer1 style='font-family:monospace' "
        "size=" << k << " "
        "value='" << kmerString1 << "'" <<
        " title='Enter a " << k << "-base k-mer.'>"

        "</table><input type=submit value='Get k-mer pair information'></form>";



    // If the k-mer strings are empty, do nothing.
    if(kmerString0.empty() or kmerString1.empty()) {
        return;
    }

    // Check the length.
    if(kmerString0.size() != k) {
        html << "<p>The left k-mer is " << kmerString0.size() << " bases long. "
            "This assembly uses " << k << "-base markers. "
            "Specify a " << k << "-mer." << endl;
        return;
    }
    if(kmerString1.size() != k) {
        html << "<p>The right k-mer is " << kmerString1.size() << " bases long. "
            "This assembly uses " << k << "-base markers. "
            "Specify a " << k << "-mer." << endl;
        return;
    }

    // Construct the k-mers.
    Kmer kmer0;
    for(uint64_t i=0; i<kmerString0.size(); i++) {
        const char c = kmerString0[i];
        const Base b = Base::fromCharacterNoException(c);
        if(not b.isValid()) {
            throw runtime_error("Invalid base character " + string(1, c) + " at left k-mer position " + to_string(i));
        }
        kmer0.set(i, b);
    }
    Kmer kmer1;
    for(uint64_t i=0; i<kmerString1.size(); i++) {
        const char c = kmerString1[i];
        const Base b = Base::fromCharacterNoException(c);
        if(not b.isValid()) {
            throw runtime_error("Invalid base character " + string(1, c) + " at right k-mer position " + to_string(i));
        }
        kmer1.set(i, b);
    }

    // Check if they are markers.
    SHASTA2_ASSERT(kmerChecker);
    if(not kmerChecker->isMarker(kmer0)) {
        html << "<p>The left k-mer is not a marker.";
        return;
    }
    if(not kmerChecker->isMarker(kmer1)) {
        html << "<p>The right k-mer is not a marker.";
        return;
    }

    // Create the MarkerKmerPair and write it out.
    const uint32_t maxPositionOffset = 5000;
    const MarkerKmerPair markerKmerPair(*markerKmers, kmer0, kmer1, maxPositionOffset);
    markerKmerPair.writeSummary(html, k);
    markerKmerPair.writeSequences(html);
    markerKmerPair.writeCommonOrientedReads(html);
    markerKmerPair.writeAlignment(html);
    markerKmerPair.writePairAlignmentDistances(html);

}



void Assembler::exploreFindMarkerKmers(const vector<string>& request, ostream& html)
{
	html << "<h2>Find markers in sequence</h2>";

    // Get the request parameters.
    string sequenceString;
    getParameterValue(request, "sequence", sequenceString);
    boost::trim(sequenceString);

    // Write the form.
    html <<
        "<p><form><table>"
        "<tr><th class=left>Sequence"
        "<td><textarea name=sequence style='font-family:monospace' "
        "rows=\"4\" cols=\"" << 100 << "\">" <<
		sequenceString << "</textarea>"
        "</table><br><input type=submit value='Get marker k-mers in this sequence'></form>";

    // If the sequence is empty, do nothing.
    if(sequenceString.empty()) {
    	return;
    }

    // Convert it to a vector of Bases.
    vector<Base> sequence;
    for(const char c: sequenceString) {
    	const Base b = Base::fromCharacterNoException(c);
    	if(not b.isValid()) {
    		html << "Invalid base character " << c << endl;
    		return;
    	}
    	sequence.push_back(b);
    }

    // Check if the sequence is shorter than a marker k-mer.
    const uint64_t k = assemblerInfo->k;
    if(sequence.size() < k) {
    	html << "<br>Sequence is " << sequence.size() << " bases long, which is shorter than marker length " << k << ".";
    	return;
    }



    // Find marker k-mers in this sequence.
    // Loop over starting positions of the k-mer.
    html << "<h3>Marker k-mers in this sequence</h3>"
    	"<table><tr><th>Begin<th>End<th>K-mer";
    uint64_t markerKmerCount = 0;
    for(uint64_t begin=0; ; ++begin) {
    	const uint64_t end = begin + k;
    	if(end > sequence.size()) {
    		break;
    	}

    	// Gather the k-mer starting here.
    	Kmer kmer;
    	for(uint64_t i=0; i<k; i++) {
    		kmer.set(i, sequence[begin + i]);
    	}

    	if(not kmerChecker->isMarker(kmer)) {
    		continue;
    	}
    	++markerKmerCount;

    	html <<
    		"<tr>"
        	"<td class=centered>" << begin <<
	    	"<td class=centered>" << end <<
			"<td style='font-family:monospace'>"
			"<a href='exploreMarkerKmer?kmer=";
    	kmer.write(html, k);
    	html << "'>";
    	kmer.write(html, k);
    	html << "</a>";
    }
    html << "</table>";
    html << "<br>Found " << markerKmerCount << " markers in " <<
    	sequence.size() << " bases of sequence."
		"<br>Marker density is " << double(markerKmerCount) / double(sequence.size() - k + 1);
}
