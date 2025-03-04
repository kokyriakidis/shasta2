// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
using namespace shasta;



void Assembler::exploreReadRaw(
    const vector<string>& request,
    ostream& html)
{

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    string requestReadName;
    getParameterValue(request, "readName", requestReadName);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = getParameterValue(request, "endPosition", endPosition);


    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Read name"
        "<td><input type=text name=readName" <<
        (requestReadName.empty() ? "" : " value='" + requestReadName + "'") << ">"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";



    // Check that one and only one of readId and readName was entered.
    if(readIdIsPresent and not requestReadName.empty()) {
        html << "Specify either a numeric read id or a read name, but not both.";
        return;
    }
    if(not readIdIsPresent and requestReadName.empty()) {
        html << "Specify a numeric read id or a read name.";
        return;
    }

    // If a read name was specified, get the read id.
    if(not requestReadName.empty()) {
        readId = getReads().getReadId(requestReadName);
        if(readId == invalidReadId) {
            html << "A read with that name was not found. See ReadSummary.csv.";
            return;
        }
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }


    // Access the read information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const auto readName = reads->getReadName(readId);
    const span<const Marker> orientedReadMarkers = markers()[orientedReadId.getValue()];

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(getReads().getReadSequenceLength(readId));
    } else {
        endPosition++; // To include the base at `endPosition`.
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }



    // Page title.
    html << "<h2 title='Read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Oriented read " << orientedReadId << "</h2>";



    // Read information.
    html << "<table>";

    html << "<tr><th class=left>Read id<td>" << readId;

    html << "<tr><th class=left>Read name<td>";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(html));

    html << "<tr><th class=left>Length<td>" << getReads().getReadSequenceLength(readId);

    html << "<tr><th class=left>Length displayed<td>" << endPosition - beginPosition;

    html << "</table>";


    // Position scale labels.
    html << "<p>For precise alignment of the following section, use Firefox to display this page.\n";
    html << "<div style='font-family:Courier New;font-size:10pt;margin:0'>";
    for(size_t position=beginPosition; position<endPosition; ) {
        if((position%10)==0) {
            const string label = to_string(position);
            html << label;
            for(size_t i=0; i<10-label.size(); i++) {
                html << "&nbsp;";
            }
            position += 10;
        } else {
            html << "&nbsp;";
            ++position;
        }
    }
    html<< "<br>";

    // Position scale
    for(size_t position=beginPosition; position<endPosition; position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html<< "<br>";

    // Sequence.
    for(uint32_t position=beginPosition; position!=endPosition; position++) {
        html << getReads().getOrientedReadBase(orientedReadId, position);
    }
    html<< "<br>";



    // Display the markers.

    // If here are no markers, there is nothing to do.
    if(orientedReadMarkers.empty()) {
        html << "</div><p>This read has no markers.";
        return;
    }

    // Because markers can overlap, we have to display them on more than one row.
    // This vector will contain, for each row, the list of marker ordinals
    // to be displayed in this row.
    vector< vector<uint64_t> > markersByRow;

    const uint64_t k = assemblerInfo->k;
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const uint64_t position = orientedReadMarkers[ordinal].position;

        // If this marker begins before our beginPosition, it will not be displayed.
        if(position < beginPosition) {
            continue;
        }

        // If this marker ends after our endPosition, it will not be displayed.
        if(position + k > endPosition) {
            continue;
        }

        // Try all rows.
        for(uint64_t row=0; ; row++) {
            if(row >= markersByRow.size()) {
                markersByRow.resize(row + 1);
            }
            vector<uint64_t>& markersOnThisRow = markersByRow[row];
            if(markersOnThisRow.empty() or position > orientedReadMarkers[markersOnThisRow.back()].position + k) {
                markersOnThisRow.push_back(ordinal);
                break;
            }
        }
    }

    // Display the markers on each row.
    for(const vector<uint64_t>&markersOnThisRow: markersByRow) {

        // Loop over the markers on this row.
        uint64_t oldPosition = 0;
        for(const uint64_t ordinal: markersOnThisRow) {
            const Marker& marker = orientedReadMarkers[ordinal];
            const uint64_t position = marker.position - beginPosition;
            const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);

            // Write the required number of spaces.
            SHASTA_ASSERT((position==0) or (position > oldPosition));  // There must be at least a blank.
            for(uint64_t i=oldPosition; i<position; i++) {
                html << "&nbsp;";
            }
            oldPosition = position + k;

            // Write this marker as text.
            html << "<span title='Marker " << ordinal <<
                ", position " << marker.position <<
                "'>";
            for(uint64_t i=0; i<k; i++) {
                html << kmer[i];
            }
            html << "</span>";
        }
        html << "<br>";
    }

    html << "</div>";

}



void Assembler::exploreReadSequence(const vector<string>& request, ostream& html)
{

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = getParameterValue(request, "endPosition", endPosition);


    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


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
    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }


    // Access the read information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const auto readName = reads->getReadName(readId);
    const span<const Marker> orientedReadMarkers = markers()[orientedReadId.getValue()];

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(getReads().getReadSequenceLength(readId));
    } else {
        endPosition++; // To include the base at `endPosition`.
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }



    // Page title.
    html << "<h2 title='Read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Oriented read " << orientedReadId << "</h2>";



    // Read information.
    html << "<table>";

    html << "<tr><th class=left>Read id<td>" << readId;

    html << "<tr><th class=left>Read name<td>";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(html));

    html << "<tr><th class=left>Length<td>" << getReads().getReadSequenceLength(readId);

    html << "<tr><th class=left>Length displayed<td>" << endPosition - beginPosition;

    html << "</table>";


    // Position scale labels.
    html << "<p><div style='font-family:Courier New;font-size:10pt;margin:0'>";
    for(size_t position=beginPosition; position<endPosition; ) {
        if((position%10)==0) {
            const string label = to_string(position);
            html << label;
            for(size_t i=0; i<10-label.size(); i++) {
                html << "&nbsp;";
            }
            position += 10;
        } else {
            html << "&nbsp;";
            ++position;
        }
    }
    html<< "<br>";

    // Position scale
    html << "<nobr>";
    for(size_t position=beginPosition; position<endPosition; position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html << "</nobr>";
    html<< "<br>";

    // Sequence.
    for(uint32_t position=beginPosition; position!=endPosition; position++) {
        html << getReads().getOrientedReadBase(orientedReadId, position);
    }
    html<< "<br>";


    // Display the markers.

    // If here are no markers, there is nothing to do.
    if(orientedReadMarkers.empty()) {
        html << "</div><p>This read has no markers.";
        return;
    }

    // Because markers can overlap, we have to display them on more than one row.
    // This vector will contain, for each row, the list of marker ordinals
    // to be displayed in this row.
    vector< vector<uint64_t> > markersByRow;

    const uint64_t k = assemblerInfo->k;
    for(uint64_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const uint64_t position = orientedReadMarkers[ordinal].position;

        // If this marker begins before our beginPosition, it will not be displayed.
        if(position < beginPosition) {
            continue;
        }

        // If this marker ends after our endPosition, it will not be displayed.
        if(position + k > endPosition) {
            continue;
        }

        // Try all rows.
        for(uint64_t row=0; ; row++) {
            if(row >= markersByRow.size()) {
                markersByRow.resize(row + 1);
            }
            vector<uint64_t>& markersOnThisRow = markersByRow[row];
            if(markersOnThisRow.empty() or position > orientedReadMarkers[markersOnThisRow.back()].position + k) {
                markersOnThisRow.push_back(ordinal);
                break;
            }
        }
    }

    // Display the markers on each row.
    for(const vector<uint64_t>&markersOnThisRow: markersByRow) {

        // Loop over the markers on this row.
        uint64_t oldPosition = 0;
        for(const uint64_t ordinal: markersOnThisRow) {
            const Marker& marker = orientedReadMarkers[ordinal];
            const uint64_t position = marker.position - beginPosition;
            const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);

            // Write the required number of spaces.
            SHASTA_ASSERT((position==0) or (position > oldPosition));  // There must be at least a blank.
            for(uint64_t i=oldPosition; i<position; i++) {
                html << "&nbsp;";
            }
            oldPosition = position + k;




            // For assembly mode 3 we write the first half of the marker
            // a different color than the second half.
            // This helps visualize the anchors.
            SHASTA_ASSERT((k % 2) == 0);
            html << "<span title='Marker " << ordinal <<
                ", position " << marker.position <<
                "'>";

            // The first half.
            html << "<span style='color:Green'>";
            for(uint64_t i=0; i<k/2; i++) {
                html << kmer[i];
            }
            html << "</span>";

            // The second half.
            html << "<span style='color:Red'>";
            for(uint64_t i=k/2; i<k; i++) {
                html << kmer[i];
            }
            html << "</span>";

            html << "</span>";
        }
        html << "<br>";
    }

    html << "</div>";

}



void Assembler::exploreLookupRead(const vector<string>& request, ostream& html)
{

    string requestReadName;
    getParameterValue(request, "readName", requestReadName);

    // Write the form.
    html <<
        "<form>"
        "<table>"
        "<tr>"
        "<th class=left>Read name"
        "<td><input type=text name=readName" <<
        (requestReadName.empty() ? "" : " value='" + requestReadName + "'") << ">"
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not requestReadName.empty()) {
        const ReadId readId = getReads().getReadId(requestReadName);
        if(readId == invalidReadId) {
            html << "A read with that name was not found. See ReadSummary.csv.";
            return;
        }

        html << "<br>Read id for this assembly is " << readId;
    }
}



