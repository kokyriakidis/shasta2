#include "ExternalAnchors.hpp"
#include "invalid.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta;



ExternalAnchors::ExternalAnchors(const string& name)
{
    data.createNew(name, pageSize);
    names.createNew(name + "-Names", pageSize);
}



void ExternalAnchors::beginNewAnchor(const string& anchorName)
{
    data.appendVector();

    names.appendVector();
    for(const char c: anchorName) {
        names.append(c);
    }
}



void ExternalAnchors::addOrientedRead(
    ReadId readId,
    Strand strand,
    uint32_t position)
{
    data.append(OrientedRead(readId, strand, position));
}



ExternalAnchors::OrientedRead::OrientedRead()
{}



ExternalAnchors::OrientedRead::OrientedRead(
    ReadId readId,
    Strand strand,
    uint32_t position) :
    orientedReadId(readId, strand), position(position)
{}



// This is used to access an existing ExternalAnchors.
ExternalAnchors::ExternalAnchors(const string& name, const AccessExisting&)
{
    data.accessExistingReadOnly(name);
    names.accessExistingReadOnly(name + "-Names");
}




// Write information about the i-th external anchor.
void ExternalAnchors::write(
    ostream& s,
    uint64_t i,
    uint64_t k,
    const Reads& reads,
    const Markers& markers) const
{
    const span<const OrientedRead>& orientedReads = data[i];
    const span<const char> name = names[i];

    s << "External anchor at position " << i << " in external anchor files." << endl;
    s << "External anchor name ";
    copy(name.begin(), name.end(), ostream_iterator<char>(s));
    s << endl;

    for(const OrientedRead& orientedRead: orientedReads) {
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;
        const uint32_t position = orientedRead.position;
        const span<const Marker> orientedReadMarkers = markers[orientedReadId.getValue()];

        // Find the ordinal at this position, if there is one.
        uint32_t ordinal = invalid<uint32_t>;
        Marker targetMarker;
        targetMarker.position = position;
        const auto it = std::lower_bound(orientedReadMarkers.begin(), orientedReadMarkers.end(), targetMarker);
        if((it != orientedReadMarkers.end()) and (it->position == position)) {
            ordinal = uint32_t(it - orientedReadMarkers.begin());
        }

        // Get the Kmer at this position, without relying on the ordinal.
        Kmer kmer;
        bool kmerIAvailable = false;
        try {
            kmer = reads.getKmer(k, orientedReadId, position);
            kmerIAvailable = true;
        } catch (const std::exception&) {
        }

        s << orientedReadId << " position " << position << ", ordinal ";
        if(ordinal == invalid<uint32_t>) {
            s << "(no ordinal)";
        } else {
            s << ordinal;
        }

        s << ", k-mer ";
        if(kmerIAvailable) {
            kmer.write(s, k);
        } else {
            s << "(not available)";
        }
        s << endl;

    }
}
