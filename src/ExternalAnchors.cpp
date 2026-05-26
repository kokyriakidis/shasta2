#include "ExternalAnchors.hpp"
#include "invalid.hpp"
#include "Reads.hpp"
using namespace shasta2;



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
    const Reads& reads) const
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
        const Kmer kmer = reads.getKmer(k, orientedReadId, position);

        s << orientedReadId << " position " << position << ", k-mer ";
        kmer.write(s, k);
        s << endl;
    }
}
