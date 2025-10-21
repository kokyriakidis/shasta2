#include "ExternalAnchors.hpp"
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

