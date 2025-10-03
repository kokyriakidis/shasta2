#include "ExternalAnchors.hpp"
using namespace shasta;



ExternalAnchors::ExternalAnchors(const string& name)
{
    data.createNew(name, pageSize);
}



void ExternalAnchors::beginNewAnchor()
{
    data.appendVector();
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
