#pragma once

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "string.hpp"



// Class ExternaAnchors can be used to define Anchors to be used
// by shasta2 instead of the ones using the standard
// anchor generation process.

namespace shasta {
    class ExternalAnchors;
}



class shasta::ExternalAnchors {
public:
    // This will generate two files with the specified name
    // and extensions .toc and .data.
    ExternalAnchors(const string& name);

    // This is used to access an existing ExternalAnchors.
    class AccessExisting {};
    ExternalAnchors(const string& name, const AccessExisting&);

    // This is called to begin the definition of a new Anchor.
    void beginNewAnchor();

    // Add an  oriented read to the last anchor.
    // Read orientation is specified via Strand,
    // an integer which can be 0 (no reverse complementing)
    // or 1 (reverse complemented read).
    // The position is the position in base space of the
    // first Anchor Base in the oriented read.
    void addOrientedRead(ReadId, Strand, uint32_t position);

    class OrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t position;
        OrientedRead();
        OrientedRead(ReadId, Strand, uint32_t position);
    };

    const uint64_t pageSize = 4096;
    MemoryMapped::VectorOfVectors<OrientedRead, uint64_t> data;
};
