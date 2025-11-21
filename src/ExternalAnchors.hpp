#pragma once

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "string.hpp"



// Class ExternaAnchors can be used to define Anchors to be used
// by shasta2 instead of the ones using the standard
// anchor generation process.

namespace shasta2 {
    class ExternalAnchors;
    class Markers;
    class Reads;
}



class shasta2::ExternalAnchors {
public:
    // This will generate four files:
    // name.toc name.data, name-Names.toc, name-names.data.
    ExternalAnchors(const string& name);

    // This is called to begin the definition of a new Anchor.
    void beginNewAnchor(const string& anchorName);

    // Add an  oriented read to the last anchor.
    // Read orientation is specified via Strand,
    // an integer which can be 0 (no reverse complementing)
    // or 1 (reverse complemented read).
    // The position is the position in base space of the
    // first Anchor Base in the oriented read.
    void addOrientedRead(ReadId, Strand, uint32_t position);

    // THE PYTHON API ENDS HERE.

    // This is used to access an existing ExternalAnchors.
    class AccessExisting {};
    ExternalAnchors(const string& name, const AccessExisting&);


    class OrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t position;
        OrientedRead();
        OrientedRead(ReadId, Strand, uint32_t position);
    };

    const uint64_t pageSize = 4096;

    // The definitions of the external anchors.
    MemoryMapped::VectorOfVectors<OrientedRead, uint64_t> data;

    // External anchor names, only used for diagnostics.
    MemoryMapped::VectorOfVectors<char, uint64_t> names;

    // Write information about the i-th external anchor.
    void write(ostream&, uint64_t i, uint64_t k, const Reads&, const Markers&) const;
};
