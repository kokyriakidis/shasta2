#pragma once

// Shasta.
#include "Anchor.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include <cstdint.hpp>
#include <span.hpp>

namespace shasta {
    class Journeys;
    using Journey = span<const AnchorId>;
}



// The journey of each oriented read is the sequence of AnchorIds
// encountered by the oriented read.
class shasta::Journeys :
    public MultithreadedObject<Journeys>,
    public MappedMemoryOwner {
public:

    // Initial creation.
    // This sets the positionInJourney for every AnchorMarkerInfo
    // stored in the Anchors, and for this reason the Anchors
    // are not passed in as const.
    Journeys(
        uint64_t orientedReadCount,
        shared_ptr<Anchors>,
        uint64_t threadCount,
        const MappedMemoryOwner&);

    // Access from binary data.
    Journeys(const MappedMemoryOwner&);

    // Return the Journey for an oriented read.
    Journey operator[](OrientedReadId orientedReadId) const
    {
        return journeys[orientedReadId.getValue()];
    }

    uint64_t size() const
    {
        return journeys.size();
    }

    void writeAnchorGapsByRead(
        const Reads&,
        const Markers&,
        const Anchors&
        ) const;

    bool isOpen() const
    {
        return journeys.isOpen();
    }

private:

    MemoryMapped::VectorOfVectors<AnchorId, uint64_t> journeys;

    // This is only stored during initial creation.
    shared_ptr<Anchors> anchorsPointer;

    void threadFunction1(uint64_t threadId);
    void threadFunction2(uint64_t threadId);
    void threadFunction12(uint64_t pass);
    void threadFunction3(uint64_t threadId);
    void threadFunction4(uint64_t threadId);

    // Temporary storage of journeys with ordinals.
    MemoryMapped::VectorOfVectors<pair<uint64_t, uint32_t>, uint64_t> journeysWithOrdinals;

};
