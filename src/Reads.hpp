#pragma once

// Shasta
#include "Base.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "shastaTypes.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "span.hpp"

namespace shasta2 {
    class Reads;
    class OrientedReadId;
}



// The reads used for this assembly, indexed by ReadId.

class shasta2::Reads {
public:
  
    // Default Constructor
    Reads(): n50(0) {};

    void createNew(
        const string& readsDataName,
        const string& readNamesDataName,
        const string& readIdsSortedByNameDataName,
        uint64_t largeDataPageSize
    );

    void access(
        const string& readsDataName,
        const string& readNamesDataName,
        const string& readIdsSortedByNameDataName
    );

    inline ReadId readCount() const {
        return ReadId(reads.size());
    }

    inline LongBaseSequenceView getRead(ReadId readId) const {
        return reads[readId];
    }

    inline span<const char> getReadName(ReadId readId) const {
        return readNames[readId];
    }

    // Get a ReadId given a read name.
    // This uses a binary search in readIdsSortedByName.
    ReadId getReadId(const string& readName) const;
    ReadId getReadId(const span<const char>& readName) const;

    Base getOrientedReadBase(
        OrientedReadId orientedReadId,
        uint32_t position
    ) const;

    // Return a vector containing the sequence of an oriented read.
    vector<Base> getOrientedReadSequence(OrientedReadId) const;

    // Return the length of the sequence of a read.
    uint64_t getReadSequenceLength(ReadId) const;

    // Return a meta data field for a read, or an empty string
    // if that field is missing. This treats the meta data
    // as a space separated sequence of Key=Value,
    // without embedded spaces in each Key=Value pair.
    span<const char> getMetaData(ReadId, const string& key) const;

    // Assertions for data integrity.
    inline void checkSanity() const {
        SHASTA2_ASSERT(readNames.size() == reads.size());
    }
    inline void checkReadsAreOpen() const {
        SHASTA2_ASSERT(reads.isOpen());
    }
    inline void checkReadNamesAreOpen() const {
        SHASTA2_ASSERT(readNames.isOpen());
    }

    inline void checkReadId(ReadId readId) const {
        if (readId >= reads.size()) {
            throw runtime_error("Read id " + to_string(readId) +
                " is not valid. Must be between 0 and " + to_string(reads.size()) +
                " inclusive.");
        }
    }

    void computeReadLengthHistogram();
    void writeReadLengthHistogram(const string& fileName);
    
    inline uint64_t getTotalBaseCount() const {
        return reads.totalBaseCount();
    }
    inline uint64_t getN50() const {
        return n50;
    }

    inline const vector<uint64_t>& getReadLengthHistogram() const {
        return histogram;
    }

    void remove();

    Kmer getKmer(uint64_t k, OrientedReadId, uint32_t position) const;

private:
    LongBaseSequences reads;

    // The names of the reads from the input fasta or fastq files.
    // Indexed by ReadId.
    // Note that we don't enforce uniqueness of read names.
    // We don't use read names to identify reads.
    // These names are only used as an aid in tracing each read
    // back to its origin.
    MemoryMapped::VectorOfVectors<char, uint64_t> readNames;

    // The read ids, sorted by name.
    // This is used to find the read id corresponding to a name.
    MemoryMapped::Vector<ReadId> readIdsSortedByName;
public:
    void computeReadIdsSortedByName();
private:

    // Class used to sort ReadId's by name, and also to look up the
    // ReadId corresponding to a given name.
    class OrderReadsByName {
    public:
        OrderReadsByName(const MemoryMapped::VectorOfVectors<char, uint64_t>& readNames) :
            readNames(readNames) {}

        // This one is used by std::sort.
        bool operator()(const ReadId& readId0, const ReadId& readId1) const {
            const auto name0 = readNames[readId0];
            const auto name1 = readNames[readId1];
            return std::lexicographical_compare(name0.begin(), name0.end(), name1.begin(), name1.end());
        }

        // This one is used by std::lower_bound.
        bool operator()(const ReadId& readId0, const span<const char>& name1) const {
            const auto name0 = readNames[readId0];
            return std::lexicographical_compare(name0.begin(), name0.end(), name1.begin(), name1.end());
        }
    private:
        const MemoryMapped::VectorOfVectors<char, uint64_t>& readNames;
    };

    
    // Read statistics.
    vector<uint64_t> histogram;
    vector< pair<uint64_t, uint64_t> > binnedHistogram;
    uint64_t totalBaseCount;
    uint64_t n50;    

    friend class ReadLoader;
};

