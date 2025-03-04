#ifndef SHASTA_READS
#define SHASTA_READS

// Shasta
#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "ReadFlags.hpp"
#include "shastaTypes.hpp"
#include "SHASTA_ASSERT.hpp"
#include "span.hpp"

namespace shasta {
    class Reads;
    class OrientedReadId;
}



// The reads used for this assembly, iIndexed by ReadId.

class shasta::Reads {
public:
  
    // Default Constructor
    Reads(): totalBaseCount(0), n50(0) {};

    void createNew(
        const string& readsDataName,
        const string& readNamesDataName,
        const string& readMetaDataDataName,
        const string& readFlagsDataName,
        const string& readIdsSortedByNameDataName,
        uint64_t largeDataPageSize
    );

    void access(
        const string& readsDataName,
        const string& readNamesDataName,
        const string& readMetaDataDataName,
        const string& readFlagsDataName,
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

    inline span<const char> getReadMetaData(ReadId readId) const {
        return readMetaData[readId];
    }

    inline const ReadFlags& getFlags(ReadId readId) const {
        return readFlags[readId];
    }
    inline const MemoryMapped::Vector<ReadFlags>& getFlags() const {
        return readFlags;
    }

    
    Base getOrientedReadBase(
        OrientedReadId orientedReadId,
        uint32_t position
    ) const;

    // Return a vector containing the raw sequence of an oriented read.
    vector<Base> getOrientedReadRawSequence(OrientedReadId) const;

    // Return the length of the raw sequence of a read.
    // If using the run-length representation of reads, this counts each
    // base a number of times equal to its repeat count.
    uint64_t getReadRawSequenceLength(ReadId) const;

    // Get a vector of the raw read positions
    // corresponding to each position in the run-length
    // representation of an oriented read.
    vector<uint32_t> getRawPositions(OrientedReadId) const;

    // Return a meta data field for a read, or an empty string
    // if that field is missing. This treats the meta data
    // as a space separated sequence of Key=Value,
    // without embedded spaces in each Key=Value pair.
    span<const char> getMetaData(ReadId, const string& key) const;


    // Setters for readFlags.
    inline void setPalindromicFlag(ReadId readId, bool value) {
        readFlags[readId].isPalindromic = value;
    }
    inline void setChimericFlag(ReadId readId, bool value) {
        readFlags[readId].isChimeric = value;
    }
    inline void setStrandFlag(ReadId readId, bool value) {
        readFlags[readId].strand = value;
    }

    // Function to write one or all reads in Fasta format.
    void writeReads(const string& fileName);
    void writeRead(ReadId, const string& fileName);
    void writeOrientedRead(ReadId, Strand, const string& fileName);

    void writeRead(ReadId, ostream&);
    void writeOrientedRead(OrientedReadId, ostream&);
    void writeOrientedRead(OrientedReadId, const string& fileName);


    // Assertions for data integrity.
    inline void checkSanity() const {
        SHASTA_ASSERT(readNames.size() == reads.size());
        SHASTA_ASSERT(readMetaData.size() == reads.size());
    }

    inline void checkReadsAreOpen() const {
        SHASTA_ASSERT(reads.isOpen());
    }

    inline void checkReadNamesAreOpen() const {
        SHASTA_ASSERT(readNames.isOpen());
    }

    inline void checkReadMetaDataAreOpen() const {
        SHASTA_ASSERT(readMetaData.isOpen());
    }

    inline void checkReadFlagsAreOpen() const {
        SHASTA_ASSERT(readFlags.isOpen);
    }

    inline void checkReadFlagsAreOpenForWriting() const {
        SHASTA_ASSERT(readFlags.isOpenWithWriteAccess);
    }

    inline void checkReadId(ReadId readId) const {
        if (readId >= reads.size()) {
            throw runtime_error("Read id " + to_string(readId) +
                " is not valid. Must be between 0 and " + to_string(reads.size()) +
                " inclusive.");
        }
    }

    inline void assertReadsAndFlagsOfSameSize() const {
        SHASTA_ASSERT(reads.size() == readFlags.size());
    }

    void computeReadLengthHistogram();
    void writeReadLengthHistogram(const string& fileName);
    
    inline uint64_t getTotalBaseCount() const {
        return totalBaseCount;
    }
    inline uint64_t getN50() const {
        return n50;
    }

    inline const vector<uint64_t>& getReadLengthHistogram() const {
        return histogram;
    }

    void rename();

    void copyDataForReadsLongerThan(
        const Reads& rhs,
        uint64_t newMinReadLength,
        uint64_t& discardedShortReadCount,
        uint64_t& discardedShortReadBases
    );

    // Find duplicate reads, as determined by name (not sequence).
    // This also sets the isDuplicate and discardDueToDuplicates read flags
    // and summarizes what it found Duplicates.csv.
    void findDuplicates(const string& handleDuplicates);

    void remove();

    uint64_t representation; // 0 = raw sequence, 1 = RLE sequence
private:
    LongBaseSequences reads;

    // The names of the reads from the input fasta or fastq files.
    // Indexed by ReadId.
    // Note that we don't enforce uniqueness of read names.
    // We don't use read names to identify reads.
    // These names are only used as an aid in tracing each read
    // back to its origin.
    MemoryMapped::VectorOfVectors<char, uint64_t> readNames;

    // Read meta data. This is the information following the read name
    // in the header line for fasta and fastq files.
    // Indexed by ReadId.
    MemoryMapped::VectorOfVectors<char, uint64_t> readMetaData;

    MemoryMapped::Vector<ReadFlags> readFlags;



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

#endif
