#ifndef SHASTA_ASSEMBLER_HPP
#define SHASTA_ASSEMBLER_HPP

// Shasta.
#include "AssemblerOptions.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "Marker.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"

namespace shasta {

    class Assembler;
    class AssemblerInfo;
    class AssemblerOptions;
    class AssembledSegment;
    class Histogram2;
    class KmerChecker;
    class KmersOptions;
    class LongBaseSequences;
    class MarkerKmers;
    class Mode2AssemblyOptions;
    class Mode3AssemblyOptions;
    class Mode3Assembler;
    class Reads;

    namespace mode3 {
        class Anchors;
    }

    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.


    extern template class MultithreadedObject<Assembler>;
}


class DisjointSets;





// Class used to store various pieces of assembler information in shared memory.
class shasta::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;

    // The method used to generate kmers (--Kmers.generationMethod).
    uint64_t kmerGenerationMethod;

    // The page size in use for this run.
    size_t largeDataPageSize;

    // Assembly mode (0=haploid, 2=phased).
    uint64_t assemblyMode;



    // Statistics on the number of reads discarded on input.
    // These are incremented during each call to addReadsFromFasta.

    // The number of reads and raw bases discarded because the read
    // contained invalid bases.
    uint64_t discardedInvalidBaseReadCount = 0;
    uint64_t discardedInvalidBaseBaseCount = 0; // Only counts the valid bases in those reads.

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount = 0;
    uint64_t discardedShortReadBaseCount = 0;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount = 0;
    uint64_t discardedBadRepeatCountBaseCount = 0;



    // Statistics for the reads kept in the assembly
    // and not discarded on input.
    // These are computed and stored by histogramReadLength.
    size_t readCount = 0;
    size_t baseCount = 0;
    size_t readN50 = 0;
    uint64_t minReadLength = 0;

    // Other read statistics.
    size_t palindromicReadCount = 0;
    size_t chimericReadCount = 0;
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;

};



class shasta::Assembler :
    public MultithreadedObject<Assembler>,
    public MappedMemoryOwner,
    public HttpServer {
public:


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    1GB huge pages are also supported.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Constructor to be called one to create a new run.
    Assembler(
        const string& largeDataFileNamePrefix,
        bool createNew,
        size_t largeDataPageSize);

    // Add reads.
    // The reads in the specified file are added to those already previously present.
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        bool noCache,
        size_t threadCount);

    // Create a histogram of read lengths.
    void histogramReadLength(const string& fileName);


    // Functions related to markers.
    // See the beginning of Marker.hpp for more information.
    void findMarkers(size_t threadCount);
    void accessMarkers();
    void writeMarkers(ReadId, Strand, const string& fileName);







    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal of the marker in the oriented read.
    MarkerGraphVertexId getGlobalMarkerGraphVertex(
        ReadId,
        Strand,
        uint32_t ordinal) const;





    // Reads.
    shared_ptr<Reads> reads;
public:
    const Reads& getReads() const {
        SHASTA_ASSERT(reads);
        return *reads;
    }

    uint64_t adjustCoverageAndGetNewMinReadLength(uint64_t desiredCoverage);


    void computeReadIdsSortedByName();

    // Find duplicate reads, as determined by name (not sequence).
    // This also sets the isDuplicate and discardDueToDuplicates read flags
    // and summarizes what it found Duplicates.csv.
    void findDuplicateReads(const string& handleDuplicates);


private:



    // Various pieces of assembler information stored in shared memory.
    // See class AssemblerInfo for more information.
public:
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;


    // The KmerChecker can find out if a given KmerId is a marker.
    shared_ptr<KmerChecker> kmerChecker;
    public:
    void createKmerChecker(
        const KmersOptions& kmersOptions,
        uint64_t threadCount);
    void accessKmerChecker();

    // This one should eventually go away, but there are several scripts
    // that depend on it.
    void accessKmers()
    {
        accessKmerChecker();
    }


private:


    // Hash a KmerId in such a way that it has the same hash as its reverse
    // complement. This is used by alignment method 3 to downsample markers.
    uint32_t hashKmerId(KmerId) const;

    // The markers on all oriented reads. Indexed by OrientedReadId::getValue().
public:
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> markers;
private:
    void checkMarkersAreOpen() const;

    // Get markers sorted by KmerId for a given OrientedReadId.
    void getMarkersSortedByKmerId(
        OrientedReadId,
        vector<MarkerWithOrdinal>&) const;

    // Given a marker by its OrientedReadId and ordinal,
    // return the corresponding global marker id.
public:
    MarkerId getMarkerId(OrientedReadId, uint32_t ordinal) const;
private:
    MarkerId getReverseComplementMarkerId(OrientedReadId, uint32_t ordinal) const;
    MarkerId getMarkerId(const MarkerDescriptor& m) const
    {
        return getMarkerId(m.first, m.second);
    }
    MarkerId getReverseComplementMarkerId(const MarkerDescriptor& m) const
    {
        return getReverseComplementMarkerId(m.first, m.second);
    }

    // Inverse of the above: given a global marker id,
    // return its OrientedReadId and ordinal.
    // This requires a binary search in the markers toc.
    // This could be avoided, at the cost of storing
    // an additional 4 bytes per marker.
public:
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;


    // KmerIds for all markers. Indexed by OrientedReadId::getValue().
    // Only stored during alignment computation, and then freed.
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> markerKmerIds;
    void computeMarkerKmerIds(uint64_t threadCount);
    void cleanupMarkerKmerIds();
private:
    void computeMarkerKmerIdsThreadFunction(size_t threadId);


    // Pairs (KmerId, ordinal), sorted by KmerId, for each oriented read.
    // Indexed by orientedReadId.getValue().
    // Used by alignment method 4.
public:
    MemoryMapped::VectorOfVectors< pair<KmerId, uint32_t>, uint64_t> sortedMarkers;
    void computeSortedMarkers(uint64_t threadCount);
    bool accessSortedMarkers();
private:
    void computeSortedMarkersThreadFunction(size_t threadId);
    // void computeSortedMarkersThreadFunction1(size_t threadId);
    // void computeSortedMarkersThreadFunction2(size_t threadId);



    // Low frequency markers for each oriented read.
    // This stores, for each oriented read, the ordinals corresponding
    // to marker with low frequency (up to maxMarkerFrequency), sorted by KmerId.
    // Used by alignment method 5. It is only stored durign alignment
    // computation.
public:
    MemoryMapped::VectorOfVectors<uint32_t, uint64_t> lowFrequencyMarkers;
    void computeLowFrequencyMarkers(uint64_t maxMarkerFrequency, uint64_t threadCount);
    void computeLowFrequencyMarkers(
        const span<const KmerId>&,  // The marker k-mers for the oriented reads (sorted by ordinal)
        uint64_t maxMarkerFrequency,
        vector<uint32_t>&);         // The ordinals of the low frequency markers, sorted by KmerId.
private:
    void computeLowFrequencyMarkersThreadFunctionPass1(uint64_t threadId);
    void computeLowFrequencyMarkersThreadFunctionPass2(uint64_t threadId);
    void computeLowFrequencyMarkersThreadFunctionPass12(uint64_t pass);
    class ComputeLowFrequencyMarkersData {
    public:
        uint64_t maxMarkerFrequency;
    };
    ComputeLowFrequencyMarkersData computeLowFrequencyMarkersData;



    // Low level functions to get marker Kmers/KmerIds of an oriented read.
    // They are obtained from the reads and not from CompressedMarker::kmerId,
    // which will soon go away.

    // Get the marker Kmer for an oriented read and ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand0(ReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand1(ReadId, uint32_t ordinal) const;

    // Get the marker KmerId for an oriented read and ordinal.
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint32_t ordinal) const;

    // Get all marker Kmers for an oriented read.
    void getOrientedReadMarkerKmers(OrientedReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand0(ReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand1(ReadId, const span<Kmer>&) const;

    // Get all marker KmerIds for an oriented read.
    void getOrientedReadMarkerKmerIds(OrientedReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand0(ReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand1(ReadId, const span<KmerId>&) const;

    // Get all MarkerWithOrdinals for an oriented read (includes position, KmerId, and ordinal).
    void getOrientedReadMarkers(OrientedReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand0(ReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand1(ReadId, const span<MarkerWithOrdinal>&) const;

    // Get all marker Kmers/KmerIds for a read in both orientations.
    void getReadMarkerKmers(
        ReadId,
        const span<Kmer>& Kmers0,
        const span<Kmer>& Kmers1) const;
    void getReadMarkerKmerIds(
        ReadId,
        const span<KmerId>& kmerIds0,
        const span<KmerId>& kmerIds1) const;

    // Get the Kmer/KmerId for an oriented read at a given marker ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint64_t ordinal) const;
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint64_t ordinal) const;

    // Given a MarkerId, compute the MarkerId of the
    // reverse complemented marker.
    MarkerId findReverseComplement(MarkerId) const;


    // The MarkerKmers keep track of the locations in the oriented reads
    // where each marker k-mer appears. It is only used for alignment-free assembly
    // (--Assembly.mode 3 --Assembly.mode3.anchorCreationMethod FromMarkerKmers)
    // but can also be created using CreateMarkerKmers.py.
    shared_ptr<MarkerKmers> markerKmers;
public:
    void createMarkerKmers(uint64_t threadCount);
    void accessMarkerKmers();
private:

    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&) const;
    void writeHtmlEnd(ostream&) const;
    static void writeStyle(ostream& html);


    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;

    static void writePngToHtml(
        ostream& html,
        const string& pngFileName,
        const string useMap = ""
        );
    static void writeGnuPlotPngToHtml(
        ostream& html,
        int width,
        int height,
        const string& gnuplotCommands);

    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreReadRaw(const vector<string>&, ostream&);
    void exploreLookupRead(const vector<string>&, ostream&);
    void exploreReadSequence(const vector<string>&, ostream&);
    void exploreReadMarkers(const vector<string>&, ostream&);
    void exploreMarkerKmers(const vector<string>&, ostream&);
    static void addScaleSvgButtons(ostream&, uint64_t sizePixels);

public:
    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;

        const AssemblerOptions* assemblerOptions = 0;
    };
    HttpServerData httpServerData;




    void writeMakeAllTablesCopyable(ostream&) const;


    // Access all available assembly data, without thorwing an exception
    // on failures.
public:
    void accessAllSoft();


    // Mode 3 assembly.
    shared_ptr<Mode3Assembler> mode3Assembler;
    void accessMode3Assembler();

    // Top level function for Mode 3 assembly.
    void mode3Assembly(
        uint64_t threadCount,
        shared_ptr<mode3::Anchors>,
        const Mode3AssemblyOptions&,
        bool debug
    );
    // Same, but use existing Anchors. Python callable.
    void mode3Reassembly(
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool debug
    );

    // Alignment-free version of mode 3 assembly.
    void alignmentFreeAssembly(
        const Mode3AssemblyOptions&,
        uint64_t threadCount);

    // Http server functions related to Mode 3 assembly.
    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAssembly(const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(const vector<string>&, ostream&);


public:
    void test();
};

#endif
