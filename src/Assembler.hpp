#pragma once

// Shasta.
#include "AssemblerOptions.hpp"
#include "HttpServer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"

namespace shasta {

    class Assembler;
    class AssemblyGraphPostprocessor;
    class Anchors;
    class AssemblerInfo;
    class Journeys;
    class KmerChecker;
    class KmersOptions;
    class LongBaseSequences;
    class Markers;
    class MarkerKmers;
    class Mode3Assembler;
    class Reads;


    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.


    extern template class MultithreadedObject<Assembler>;
}



// Class used to store various pieces of assembler information in shared memory.
class shasta::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;

    // The page size in use for this run.
    size_t largeDataPageSize;
};



class shasta::Assembler :
    public MultithreadedObject<Assembler>,
    public MappedMemoryOwner,
    public HttpServer {
public:


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor for a new assembler also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Construct a new Assembler.
    Assembler(
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize);

    // Construct an Assembler from binary data. This accesses the AssemblerInfo and the Reads.
    Assembler(const string& largeDataFileNamePrefix);



    // Various pieces of assembler information stored in shared memory.
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;

    // This runs the entire assembly, under the following assumptions:
    // - The current directory is the run directory.
    // - The Data directory has already been created and set up, if necessary.
    // - The input file names are either absolute,
    //   or relative to the run directory, which is the current directory.
    void assemble(
        const AssemblerOptions& assemblerOptions,
        vector<string> inputFileNames);


    // Reads.
    shared_ptr<Reads> readsPointer;
    const Reads& reads() const {
        SHASTA_ASSERT(readsPointer);
        return *readsPointer;
    }
    void computeReadIdsSortedByName();
    void addReads(
        const vector<string>& fileNames,
        uint64_t minReadLength,
        size_t threadCount);
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        size_t threadCount);
    void histogramReadLength(const string& fileName);



    // The KmerChecker is used to find out if a given KmerId is a marker.
    shared_ptr<KmerChecker> kmerChecker;
    public:
    void createKmerChecker(
        uint64_t k,
        double markerDensity,
        uint64_t threadCount);
    void accessKmerChecker();



    // The markers on all oriented reads.
    shared_ptr<Markers> markersPointer;
    const Markers& markers() const
    {
        SHASTA_ASSERT(markersPointer);
        return *markersPointer;
    }
    void checkMarkersAreOpen() const;
    void createMarkers(size_t threadCount);
    void accessMarkers();



    // The MarkerKmers keep track of the locations in the oriented reads
    // where each marker k-mer appears.
    shared_ptr<MarkerKmers> markerKmers;
    void createMarkerKmers(uint64_t threadCount);
    void accessMarkerKmers();



    // Anchors.
    shared_ptr<Anchors> anchorsPointer;
    const Anchors& anchors() const
    {
        SHASTA_ASSERT(anchorsPointer);
        return *anchorsPointer;
    }
    void createAnchors(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t maxHomopolymerLength,
        uint64_t threadCount);
    void accessAnchors(bool writeAccess);



    // Journeys.
    shared_ptr<Journeys> journeysPointer;
    const Journeys& journeys() const
    {
        SHASTA_ASSERT(journeysPointer);
        return *journeysPointer;
    }
    void createJourneys(uint64_t threadCount);
    void accessJourneys();


    // AssemblyGraph.
    void createAssemblyGraph(
        const AssemblerOptions&,
        uint64_t threadCount);



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
    void accessAllSoft();


    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(const vector<string>&, ostream&);
    void exploreSegments(const vector<string>&, ostream&);
    void exploreSegment(const vector<string>&, ostream&);
    void exploreVertexTangle(const vector<string>&, ostream&);
    void exploreEdgeTangle(const vector<string>&, ostream&);
    void exploreTangleMatrix(const vector<string>&, ostream&);
    void exploreLocalAssembly(const vector<string>&, ostream&);
    void exploreLocalAssembly1(const vector<string>&, ostream&);

    // Get the AssemblyGraph for a given assembly stage.
    AssemblyGraphPostprocessor& getAssemblyGraph(
        const string& assemblyStage,
        const AssemblerOptions&);
    std::map<string, shared_ptr<AssemblyGraphPostprocessor> > assemblyGraphTable;

};

