#pragma once

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "SegmentStepSupport.hpp"
#include "Tangle.hpp"

// Standard library.
#include "fstream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    class StrandSplitter;

    class AssemblyGraph;
}



// This takes as input a self-complementary tangle
// in an AssemblyGraph and attempts to split the strands.
class shasta2::StrandSplitter {
public:


    // The last two arguments are only used for debug output.
    StrandSplitter(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
        uint64_t tangleId,
        const string& debugOutputBaseName);

    ~StrandSplitter();

private:

    // EXPOSE WHEN CODE STABILIZES.
    const double maxCoverage = 16.;

    AssemblyGraph& assemblyGraph;
    uint64_t id(Segment) const;

    // If debug is set to true, debug output is written to html.
    bool debug = true;
    uint64_t tangleId;
    string debugOutputBaseName;
    ofstream html;
    void writeInitialDebugOutput();

    Tangle tangle;

    bool isEntrance(Segment) const;
    bool isExit(Segment) const;


    // The Tangle Segments.
    void gatherSegments();

    // All Tangle Segments (entrances, exits, and internal segments).
    // Sorted by id.
    vector<Segment> allTangleSegments;

    // Low coverage segments are the ones with coverage
    // no greater than maxCoverage. These are the ones
    // that are considered reliably single-copy and are used
    // for strand separation.
    // They don't necesserily include entrances and exits,
    // but they can.

    // Reverse complemented pairs of low coverage segments.
    // In each pair, the id of thefirst segment is less
    // than the id of the second segment.
    // Sorted by id of the first segment in the pair.
    vector< pair<Segment, Segment> > lowCoverageSegmentPairs;
    void writeLowCoverageSegmentPairs();

    // Low coverage segments, ordered by their appearance in segmentPairs.
    // This way, if a segment has index i in this vector, its
    // reverse complement has index i^1.
    vector<Segment> lowCoverageSegments;



    // Gather occurrences of reads in the first Segment of each pair.
    void findReadOccurrences();
    class ReadOccurrence {
    public:
        uint64_t segmentPairIndex = invalid<uint64_t>;
        Strand strand = invalid<Strand>;
        uint64_t frequency = 0;
        ReadOccurrence() {}
        ReadOccurrence(uint64_t segmentPairIndex, Strand strand) :
            segmentPairIndex(segmentPairIndex), strand(strand) {}
        bool operator==(const ReadOccurrence& that) const
        {
            return tie(segmentPairIndex, strand) == tie(that.segmentPairIndex, that.strand);
        }
        bool operator<(const ReadOccurrence& that) const
        {
            return tie(segmentPairIndex, strand) < tie(that.segmentPairIndex, that.strand);
        }
    };
    std::map<ReadId, vector<ReadOccurrence> > readOccurrenceMap;



    // For strand separation we construct an undirected graph
    // with a vertex for each of the segments in the lowCoverageSegments.
    // The vertex_descriptor of this graph is the index of the Segment
    // in the segments vector.
    // An edge segment0 to segment1 is created if there are reads that appear in the
    // same strand on segment0 and segment1.
    class StrandSeparationVertex {
    public:
        Segment segment;
        uint64_t component = invalid<uint64_t>;
        StrandSeparationVertex(Segment segment = assemblyGraphNullEdge) : segment(segment) {}
    };
    class StrandSeparationEdge {
    public:
        uint64_t frequency;
        bool isCrossStrandEdge = false;
        StrandSeparationEdge(uint64_t frequency) : frequency(frequency) {}
    };
    using StrandSeparationGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::vecS,
        boost::undirectedS,
        StrandSeparationVertex,
        StrandSeparationEdge>;
    class StrandSeparationGraph: public StrandSeparationGraphBaseClass {
    public:
        void addToEdge(
            uint64_t segmentIndex0,
            uint64_t segmentIndex1,
            uint64_t frequency);

        // Pairs of reverse complemented edges in the Graph,
        // sorted by decreasing frequency.
        class EdgePair {
        public:
            StrandSeparationGraph::edge_descriptor e;
            StrandSeparationGraph::edge_descriptor eRc;
            uint64_t frequency;
            bool operator<(const EdgePair& that) const
            {
                return frequency > that.frequency;
            }
        };
        vector<EdgePair> edgePairs;
        void findEdgePairs();
        void writeGraphviz(const string& fileName, const AssemblyGraph&) const;
    };
    StrandSeparationGraph strandSeparationGraph;
    void createStrandSeparationGraph();

    // Use the StrandSeparationGraph to separate strands.
    // If successful, this stores the strandSegments vectors.
    // strandSegments[0] are the segment in the first strand
    // which are used in the rest of the process.
    // They are stored sorted by id.
    bool separateStrands();
    array< vector<Segment>, 2> strandSegments;
    bool isStrand0Segment(Segment) const;
    bool isStrand1Segment(Segment) const;



    // In the ConnectionGraph, each vertex represents a Segment.
    // There is a vertex for each strand 0 Segment.
    // Edges correspond to connections already present in the
    // AssemblyGraph or additional connections that can be
    // made to split our Tangle.

    class ConnectionVertex {
    public:
        Segment segment;
        bool isEntrance;
        bool isExit;
        ConnectionVertex(
            Segment segment,
            bool isEntrance,
            bool isExit) :
            segment(segment),
            isEntrance(isEntrance),
            isExit(isExit)
        {}
    };

    class ConnectionEdge {
    public:

        // If this is true, a connection between these two segments
        // is already present in the AssemblyGraph.
        bool isDirectConnection;

        // The remaining fields are only filled in if directConnection is false;
        SegmentPairInformation segmentPairInformation;
        ConnectionEdge() : isDirectConnection(true) {}
        ConnectionEdge(const SegmentPairInformation& segmentPairInformation) :
            isDirectConnection(false),
            segmentPairInformation(segmentPairInformation) {}
    };

    using ConnectionGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        ConnectionVertex,
        ConnectionEdge>;
    class ConnectionGraph: public ConnectionGraphBaseClass {
    public:
        void addVertex(Segment, bool isEntrance, bool isExit);
        std::map<Segment, vertex_descriptor> vertexMap;
        void writeGraphviz(const string& fileName, const AssemblyGraph&) const;
    };
    ConnectionGraph connectionGraph;
    void createConnectionGraph();
};
