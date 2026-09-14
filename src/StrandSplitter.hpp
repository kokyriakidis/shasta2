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
    // with a vertex for each of the segments in the segmentsVector.
    // The vertex_descriptor of this graph is the index of the Segment
    // in the segments vector.
    // An edge s0-s1 is created if there are reads that appear in the
    // same strand on s0 and s1.
    class Vertex {
    public:
        Segment segment;
        uint64_t component = invalid<uint64_t>;
        Vertex(Segment segment = assemblyGraphNullEdge) : segment(segment) {}
    };
    class Edge {
    public:
        uint64_t frequency;
        bool isCrossStrandEdge = false;
        Edge(uint64_t frequency) : frequency(frequency) {}
    };
    using GraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::vecS,
        boost::undirectedS,
        Vertex,
        Edge>;
    class Graph: public GraphBaseClass {
    public:
        void addToEdge(
            uint64_t segmentIndex0,
            uint64_t segmentIndex1,
            uint64_t frequency);

        // Pairs of reverse complemented edges in the Graph,
        // sorted by decreasing frequency.
        class EdgePair {
        public:
            Graph::edge_descriptor e;
            Graph::edge_descriptor eRc;
            uint64_t frequency;
            bool operator<(const EdgePair& that) const
            {
                return frequency > that.frequency;
            }
        };
        vector<EdgePair> edgePairs;
        void findEdgePairs();
    };
    Graph graph;
    void createGraph();

    // Use the Graph to separate strands.
    // if successful, this stores the segments attribute to each thread.
    // They are stored sorted by id.
    bool separateStrands();
    array< vector<Segment>, 2> strandSegments;
    bool isStrand0Segment(Segment) const;
    bool isStrand1Segment(Segment) const;

    // A forward hanging segment is a strand 0 segment or an entrance that is not
    // immediately followed by at least another strand0 segment
    // or an exit.
    // A backward orphan segment is a strand 0 segment or an exit that is not
    // immediately preceded by at least another strand0 segment
    // or an entrance.
    array< vector<Segment>, 2> hangingSegments;  // 0 = forward, 1 = backward.
    void findHangingSegments();

    // Candidate connections between strand 0 segments are found using
    // forward BFS from the forward hanging segments
    // and backward BFS from the backward hanging segments.
    // The BFSs are not allowed to use strand 1 segments that are not
    // entrances or exits,
    // and stop when a strand0 segment or an entrance or an exit is found.
    // In addition, there are direct connections, which are
    // connections segment0->segment1 where the target vertex
    // of segment0 is the same as the source vertex of segment1.
    class CandidateConnection : public pair<Segment, Segment> {
    public:
        bool isDirectConnection;

        // The remaining fields are only filled in if directConnection is false;
        SegmentPairInformation segmentPairInformation;
        bool canConnect = false;
        bool canConnectDeep = false;
        CandidateConnection(Segment segment0, Segment segment1, bool isDirectConnection) :
            pair<Segment, Segment>(segment0, segment1), isDirectConnection(isDirectConnection) {}
    };
    vector<CandidateConnection> candidateConnections;
    void findCandidateConnections();
};
