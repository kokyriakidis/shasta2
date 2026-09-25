#pragma once

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "graphvizToHtml.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include "string.hpp"
#include "tuple.hpp"
#include "vector.hpp"



namespace shasta2 {
    namespace StrandSeparation {
        class StrandContact;
        class SegmentPair;
        class ReadOccurrence;

        class BipartiteGraphVertex;
        class BipartiteGraphEdge;
        using BipartiteGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            BipartiteGraphVertex,
            BipartiteGraphEdge>;
        class BipartiteGraph;
    }

    class AssemblyGraph;
}



// A class to store a pair of Segments that are
// the reverse complement of each pair.
// They are stored with the lower id Segment first.
class shasta2::StrandSeparation::SegmentPair {
public:
    Segment segment0;
    Segment segment1;
};



// Class used to count occurrences of reads in Segments.
class shasta2::StrandSeparation::ReadOccurrence {
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



// A vertex of the BipartiteGraph can represent an OrientedReadId
// or a Segment.
class shasta2::StrandSeparation::BipartiteGraphVertex {
public:
    bool isSegment;

    OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);

    // Index of the segment in the StrandContact::lowCoverageSegments vector.
    uint64_t segmentIndex = invalid<uint64_t>;

    BipartiteGraphVertex(uint64_t segmentIndex) :
        isSegment(true),
        segmentIndex(segmentIndex)
    {}

    BipartiteGraphVertex(OrientedReadId orientedReadId) :
        isSegment(false),
        orientedReadId(orientedReadId) {}
};



class shasta2::StrandSeparation::BipartiteGraphEdge {
public:

    // The number of times (steps) this OrientedReadId appears
    // in this Segment.
    uint64_t frequency;
};



class shasta2::StrandSeparation::BipartiteGraph : public BipartiteGraphBaseClass {
public:
    // Map a low coverage segment index to a vertex_descriptor.
    vector<vertex_descriptor> segmentIndexToVertexMap;

    // Map an OrientedReadId to a vertex_descriptor.
    std::map<OrientedReadId, vertex_descriptor> orientedReadIdToVertexMap;

    // Map vertices to indexes.
    std::map<BipartiteGraph::vertex_descriptor, uint64_t> vertexIndexMap;

    // Map indexes to vertices.
    vector<BipartiteGraph::vertex_descriptor> vertexTable;

    // Get the reverse complement of a vertex.
    vertex_descriptor reverseComplement(vertex_descriptor) const;

    // Graphviz output.
    void writeGraphviz(
        const string& fileName,
        const vector<Segment>& lowCoverageSegments,
        const vector< vector<uint64_t> >& componentsIndexes,
        const AssemblyGraph&) const;

    // The pairs of reverse complemented edges.
    // Sorted by decreasing frequency.
    class EdgePair {
    public:
        edge_descriptor e;
        edge_descriptor eRc;
        uint64_t frequency;
        bool operator<(const EdgePair& that) const
        {
            return frequency > that.frequency;
        }
    };
    vector<EdgePair> edgePairs;


    // One attempt at strand separation.
    // This add edge pairs in the given order to a disjoint sets
    // data structure. However an edge pair is discarded if it
    // would cause contacts between strands.
    // This leaves the BipartiteGraph unchanged and returns
    // the connected components computed in this way,
    // sorted by decreasing size. Reverse complemented
    // connected components are guaranteed to be consecutive.
    // Vertex indexes in each component are sorted.
    // They can be converted to vertex_descriptors
    // via the vertexTable.
    void strandSeparation(
        const vector<uint64_t>& edgePairsIndexes,
        vector< vector<uint64_t> >& componentsIndexes) const;

    // High level function for strand separation.
    void strandSeparation(vector< vector<uint64_t> >& componentsIndexes);

};



class shasta2::StrandSeparation::StrandContact {
public:

    StrandContact(
        AssemblyGraph&,
        const vector<Segment>& allSegmentsById, // All Segments,sorted by id.
        const string& debugOutputBaseName,  // Only used for debug output.
        uint64_t strandContactId            // Only used for debug output.
        );

private:

    // EXPOSE WHEN CODE STABILIZES.
    const double coverageThreshold = 16.;
    const double strandFractionThreshold = 0.8;

    // Constructor arguments.
    AssemblyGraph& assemblyGraph;
    vector<Segment> allSegmentsById;
    const string& debugOutputBaseName;
    uint64_t strandContactId;

    ofstream html;

    uint64_t id(Segment) const;



    // SegmentPairs are stored sorted by id
    // of the first Segment in the SegmentPair.
    // The lowCoverageSegmentPairs are the ones with
    // coverage up to coverageThreshold. They are
    // the ones used for strand separation.
    vector<SegmentPair> allSegmentPairs;
    vector<SegmentPair> lowCoverageSegmentPairs;
    vector<SegmentPair> highCoverageSegmentPairs;

    // All the Segments, stored in the same order as they
    // appear in the lowCoverageSegmentPairs. This means that pairs
    // of reverse complemented Segments have consecutive indexes
    // in this vector.
    // This means that they are not sorted by id.
    vector<Segment> allSegments;

    // The low coverage Segments, stored in the same order as they
    // appear in the lowCoverageSegmentPairs. This means that pairs
    // of reverse complemented Segments have consecutive indexes
    // in this vector.
    // This means that they are not sorted by id.
    vector<Segment> lowCoverageSegments;

    void gatherSegments();
    void writeAllSegmentsById();
    void writeSegments();



    // Count occurrences of reads in the first Segment of each low coverage Segment pair.
    std::map<ReadId, vector<ReadOccurrence> > readOccurrenceMap;
    void countReadOccurrences();
    void writeReadOccurrences();


    // The BipartiteGraph has a vertex for each low coverage Segment
    // and a vertex for each OrientedReadId that appears
    // in the low coverage Segments.
    // An edge (segment--orientedReadId) contains the number of times (steps)
    // that orientedReadId occures in that segment.
    BipartiteGraph bipartiteGraph;
    void createBipartiteGraph();
    void writeBipartiteGraph(const vector< vector<uint64_t> >& componentsIndexes);
    void writeComponents(const vector< vector<uint64_t> >& componentsIndexes);

    // The low coverage Segments and OrientedReadIds of the first two components,
    // after strand separation.
    // These are assumed to define strand0 and strand1.
    // The two strandSegments vectors are sorted by id.
    array<vector<Segment>, 2> strandSegments;
    array<vector<OrientedReadId>, 2> strandOrientedReadIds;
    void gatherStrands(const vector< vector<uint64_t> >& componentsIndexes);



    // Classify segments.
    // There is an entry for each Segment in the allSegments vector.
public:
    enum class SegmentClassification {
        Invalid,
        LowCoverageStrand0,
        LowCoverageStrand1,
        LowCoverageUnclassified,
        HighCoverageStrand0,
        HighCoverageStrand1,
        HighCoverageAmbiguous,
        MaxValue
    };
    static string color(SegmentClassification);
    static string description(SegmentClassification);
private:
    vector<SegmentClassification> segmentClassifications;
    void classifySegments();

    void updateAssemblyGraph();
};
