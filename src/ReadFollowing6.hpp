#pragma once

// ReadFollowing6 code is patterned after StrandSeparation1.
// It uses a bipartite graph.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"
#include "tuple.hpp"
#include "vector.hpp"

namespace shasta2 {
    namespace ReadFollowing6 {
        class ReadOccurrence;
        class SegmentInfo;
        class Tangle;

        class BipartiteGraphVertex;
        class BipartiteGraphEdge;
        using BipartiteGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::undirectedS,
            BipartiteGraphVertex,
            BipartiteGraphEdge>;
        class BipartiteGraph;
    }
}



// Information about a single Segment in a Tangle.
class shasta2::ReadFollowing6::SegmentInfo {
public:
    Segment segment;
    uint64_t id;
    bool isEntrance = false;
    bool isExit = false;
    uint64_t length;
    double coverage;
    BipartiteGraphBaseClass::vertex_descriptor v = BipartiteGraphBaseClass::null_vertex();
};



// Class used to count occurrences of reads in Segments.
class shasta2::ReadFollowing6::ReadOccurrence {
public:
    uint64_t segmentIndex = invalid<uint64_t>;
    uint64_t frequency = 0;

    ReadOccurrence() {}
    ReadOccurrence(uint64_t segmentIndex) :
        segmentIndex(segmentIndex) {}

    bool operator==(const ReadOccurrence& that) const
    {
        return segmentIndex == that.segmentIndex;
    }

    bool operator<(const ReadOccurrence& that) const
    {
        return segmentIndex < that.segmentIndex;
    }
};



// A vertex of the BipartiteGraph can represent an OrientedReadId
// or a Segment.
class shasta2::ReadFollowing6::BipartiteGraphVertex {
public:
    bool isSegment;

    OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);

    uint64_t segmentIndex = invalid<uint64_t>;

    BipartiteGraphVertex(uint64_t segmentIndex) :
        isSegment(true),
        segmentIndex(segmentIndex)
    {}

    BipartiteGraphVertex(OrientedReadId orientedReadId) :
        isSegment(false),
        orientedReadId(orientedReadId)
    {}

    BipartiteGraphVertex()
    {}
};



class shasta2::ReadFollowing6::BipartiteGraphEdge {
public:

    // The number of times (steps) this OrientedReadId appears
    // in this Segment.
    uint64_t frequency;

    // The weight is the ratio of frequency to the total frequency
    // for the same OrientedReadId.
    double weight;
};



class shasta2::ReadFollowing6::BipartiteGraph : public BipartiteGraphBaseClass {
public:

    // Add a vertex representing a Segment.
    vertex_descriptor addVertex(uint64_t segmentIndex);

    // Add a vertex representing an OrientedReadId.
    vertex_descriptor addVertex(OrientedReadId);

};



// Note this is unrelated to shasta2::Tangle.
class shasta2::ReadFollowing6::Tangle {
public:

    // Tangle constructor.
    // The tangleVertices must be sorted by id.
    // The Tangle must not be self-complementary - that is,
    // if it includes a vertex, it cannot also include its reverse complement.
    // The code will work on the tangle passed in and also make the
    // necessary changes on its reverse complement to keep the AssemblyGraph
    // strand-symmetric.
    // The debugOutputBaseName and tangleId are only used for debug output.
    Tangle(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
        ostream& html,
        uint64_t tangleId,
        const string& debugOutputBaseName);

private:

    // Constructor arguments.
    AssemblyGraph& assemblyGraph;
    const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices;
    ostream& html;
    uint64_t tangleId;
    string debugOutputBaseName;

    uint64_t id(AssemblyGraphBaseClass::vertex_descriptor) const;
    uint64_t id(Segment) const;

    // The tangleVertices must be sorted by id.
    // If the Tangle includes a vertex, it cannot also include its reverse complement.
    void checkTangleVertices();

    // Segments in the Tangle.
    // This includes Segments internal to the Tangle plus
    // entrances and exits.
    // Ordered by the id.
    vector<SegmentInfo> segmentInfos;
    void gatherSegments();
    void writeSegments();

    // Count occurrences of reads in the Segments.
    std::map<OrientedReadId, vector<ReadOccurrence> > readOccurrenceMap;
    void countReadOccurrences();
    void writeReadOccurrences();

    // The BipartiteGraph has a vertex for Segment
    // and a vertex for each OrientedReadId that appears
    // in more than one Segment.
    // An edge (Segment--OrientedReadId) contains the number of times (steps)
    // that OrientedReadId occures in that Segment.
    // The vertex corresponding to a Segment is stored in the SegmentInfo.
    // The vertex corresponding to an OrientedReadId is stored in the
    // orientedReadIdVertexMap.
    BipartiteGraph bipartiteGraph;
    std::map<OrientedReadId, BipartiteGraph::vertex_descriptor> orientedReadIdVertexMap;
    void createBipartiteGraph();
    void writeBipartiteGraph();

};
