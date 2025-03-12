#pragma once


// Shasta.
#include "AnchorPair.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {

    class AssemblyGraph;
    class AssemblyGraphVertex;
    class AssemblyGraphEdge;
    class AssemblyGraphStep;

    using AssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraphVertex,
        AssemblyGraphEdge>;

    class AnchorGraph;

}



class shasta::AssemblyGraphVertex {
public:
    AnchorId anchorId;

    AssemblyGraphVertex(AnchorId anchorId) : anchorId(anchorId) {}
};



class shasta::AssemblyGraphStep {
public:
    AnchorPair anchorPair;
    uint32_t averageOffset;
    uint32_t minOffset;
    uint32_t maxOffset;
};



class shasta::AssemblyGraphEdge : public vector<AssemblyGraphStep> {
public:
    uint64_t id = invalid<uint64_t>;
    uint32_t averageOffset = invalid<uint32_t>;
    uint32_t minOffset = invalid<uint32_t>;
    uint32_t maxOffset = invalid<uint32_t>;
    void computeOffsets();

    AnchorId anchorIdA() const
    {
        return (*this).front().anchorPair.anchorIdA;
    }
    AnchorId anchorIdB() const
    {
        return (*this).back().anchorPair.anchorIdB;
    }
};



class shasta::AssemblyGraph: public AssemblyGraphBaseClass {
public:

    AssemblyGraph(
        const Anchors&,
        const AnchorGraph&);

    void write(const string& name) const;
    void writeGfa(const string& fileName) const;
    void writeSegments(const string& fileName) const;
    void writeSegmentDetails(const string& fileName) const;

private:
    uint64_t nextEdgeId = 0;
};
