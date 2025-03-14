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

    // This is used for the BFS in AssemblyGraph::transitiveReduction.
    uint64_t bfsDistance = invalid<uint64_t>;
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

    // The length of an AssemblyGraphEdge is the estimated length of its sequence,
    // equal to averageOffset, which is the sum of the averageOffsets
    // of all the AssemblyGraphSteps of this edge.
    uint64_t length() const
    {
        return averageOffset;
    }

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


private:
    uint64_t nextEdgeId = 0;

    void transitiveReduction();
    void compress();

    void write(const string& name) const;
    void writeGfa(const string& fileName) const;
    void writeSegments(const string& fileName) const;
    void writeSegmentDetails(const string& fileName) const;
};
