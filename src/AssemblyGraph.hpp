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
};



class shasta::AssemblyGraphEdge : public vector<AssemblyGraphStep> {
public:
};



class shasta::AssemblyGraph: public AssemblyGraphBaseClass {
public:

    AssemblyGraph(const AnchorGraph&);
};
