#pragma once

// Local read following, one Tangle at a time.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    namespace ReadFollowing5 {

    class Graph;
    class Vertex;
    class Edge;

    using GraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge>;


    }

    class OrientedReadId;
    class Tangle;
}


// There is a vertex for each Segment internal to the Tangle,
// plus a vertex for each entrance and exit of the tangle.
// If a Segment is both an entrance and an exit
// there will be two vertices associated with it.
class shasta2::ReadFollowing5::Vertex {
public:
    Segment segment;
    bool isEntrance;
    bool isExit;

    vector<OrientedReadId> initialSupport;
    vector<OrientedReadId> finalSupport;

    Vertex(
        const AssemblyGraph&,
        Segment segment,
        bool isEntrance,
        bool isExit);
};



class shasta2::ReadFollowing5::Edge {
public:
    uint64_t commonCount;

    Edge(uint64_t commonCount) : commonCount(commonCount) {}

};



class shasta2::ReadFollowing5::Graph: public GraphBaseClass {
public:
    Graph(const AssemblyGraph&, uint64_t tangleId, const Tangle&);

private:
    const AssemblyGraph& assemblyGraph;
    const Tangle& tangle;
    uint64_t tangleId;
    bool isSelfComplementaryTangle;

    void createVertices();
    void createEdges();

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
};

