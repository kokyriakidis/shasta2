#pragma once

// Local read following, one Tangle at a time.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"

// Standard library.
#include "fstream.hpp"
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
    double logP;
    double weight;

    Edge(uint64_t commonCount, double logP) :
        commonCount(commonCount),
        logP(logP),
        weight(std::pow(10., -0.1 * logP))
    {}

};



class shasta2::ReadFollowing5::Graph: public GraphBaseClass {
public:
    Graph(const AssemblyGraph&, uint64_t tangleId, const Tangle&, ostream& html);

private:
    const AssemblyGraph& assemblyGraph;
    const Tangle& tangle;
    uint64_t tangleId;
    bool isSelfComplementaryTangle;

    vector<vertex_descriptor> entranceVertices;
    vector<vertex_descriptor> exitVertices;

    void createVertices();
    void createEdges();

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;



    // Shortest paths on the graph between each entrance/exit pair.
    class ShortestPath {
    public:
        Segment entrance;
        Segment exit;
        double distance;

        // The path Segments, excluding the entrance and the exit.
        vector<Segment> segments;
    };
    vector< vector<ShortestPath> > shortestPaths;
    void findShortestPaths();
    void writeShortestPaths(ostream& html);



};

