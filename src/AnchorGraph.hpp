#pragma once

// In the AnchorGraph, each vertex corresponds to an AnchorId
// and each edge corresponds to an AnchorPair.
// It uses boost::vecS as its second template argument,
// and as a result its vertex descriptors are AnchorIds.

// Shasta.
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"



namespace shasta2 {

    class AnchorGraph;
    class AnchorGraphEdge;
    using AnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        boost::no_property,
        AnchorGraphEdge>;

    class Anchors;
    class Journeys;
    class ReadLengthDistribution;
}



class shasta2::AnchorGraphEdge {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    uint64_t id = invalid<uint64_t>;
    bool useForAssembly = false;

    AnchorGraphEdge(const AnchorPair& anchorPair, uint64_t offset, uint64_t id) :
        anchorPair(anchorPair),
        offset(offset),
        id(id)
    {}

    AnchorGraphEdge() {}

    uint64_t coverage() const {return anchorPair.size();}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & id;
        ar & useForAssembly;
    }
};



class shasta2::AnchorGraph :
    public AnchorGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AnchorGraph> {
public:

    // Construct the AnchorGraph from the Journeys.
    // Only include edges with at least the specified minCoverage.
    AnchorGraph(const Anchors&, const Journeys&, uint64_t minEdgeCoverage);

    // Constructor to create an anchor similarity graph,
    // in which an edge between two anchors is created if the
    // oriented read compositions of the two anchors are
    // sufficiently similar.
    class UseSimilarity {};
    AnchorGraph(const Anchors&, const Journeys&, uint64_t minEdgeCoverage, const UseSimilarity&);



    // This uses read following in the complete AnchorGraph
    // to create the AnchorGraph to be used for assembly.
    // This is meant to be used with strict anchor generation,
    // where most anchors correspond to a single copy.
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        const AnchorGraph& completeAnchorGraph);

    // Types and functions used by the above constructor.

    // A Subgraph of the AnchorGraph, created using a BFS starting
    // at anchorIdStart  and moving in the specified direction.
    // It skips vertices that have less the minCommon OrientedReadIds
    // in common with anchorIdStart.
    class SubgraphVertex {
    public:
        AnchorId anchorId = invalid<AnchorId>;
        uint64_t commonCount = invalid<uint64_t>;
        SubgraphVertex(
            AnchorId anchorId,
            uint64_t commonCount) :
            anchorId(anchorId),
            commonCount(commonCount)
        {}
        SubgraphVertex() {}
        // These are used by approximate topological sort.
        uint64_t color = invalid<uint64_t>;
        uint64_t rank = invalid<uint64_t>;
    };
    class SubgraphEdge {
    public:
        uint64_t coverage;
        bool isDagEdge = false; // From approximate topological sort.
        SubgraphEdge(uint64_t coverage, bool isDagEdge=false) : coverage(coverage), isDagEdge(isDagEdge) {}
    };

    using SubgraphBaseClass = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS,
        SubgraphVertex, SubgraphEdge>;
    class Subgraph: public SubgraphBaseClass {
    public:
        Subgraph(
            const Anchors&,
            const AnchorGraph& anchorGraph,
            AnchorId,
            uint64_t direction,
            uint64_t minCommonCount);
        vertex_descriptor vStart;
        uint64_t direction;
        std::map<AnchorId, vertex_descriptor> vertexMap;

        // This does an approximate topological sort, then removes
        // edges not flagged as DAG edges.
        void removeCycles();

        void transitiveReduction();

        void writeFasta() const;

        // Recursively prune leafs with commonCount less than minLeafCommonCount.
        void prune(uint64_t);

        // Find exits, ignoring isolated verties.
        void findExits(vector<vertex_descriptor>&) const;

        // If there are multiple exits, keep only vertices that are
        // reachable (backward) from all exits.
        void pruneMultipleExits();

        // Get the AnchorIds of the linear chain on vertices containing vStart.
        // They are returned in order.
        void getLinearPortion(vector<AnchorId>&) const;

        // Constructor the dominator tree (in the given direction).
        class DominatorTree{};
        Subgraph(const Subgraph&, const DominatorTree&, const Anchors&);

        // Walk up the dominator tree.
        // Note this returns a path in the dominator tree.
        void walkUp(const Subgraph& dominatorTree, vector<vertex_descriptor>& path) const;

        // Walk up the dominator tree.
        // This will assert if not called on the dominator tree.
        void walkUp(vertex_descriptor, vector<vertex_descriptor>& path) const;

        void writeGraphviz(const string& fileName, const Anchors&) const;
        void writeGraphviz(ostream&, const Anchors&) const;
        void writeHtml(ostream&, const Anchors&) const;

        void writeFasta(
            const vector<vertex_descriptor>& path,
            const string& fileName,
            const Anchors&) const;
        void writeFasta(
            const vector<vertex_descriptor>& path,
            ostream& fasta,
            const Anchors&) const;
        void writeFastaHtml(
            const vector<vertex_descriptor>& path,
            ostream& html,
            const Anchors&) const;
    };




public:
    // Constructor from binary data.
    AnchorGraph(const MappedMemoryOwner&, const string& name);

    uint64_t nextEdgeId = 0;

    void transitiveReduction(
        uint64_t transitiveReductionMaxEdgeCoverage,
        uint64_t maxDistance);
private:
    bool transitiveReductionCanRemove(edge_descriptor, uint64_t transitiveReductionMaxDistance) const;
public:

    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorGraphBaseClass>(*this);
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);
};

