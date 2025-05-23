#pragma once

// In the AnchorGraph, each vertex corresponds to an AnchorId
// and each edge corresponds to an AnchorPair.
// It uses boost::vecS as its second template argument,
// and as a result its vertex descriptors are AnchorIds.
// The AnchorGraph is not stored. It is used to create
// the initial Assembly graph and then discarded.

// Shasta.
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <utility.hpp>
#include <vector.hpp>



namespace shasta {

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



class shasta::AnchorGraphEdge {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    uint64_t id = invalid<uint64_t>;
    bool isParallelEdge = false;
    bool useForAssembly = false;
    bool addedAtDeadEnd = false;

    AnchorGraphEdge(const AnchorPair& anchorPair, uint64_t offset, uint64_t id) :
        anchorPair(anchorPair),
        offset(offset),
        id(id)
    {}

    AnchorGraphEdge() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & id;
        ar & isParallelEdge;
        ar & useForAssembly;
        ar & addedAtDeadEnd;
    }
};



class shasta::AnchorGraph :
    public AnchorGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AnchorGraph> {
public:
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        uint64_t minEdgeCoverage);

    // Constructor that splits edges that have an AnchorPair
    // with inconsistent offsets, and also does local search to
    // eliminate dead ends where possible.
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        const ReadLengthDistribution&,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        double aDrift,
        double bDrift,
        uint64_t threadCount);

    // Constructor from binary data.
    AnchorGraph(const MappedMemoryOwner&);

    uint64_t nextEdgeId = 0;

    // Dijkstra search.
    // This performs a shortest path search starting at the specified AnchorId
    // and stops when it finds a consistent AnchorPair that
    // satisfies the coverage criteria and can therefore be used to create a new edge.
    // If successful, this returns true and:
    // - If direction is 0, the last arguments is set to an AnchorPair with AnchorIdA = startAnchorId.
    // - If direction is 1, the last arguments is set to an AnchorPair with AnchorIdB = startAnchorId.
    // If not successful, this return false and the last argument is not modified.
    bool search(
        uint64_t direction,
        AnchorId startAnchorId,
        const Anchors&,
        const ReadLengthDistribution&,
        double aDrift,
        double bDrift,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        uint64_t maxDistance,
        AnchorPair&,
        uint64_t& offset
        ) const;
    bool searchForward(
        AnchorId startAnchorId,
        const Anchors&,
        const ReadLengthDistribution&,
        double aDrift,
        double bDrift,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        uint64_t maxDistance,
        AnchorPair&,
        uint64_t& offset
        ) const;
    bool searchBackward(
        AnchorId startAnchorId,
        const Anchors&,
        const ReadLengthDistribution&,
        double aDrift,
        double bDrift,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        uint64_t maxDistance,
        AnchorPair&,
        uint64_t& offset
        ) const;



    // Eliminate dead ends where possible, using shortest path searches.
    void handleDeadEnds(
        const Anchors&,
        const ReadLengthDistribution&,
        double aDrift,
        double bDrift,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        uint64_t maxDistance,
        uint64_t threadCount);

    class HandleDeadEndsData {
    public:

        // The arguments of handleDeadEnds, so they are visible to the threads.
        const Anchors* anchors;
        const ReadLengthDistribution* readLengthDistribution;
        double aDrift;
        double bDrift;
        uint64_t minEdgeCoverageNear;
        uint64_t minEdgeCoverageFar;
        uint64_t maxDistance;

        // The list of the dead end AnchorIds and their direction
        // (0 = forward dead end, 1 = backward dead end).
        vector< pair<AnchorId, uint64_t> > deadEnds;

        // The AnchorPairs found by each thread and their offset.
        vector< vector< pair<AnchorPair, uint64_t> > > threadPairs;
    };
    HandleDeadEndsData handleDeadEndsData;

    void handleDeadEndsThreadFunction(uint64_t threadId);



    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorGraphBaseClass>(*this);
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save() const;
    void load();
};

