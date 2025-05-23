#pragma once

// In the AnchorGraph, each vertex corresponds to an AnchorId
// and each edge corresponds to an AnchorPair.
// It uses boost::vecS as its second template argument,
// and as a result its vertex descriptors are AnchorIds.
// The AnchorGraph is not stored. It is used to create
// the initial Assembly graph and then discarded.

// Shasta.
#include "AnchorPair.hpp"
#include "MappedMemoryOwner.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

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
    }
};



class shasta::AnchorGraph :
    public AnchorGraphBaseClass,
    public MappedMemoryOwner {
public:
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        uint64_t minEdgeCoverage);

#if 0
    // Constructor that splits edges that have an AnchorPair
    // with inconsistent offsets.
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        double aDrift,
        double bDrift);
#endif

    // Constructor that splits edges that have an AnchorPair
    // with inconsistent offsets, and also does local search to
    // eliminate dead ends where possible.
    struct FixDeadEnds{};
    AnchorGraph(
        const Anchors&,
        const Journeys&,
        const ReadLengthDistribution&,
        uint64_t minEdgeCoverageNear,
        uint64_t minEdgeCoverageFar,
        double aDrift,
        double bDrift,
        FixDeadEnds);

    // Constructor from binary data.
    AnchorGraph(const MappedMemoryOwner&);


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
        AnchorPair&,
        uint64_t& offset
        ) const;



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

