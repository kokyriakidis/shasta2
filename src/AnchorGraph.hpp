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

