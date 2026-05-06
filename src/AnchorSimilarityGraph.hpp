#pragma once

// The AnchorSimilarityGraph is a directed graph in which
// each vertex represents an anchor. An edge anchorId0->anchorId1
// is created if the oriented read compositions of
// anchorId0 and anchorId1 are sufficiently similar.
// The connectivity of this graph is very high,
// so this works better when --min-anchor-coverage is
// low and most anchors correspond to a single copy of sequence.

// Shasta2.
#include "Anchor.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"



namespace shasta2 {

    class AnchorSimilarityGraphEdge;
    class AnchorSimilarityGraph;
    using AnchorSimilarityGraphBaseClass = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::bidirectionalS,
        boost::no_property,
        AnchorSimilarityGraphEdge>;

    class AnchorGraph;
}



class shasta2::AnchorSimilarityGraphEdge {
public:

    // The weight is computed by analyzing the AnchorPair for this edge.
    // It is low if the two anchors have very similar read composition.
    double weight = invalid<double>;

    uint64_t baseOffset = invalid<uint64_t>;
    bool isShortestPathEdge = false;

    AnchorSimilarityGraphEdge(
        double weight = invalid<double>,
        uint64_t baseOffset = invalid<uint64_t>) :
        weight(weight),
        baseOffset(baseOffset)
    {}

    bool operator==(const AnchorSimilarityGraphEdge& that) const
    {
        return
            (weight == that.weight) and
            (baseOffset == that.baseOffset) and
            (isShortestPathEdge == that.isShortestPathEdge);
    }

    template<class Archive> void serialize(
        Archive& ar,
        [[maybe_unused]] unsigned int version)
    {
        ar & weight;
        ar & baseOffset;
        ar & isShortestPathEdge;
    }
};



class shasta2::AnchorSimilarityGraph :
    public AnchorSimilarityGraphBaseClass,
    public MappedMemoryOwner {
public:

    // Construct the AnchorSimilarityGraph from the completeAnchorGraph.
    // Only include edges with at least the specified minCoverage.
    AnchorSimilarityGraph(
        const Anchors&,
        const AnchorGraph& completeAnchorGraph);

    // Constructor from binary data.
    AnchorSimilarityGraph(const MappedMemoryOwner&, const string& name);



    // Work areas for the computation of shortest path trees.
    // This implements property maps that can be used in the call to
    // dijkstra_shortest_paths_no_init. This gives better memory
    // locality than using property maps created using
    // using boost::make_iterator_property_map.
    class ShortestPathTreeWorkAreas {
    public:

        class Data {
        public:
            AnchorId predecessor;
            double distance = std::numeric_limits<double>::max();
            boost::default_color_type color = boost::default_color_type::white_color;
        };
        vector<Data> data;

        // After a call to dijkstra_shortest_paths_no_init
        // using the above DijkstraVisitor, the accessibleVertices
        // vector contains the vertices that were reached and for
        // which the data values need to be reset.
        vector<AnchorId> accessibleVertices;

        ShortestPathTreeWorkAreas(uint64_t anchorCount);
        void reset();
        void check() const;



        // Property maps.
        class PropertyMap {
        public:
            ShortestPathTreeWorkAreas& workAreas;
            PropertyMap(ShortestPathTreeWorkAreas& workAreas) :
                workAreas(workAreas) {}
        };
        class PredecessorMap : public PropertyMap {
        public:
            PredecessorMap(ShortestPathTreeWorkAreas& workAreas) :
                PropertyMap(workAreas) {}
        };
        class DistanceMap : public PropertyMap {
        public:
            DistanceMap(ShortestPathTreeWorkAreas& workAreas) :
                PropertyMap(workAreas) {}
            using value_type = double;
        };
        class ColorMap : public PropertyMap {
        public:
            ColorMap(ShortestPathTreeWorkAreas& workAreas) :
                PropertyMap(workAreas) {}
            using value_type = boost::default_color_type;
        };
        PredecessorMap predecessorMap = *this;
        DistanceMap distanceMap = *this;
        ColorMap colorMap = *this;
    };



    // Compute a shortest path tree starting at the given AnchorId.
    // These are forward paths. Backward paths can be created
    // by reverse complementing.
    // On input, predecessorMap, distanceMap, and colorMap
    // must be vectors of size anchors.size() set as follows:
    // - predecessorMap[anchorId] == anchorId
    // - distanceMap[anchorId] == std::numeric_limits<double>::max()
    // - colorMap[anchorId] == boost::default_color_type::white_color
    // For performance, these conditions are not checked.
    // On exit, the accessibleVertices vector contains the AnchorIds
    // in the shortest path tree, that is, all the AnchorIds
    // accessible from the root anchorId.
    // In addition, on exit, the ShortestPathTreeWorkAreas
    // describe the shortest path tree with root at the given AnchorId.
    // Only the entries corresponding to the accessibleVertices are
    // changed. The caller should call ShortestPathTreeWorkAreas::reset
    // before reusing the ShortestPathTreeWorkAreas.
    void createShortestPathTree(
        AnchorId,
        ShortestPathTreeWorkAreas&
        ) const;
    void createShortestPathTree(AnchorId, const Anchors&) const;
    void flagShortestPathEdges(const Anchors&);

    void checkStrandInvariant() const;


    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorSimilarityGraphBaseClass>(*this);
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);

private:

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = 6;
    const double a = 3.;
    const double b = 10.;
    const double minLogP = 10.;
    const uint64_t pruneLength = 10;


    void createVertices(const Anchors&);

    // Create all the edges.
    void createEdges(const Anchors&, const AnchorGraph& completeAnchorGraph);

    // Create the edges with source anchorIdA.
    void createEdges(
        const Anchors&,
        const AnchorGraph& completeAnchorGraph,
        AnchorId anchorIdA,
        vector<uint8_t>& color);

    // Graphviz output only includes the edges flagged as shortest path edges.
    void writeGraphviz(const string& fileName, bool shortPathEdgesOnly) const;
    void writeGraphviz(ostream&, bool shortPathEdgesOnly) const;



    // Classes to describe a shortest path tree rooted at a given AnchorId.

    class ShortestPathTreeVertex {
    public:
        AnchorId anchorId;

        // Distance (number of edges) to the root and to the most distant leaf down
        // from this vertex.
        uint64_t rank = invalid<uint64_t>;
        uint64_t longestDistanceToLeaf = invalid<uint64_t>;

        // Distance to the root, using estimated base offsets.
        uint64_t offset =  invalid<uint64_t>;

        ShortestPathTreeVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}
    };

    class ShortestPathTreeEdge {
    public:
        double logP;
        // The rest is filled in by fillInEdgeInformation().
        uint64_t common = invalid<uint64_t>;
        uint64_t missing = invalid<uint64_t>;
        uint64_t offset = invalid<uint64_t>;
        uint64_t minOffset = invalid<uint64_t>;
        uint64_t maxOffset = invalid<uint64_t>;
        ShortestPathTreeEdge(double logP) : logP(logP) {}
    };

    using ShortestPathTreeBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        ShortestPathTreeVertex,
        ShortestPathTreeEdge>;

    class ShortestPathTree : public ShortestPathTreeBaseClass {
    public:
        ShortestPathTree(
            const AnchorSimilarityGraph&,
            AnchorId rootAnchorId,
            const ShortestPathTreeWorkAreas&);
        AnchorId rootAnchorId;
        std::map<AnchorId, vertex_descriptor> vertexMap;

        uint64_t maximumRank() const;
        uint64_t maximumPathLength() const;

        // This removes vertices with longestDistanceToLeaf < pruneLength
        // and their descendants, as long as they have a sibling
        // with greater longestDistanceToLeaf.
        void prune(uint64_t pruneLength);

        // Store additionalInformation in the vertices and edges.
        // The one for the edges must be called first.
        void fillEdgeInformation(const Anchors&);
        void fillVertexInformation();

        // Find the sequence of vertices or AnchorIds
        // of a path starting at root and ending at the
        // specified vertex or AnchorId.
        void findPath(vertex_descriptor, vector<vertex_descriptor>&) const;
        void findPathAnchorIds(AnchorId, vector<AnchorId>&) const;

        void writeGraphviz(const string& fileName) const;
        void writeGraphviz(ostream&) const;
    private:
        void computeRanks();
        void computeLongestDistancesToLeaf();
        void gatherVerticesByRank(vector< vector<vertex_descriptor> >&) const;
    };

    // Graphviz output of shortest path edges (only) of the AnchorSimilarityGraph,
    // highlighting a given ShortestPathTree and a path on the ShortestPathTree.
    void writeGraphviz(
        const string& fileName,
        const ShortestPathTree&,
        const vector<ShortestPathTree::vertex_descriptor>& path) const;
    void writeGraphviz(
        ostream&,
        const ShortestPathTree&,
        const vector<ShortestPathTree::vertex_descriptor>& path) const;

};



// Get/put functions and boost::property_traits specializations
// for ShortestPathTreeWorkAreas property maps.
// There is probably a better way to do this.
namespace shasta2 {
    inline
        AnchorId get(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::PredecessorMap& predecessorMap,
        AnchorId anchorId)
    {
        return predecessorMap.workAreas.data[anchorId].predecessor;
    }
    inline
        double get(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::DistanceMap& distanceMap,
        AnchorId anchorId)
    {
        return distanceMap.workAreas.data[anchorId].distance;
    }
    inline
        boost::default_color_type get(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::ColorMap& colorMap,
        AnchorId anchorId)
    {
        return colorMap.workAreas.data[anchorId].color;
    }
    inline
        void put(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::PredecessorMap& predecessorMap,
        AnchorId anchorId,
        AnchorId predecessor)
    {
        predecessorMap.workAreas.data[anchorId].predecessor = predecessor;
    }
    inline
        void put(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::DistanceMap& distanceMap,
        AnchorId anchorId,
        double distance)
    {
        distanceMap.workAreas.data[anchorId].distance = distance;
    }
    inline
        void put(
        AnchorSimilarityGraph::ShortestPathTreeWorkAreas::ColorMap& colorMap,
        AnchorId anchorId,
        boost::default_color_type color)
    {
        colorMap.workAreas.data[anchorId].color = color;
    }
}
namespace boost {
    template<> struct property_traits<shasta2::AnchorSimilarityGraph::ShortestPathTreeWorkAreas::PredecessorMap> {
        typedef shasta2::AnchorId key_type;
        typedef shasta2::AnchorId value_type;
        typedef shasta2::AnchorId reference;
        typedef boost::lvalue_property_map_tag category;
    };
    template<> struct property_traits<shasta2::AnchorSimilarityGraph::ShortestPathTreeWorkAreas::DistanceMap> {
        typedef shasta2::AnchorId key_type;
        typedef double value_type;
        typedef double reference;
        typedef boost::lvalue_property_map_tag category;
    };
    template<> struct property_traits<shasta2::AnchorSimilarityGraph::ShortestPathTreeWorkAreas::ColorMap> {
        typedef shasta2::AnchorId key_type;
        typedef boost::default_color_type value_type;
        typedef boost::default_color_type reference;
        typedef boost::lvalue_property_map_tag category;
    };
}
