#pragma once


/******************************************************************

The StrandSeparator1 creates a single-stranded Anchors object
given a double-stranded one as input.

In the input double-stranded Anchors, anchors2, anchors must be
stored in consecutively numbered pairs of reverse complemented
anchors. This is not checked.

The output single-stranded Anchors, anchors1, must be empty
on input. This is not checked.

We use undirected bipartite graphs in which each vertex represents
an oriented read or an anchor.

******************************************************************/

// Shasta.
#include "Anchor.hpp"

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <random>



namespace shasta2 {
    class StrandSeparator1;
}


class shasta2::StrandSeparator1 {
public:


    StrandSeparator1(

        // The input double-stranded Anchors.
        // Anchors must be in consecutively numbered pairs of
        // reverse complemented anchors.
        const Anchors& anchors2,

        // The output single-stranded Anchors.
        // It must be empty on input.
        Anchors& anchors1
        );

private:
    const Anchors& anchors2;
    const uint64_t anchorCount;
    const uint64_t readCount;
    const uint64_t orientedReadCount;
    const uint64_t vertexCount;

    // A disjoint set data structure used ot represent bipartite graphs
    // in which the vertices are oriented reads or anchors.
    vector<uint64_t> rank;
    vector<uint64_t> parent;
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets;
    void initializeDisjointSets();



    // Each AnchorMarkerInfo in a positive anchor in anchors2
    // (that is, the AnchorId is even) generates a pair
    // of reverse complemented edges in the complete, double-stranded
    // bipartite graph.
    class EdgePair {
    public:
        AnchorId anchorId;
        ReadId readId;
        bool onSameStrand;
        EdgePair(AnchorId anchorId, ReadId readId, bool onSameStrand) :
            anchorId(anchorId),
            readId(readId),
            onSameStrand(onSameStrand)
        {
            SHASTA2_ASSERT((anchorId % 2) == 0); // anchorId must be a "positive" anchor.
        }

        AnchorId anchorIdA() const
        {
            return anchorId;
        }
        AnchorId anchorIdB() const
        {
            return anchorId ^ 1;
        }
        OrientedReadId orientedReadIdA() const
        {
            return OrientedReadId(readId, onSameStrand ? 0 : 1);
        }
        OrientedReadId orientedReadIdB() const
        {
            return OrientedReadId(readId, onSameStrand ? 1 : 0);
        }
    };
    vector<EdgePair> edgePairs;
    void gatherEdgePairs();

    // Get the vertexIds for a given EdgePair.
    class VertexIds {
        public:
        uint64_t anchorVertexA;
        uint64_t anchorVertexB;
        uint64_t orientedReadVertexA;
        uint64_t orientedReadVertexB;
    };
    VertexIds getVertexIds(const EdgePair& edgePair) const
    {
        const AnchorId anchorIdA = edgePair.anchorIdA();
        const AnchorId anchorIdB = edgePair.anchorIdB();
        const OrientedReadId orientedReadIdA = edgePair.orientedReadIdA();
        const OrientedReadId orientedReadIdB = edgePair.orientedReadIdB();

        VertexIds vertexIds;
        vertexIds.anchorVertexA = anchorIdA;
        vertexIds.anchorVertexB = anchorIdB;
        vertexIds.orientedReadVertexA = anchorCount + orientedReadIdA.getValue();
        vertexIds.orientedReadVertexB = anchorCount + orientedReadIdB.getValue();

        return vertexIds;
    }

    // Apply the given EdgePair to the disjointSets.
    void applyEdgePairToDisjointSets(const EdgePair&);


    // A connected component of the double-stranded or single-stranded
    // bipartite graph. The AnchorIds and OrientedReadIds are stored sorted.
    // We don't store trivial components with one vertex.
    class Component {
    public:
        vector<AnchorId> anchorIds;
        vector<OrientedReadId> orientedReadIds;

        bool isTrivial() const
        {
            return anchorIds.size() + orientedReadIds.size() <= 1;
        }

        void check() const;

        bool isDoubleStranded() const
        {
            return orientedReadIds[0].getReadId() == orientedReadIds[1].getReadId();
        }
    };



    // A set of connected components of the double-stranded or single-stranded
    // bipartite graph.
    class Components {
    public:

        vector<Component> components;

        // The componentId that each AnchorId belongs to,
        // or invalid<uint64_t> if none.
        // Indexed by AnchorId.
        vector<uint64_t> anchorComponent;

        // The componentId that each OrientedReadId belongs to,
        // or invalid<uint64_t> if none.
        // indexed by OrientedReadId::getValue().
        vector<uint64_t> orientedReadComponent;

        void initialize(uint64_t anchorCount, uint64_t orientedReadCount);
        void fill();
    };

    // Construct a Components object from the current disjointSets.
    void createComponents(Components&);



    // The connected components of the double-stranded bipartite graph
    // that uses all the edge pairs.
    Components components2;
    void computeComponents2();

    // Generate single-strand components from double-strand components.
    Components components1;
    void createComponents1();
    void createComponent1(uint64_t componentId2);

    // Random generator used for shuffles.
    std::mt19937 random;

    void generateSingleStrandedAnchors(Anchors& anchors1) const;

};

