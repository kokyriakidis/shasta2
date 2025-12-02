#pragma once



/******************************************************************

The Strand separator creates a single-stranded Anchors object
given a double-stranded one as input.

In the input double-stranded Anchors, anchors2,
anchors must be in consecutively numbered pairs of
reverse complemented anchors.

The output single-stranded Anchors, anchors1,
must be empty on input.

We use an undirected bipartite graph in which each vertex represents
an OrientedReadId or an AnchorId. The graph only contains edges
between a vertex corresponding to an OrientedReadId
and a vertex corresponding to an AnchorId.

We use two version of this graph:
- A double-stranded version, in which each AnchorMarkerInfo in the input
double-stranded anchors generates an edge.
- A single-stranded version which includes a subset of the edges.

We use "2" suffixes to refer to the double-stranded bipartite graph
and "1" suffixes to refer to the single-stranded bipartite graph.

******************************************************************/

#include "Anchor.hpp"

namespace shasta2 {
    class StrandSeparator;
}



class shasta2::StrandSeparator {
public:


    StrandSeparator(

        // The input double-stranded Anchors.
        // Anchors must be in consecutively numbered pairs of
        // reverse complemented anchors.
        const Anchors& anchors2,

        // The output single-stranded Anchors.
        // It must be empty on input.
        Anchors& anchors1
        );

private:
    const uint64_t anchorCount;
    const uint64_t readCount;
    const uint64_t orientedReadCount;

    // Each AnchorMarkerInfo in anchors2 generates an edge
    // of the double-stranded bipartite graph.
    // The edge has a keep flag which is set to true for
    // edges that are kept in the single-stranded bipartite graph.
    class Edge {
    public:
        AnchorId anchorId;
        OrientedReadId orientedReadId;
        bool keep = false;
        Edge(AnchorId anchorId, OrientedReadId orientedReadId) :
            anchorId(anchorId), orientedReadId(orientedReadId) {}
    };
    vector<Edge> edges;
    void gatherEdges(const Anchors& anchors2);



    // Compute connected components of the double-stranded
    // bipartite graph. They are sorted by decreasing number
    // of oriented reads. Components with just one
    // oriented read are not stored.
    class Component {
    public:
        vector<AnchorId> anchorIds;
        vector<OrientedReadId> orientedReadIds;
    };
    vector<Component> components2;
    void computeComponents2();

    // The index of the component that each AnchorId belongs to.
    // Indexed by the AnchorId.
    vector<uint64_t> anchorComponent2;

    // The index of the component that each OrientedReadId belongs to.
    // Indexed by OrientedReadId::getValue()
    vector<uint64_t> orientedReadComponent2;

    // Strand separation.
    void separate(const Anchors& anchors2);
};
