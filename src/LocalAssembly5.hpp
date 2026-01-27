#pragma once

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "Kmer.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "vector.hpp"



namespace shasta2 {
    class LocalAssembly5;
    class LocalAssembly5Vertex;
    class LocalAssembly5Edge;

    using LocalAssembly5BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalAssembly5Vertex,
        LocalAssembly5Edge>;

    class Anchors;
    class AnchorPair;
    class Base;
    class Marker;
    class MarkerInfo;
}



class shasta2::LocalAssembly5Vertex {
public:

    // The id of the Kmer that corresponds to this vertex.
    // This is the index of this Kmer in the LocalAssembly5::kmers vector.
    uint64_t kmerId;

    class Info {
    public:

        // Index in orientedReadInfos.
        uint64_t orientedReadIndex;

        // Ordinal.
        uint32_t ordinal;
    };
    vector<Info> infos;

    LocalAssembly5BaseClass::vertex_descriptor dominator = LocalAssembly5BaseClass::null_vertex();
    bool isOnDominatorTreePath = false;

    // These are used for approximate topological ordering.
    // The color is also used by computeAssemblyPath.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    LocalAssembly5Vertex(uint64_t kmerId = invalid<uint64_t>) :
        kmerId(kmerId) {}
};



class shasta2::LocalAssembly5Edge {
public:
    class Info {
    public:

        // Index in orientedReadInfos.
        uint64_t orientedReadIndex;

        // Ordinal.
        uint32_t ordinal0;

        // Ordinal.
        uint32_t ordinal1;
    };
    vector<Info> infos;

    // This is used for approximate topological ordering.
    // which is only used by writeGraphviz.
    bool isDagEdge = false;
};



class shasta2::LocalAssembly5 : public LocalAssembly5BaseClass {
public:
    // Hide boost::adjacency_list::Base;
    using Base = shasta2::Base;

    // This assembles between anchorIdA and anchorIdB of the given AnchorPair.
    // This can use all of the OrientedReadIds in the AnchorPair
    // and/or in the additionalOrientedReadIds that also appear in the
    // Kmers corresponding to the left OR right Anchors.
    // In contrast, LocalAssembly4 uses oriented reads that appear in the
    // Kmers corresponding to the left AND right Anchors.
    // At sufficiently high coverage, LocalAssembly4 is a better choice
    // because it is simpler and faster. LocalAssembly5 is better
    // at low coverage, when it is important to use as many reads as possible,
    // The additionalOrientedReadIds must be sorted.
    // The additionalOrientedReadIds are allowed to contain OrientedReadIds
    // that are also in the AnchorPair.
    LocalAssembly5(
        const Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

    // Assembled sequence and its coverage.
    vector<Base> sequence;
    vector<uint64_t> coverage;

private:

    // EXPOSE WHEN CODE STABILIZES.
    const double drift = 0.1;

    const Anchors& anchors;
    ostream& html;

    // The two anchors of the AnchorPair used for this assembly.
    AnchorId leftAnchorId;
    AnchorId rightAnchorId;

    // The corresponding Kmers.
    Kmer leftKmer;
    Kmer rightKmer;

    // MarkerInfos for the Marker Kmers corresponding to the left
    // and right anchors of the AnchorPair being assembled.
    vector<MarkerInfo> leftMarkerInfos;
    vector<MarkerInfo> rightMarkerInfos;
    void fillMarkerInfos();

    // The union of the OrientedReadIds in the AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> allOrientedReadIds;
    void gatherAllOrientedReads(
        const AnchorPair&,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // The oriented reads used in this local assembly.
    // These are all OrientedReadIds that are in allOrientedReadIds and
    // also on the left or right Marker.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // Whether this OrientedReadId appear in the left and right Markers.
        bool isOnLeftMarker = false;
        bool isOnRightMarker = false;
        bool isOnBothMarkers() const
        {
            return isOnLeftMarker and isOnRightMarker;
        }

        // The ordinals of this OrientedReadId in the left/right Markers, if any.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            SHASTA2_ASSERT(isOnLeftMarker);
            SHASTA2_ASSERT(isOnRightMarker);
            SHASTA2_ASSERT(rightOrdinal > leftOrdinal);
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right Markers, if any.
        // These are positions in the oriented read sequence
        // of the mid point of the marker.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            SHASTA2_ASSERT(isOnLeftMarker);
            SHASTA2_ASSERT(isOnRightMarker);
            SHASTA2_ASSERT(rightPosition > leftPosition);
            return rightPosition - leftPosition;
        }



        // The region of this oriented read that can be used in this local assembly.
        // The positions are positions of the marker midpoints.
        class LocalRegion {
        public:
            uint32_t firstOrdinal = invalid<uint32_t>;
            uint32_t lastOrdinal = invalid<uint32_t>;
            uint32_t firstPosition = invalid<uint32_t>;
            uint32_t lastPosition = invalid<uint32_t>;

            uint32_t ordinalOffset() const
            {
                return lastOrdinal - firstOrdinal;
            }

            uint32_t positionOffset() const
            {
                return lastPosition - firstPosition;
            }

            // The Kmers for all ordinals in [firstOrdinalForAssembly, lastOrdinalForAssembly],
            // in that order.
            vector<Kmer> kmers;

            const Kmer& getKmer(uint32_t ordinal) const
            {
                SHASTA2_ASSERT(ordinal >= firstOrdinal);
                SHASTA2_ASSERT(ordinal <= lastOrdinal);
                return kmers[ordinal - firstOrdinal];
            }

            // The Kmers that appear more than once in this LocalRegion, sorted by Kmer.
            vector<Kmer> nonUniqueKmers;
            void computeNonUniqueKmers();
        };
        LocalRegion localRegion;
        void fillLocalRegion(const Anchors&, uint32_t length);


        // The ordinals of the Marker Kmers to be used for assembly.
        // Initially these are all the ordinals in [firstOrdinalForAssembly, lastOrdinalForAssembly],
        // except for the ones that correspond to Kmers that are non-unique in one or more reads.
        // Later, we remove ordinals that causes cycles in the graph.
        // These are sorted.
        vector<uint32_t> ordinalsForAssembly;

    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();

    uint32_t offset;
    void estimateOffset();

    // This fills the LocalRegion of all OrientedReadInfos.
    void fillLocalRegions();

    // The Kmers that will be used in this assembly.
    vector<Kmer> kmers;

    // Gather the Kmers that will be used in this assembly
    // and fill in the OrientedReadInfos;
    void gatherKmers();


    // Create the graph using the current ordinalsForAssembly.
    vertex_descriptor vLeft;
    vertex_descriptor vRight;
    void createGraph();

    // Compute position offset for an edge by averaging the
    // position offsets of the participating orientedReads.
    uint32_t edgeOffset(edge_descriptor) const;

    // Remove vertices that are not forward accessible from vLeft
    // and backward accessible from vRight.
    void removeInaccessibleVertices();

    // Compute a dominator tree starting at vLeft and the
    // dominator tree path between vLeft and vRight.
    void computeDominatorTree();
    vector<vertex_descriptor> dominatorTreePath;

    // Given two vertices which are adjacent in the dominator tree path, vA and vB,
    // remove low coverage edges in-between without destroying reachability
    // of vB from vA.
    void removeLowCoverageEdges(vertex_descriptor vA, vertex_descriptor vB);
    void removeLowCoverageEdges();

    // Compute the assembly path.
    void computeAssemblyPath();
    vector<vertex_descriptor> assemblyPath;

    void removeIsolatedVertices();

    // Html output.
    void writeInput(
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeAllOrientedReadIds() const;
    void writeMarkerInfos() const;
    void writeMarkerInfos(const string& side, const vector<MarkerInfo>&) const;
    void writeOrientedReads() const;
    void writeGraphviz(const string& fileName);
    void writeGraphviz(ostream&);
    void writeGraph();
    void writeKmers() const;
};
