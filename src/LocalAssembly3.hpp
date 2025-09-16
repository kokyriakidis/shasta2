#pragma once

// Shasta.
#include "Kmer.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    class LocalAssembly3;
    class LocalAssembly3Vertex;
    class LocalAssembly3Edge;

    using LocalAssembly3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalAssembly3Vertex,
        LocalAssembly3Edge>;

    class Anchors;
    class AnchorPair;
    class Markers;
}



class shasta::LocalAssembly3Vertex {
public:

    // The index of the Kmer corresponding to this vertex
    // in the kmers vector.
    uint64_t kmerIndex;

    // The markers that contain this Kmer.
    class Data {
    public:
        uint64_t orientedReadIndex; // Index in orientedReadInfos vector.
        uint32_t ordinal;
    };
    vector<Data> data;
    uint64_t coverage() const
    {
        return data.size();
    }

    LocalAssembly3Vertex(uint64_t kmerIndex) : kmerIndex(kmerIndex) {}

    // These are used when computing the dominatyor tree.
    LocalAssembly3BaseClass::vertex_descriptor dominator =
        LocalAssembly3BaseClass::null_vertex();
    bool isOnDominatorTreePath = false;
};



class shasta::LocalAssembly3Edge {
public:

    // The oriented reads that transition from the source of this vertex
    // to its target.
    class Data {
    public:
        uint64_t orientedReadIndex; // Index in orientedReadInfos vector.
        uint32_t ordinal0;
        uint32_t ordinal1;
    };
    vector<Data> data;
    uint64_t coverage() const
    {
        return data.size();
    }
};



class shasta::LocalAssembly3 : public LocalAssembly3BaseClass {
public:
    using Base = shasta::Base;

    // This assembles between anchorIdA and anchorIdB
    // of the given AnchorPair. It uses all the OrientedReadIds
    // stored in the AnchorPair, and which all appear in both
    // anchorIdA and anchorIdB. In addition, it uses OrientedReadIds
    // stored in additionalOrientedReadIds that:
    // - Are not also in the AnchorPair.
    // - Appear in at least one of anchorIdA and anchorIdB.
    // The additionalOrientedReadIds must be sorted.
    LocalAssembly3(
        const Anchors&,
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

public:

    // The two anchors that bound the region assembled by this local assembly.
    AnchorId leftAnchorId;
    AnchorId rightAnchorId;



    // Information for each oriented read used in this local assembly.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // Whether this OrientedReadId appear in the left and right anchors.
        bool isOnLeftAnchor;
        bool isOnRightAnchor;
        bool isOnBothAnchors() const
        {
            return isOnLeftAnchor and isOnRightAnchor;
        }

        // The ordinals of this OrientedReadId in the left/right anchor, if any.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            SHASTA_ASSERT(isOnLeftAnchor);
            SHASTA_ASSERT(isOnRightAnchor);
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right anchor, if any.
        // These are positions in the oriented read sequence
        // of the leftmost base of the marker.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            SHASTA_ASSERT(isOnLeftAnchor);
            SHASTA_ASSERT(isOnRightAnchor);
            return rightPosition - leftPosition;
        }

        // All data members up to here are filled by gatherOrientedReads.

        // The first and last ordinal of the portion of this OrientedReadId
        // sequence that will be used in this local assembly.
        // These are filled by fillFirstLastOrdinalForAssembly.
        uint32_t firstOrdinalForAssembly;
        uint32_t lastOrdinalForAssembly;
        void fillFirstLastOrdinalForAssembly(const Markers&, uint32_t length);
        uint32_t firstPositionForAssembly(const Markers&) const;
        uint32_t lastPositionForAssembly(const Markers&) const;

        // The Kmers for each of the ordinal in [firstOrdinalForAssembly, lastOrdinalForAssembly].
        // These are filled by fillOrientedReadKmers.
        class OrientedReadKmerInfo {
        public:
            Kmer kmer;
            uint64_t kmerIndex; // The index of this Kmer in the kmers vector.
        };
        vector<OrientedReadKmerInfo> orientedReadKmerInfos;
        void fillOrientedReadKmers(const Markers&);
    };

    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads(
        const Anchors&,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // The estimated base offset between the left and right anchor.
    // It is estimated using the oriented reads that appear in both anchors.
    uint32_t offset;
    void estimateOffset();

    // Use the estimated offset to fill in the firstOrdinal and lastOrdinal
    // of each oriented read. These define the portion of this OrientedReadId
    // sequence that will be used in this local assembly.
    void fillFirstLastOrdinalForAssembly(const Markers&, double drift);

    void fillOrientedReadKmers(const Markers&);

    // A sorted vector containing all the Kmers found in all oriented reads
    // of this local assembly. This provides a perfect hash function for
    // these Kmers.
    vector<Kmer> kmers;
    void gatherKmers(const Anchors&);

    // Return the index of the given Kmer in the kmers vector.
    uint64_t getKmerIndex(const Kmer&) const;

    // The kmer indexes of the left and right anchors.
    uint64_t leftAnchorKmerIndex;
    uint64_t rightAnchorKmerIndex;

    // Map k-mer indexes to vertices.
    vector<vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();
    vertex_descriptor leftAnchorVertex;
    vertex_descriptor rightAnchorVertex;

    // Compute the average position offset for an edge.
    uint32_t edgePositionOffset(
        edge_descriptor,
        const Markers&) const;

    // Compute a dominator tree stasting at leftAnchorVertex
    // and the dominator tree path from leftAnchorVertex
    // to rightAnchorVertex.
    void computeDominatorTree();
    vector<vertex_descriptor> dominatorTreePath;

    // Assemble sequence. Use the dominatorTreePath as the assembly path.
    void assemble(
        const Anchors&,
        ostream& html,
        bool debug);
    void assemble(
        const Anchors&,
        vertex_descriptor,
        vertex_descriptor,
        ostream& html,
        bool debug);

    // Assembled sequence and its coverage.
    vector<Base> sequence;
    vector<uint64_t> coverage;

    // Html output.
    void writeInput(
        ostream& html,
        bool debug,
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeOrientedReads(
        const Anchors&,
        ostream& html) const;
    void writeKmers(ostream& html, uint64_t k) const;
    void writeOrientedReadKmers(ostream& html) const;
    void writeGraphviz(const string& fileName, const Markers&) const;
    void writeGraphviz(ostream&, const Markers&) const;
    void writeHtml(ostream&, const Markers&) const;
    void writeConsensus(
        ostream& html,
        const vector< pair<Base, uint64_t> >& consensus,
        uint64_t maxCoverage) const;
    void writeAlignment(
        ostream& html,
        const vector< vector<Base> >& inputSequences,
        const vector< pair<Base, uint64_t> >& consensus,
        const vector< vector<AlignedBase> >& alignment,
        const vector<AlignedBase>& alignedConsensus,
        uint64_t maxCoverage) const;
};
