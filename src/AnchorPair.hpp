#pragma once

// Shasta.
#include "invalid.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Standard library.
#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {
    class AnchorPair;
    class Anchors;
    class Base;
    class Journeys;
}



// An AnchorPair is a set of OrientedReadIds that visit anchorIdA
// and then, later in their Journey, anchorIdB.
// That is, each of the OrientedReadIds appear both in anchorIdA and anchorIdB,
// and the position at anchorIdB is greater than the position at anchorIdA.
// The AnchorPair can use a subset of all possible OrientedReadIds
// that satisfy the above.
class shasta2::AnchorPair {
public:

    // The constructor creates an AnchorPair between anchorIdA and anchorIdB.
    // - If adjacentInJourney is false, it includes all OrientedReadIds
    //  that visit anchorIdA and then, later in their Journey, anchorIdB.
    // - If adjacentInJourney is true, only OrientedReadIds that visit anchorIdB
    //   immediately after visiting anchorIdA are included.
    //   That is, the journey offset between anchorIdA and anchorIdB
    //   for OrientedReadIds that are included must be exactly 1.
    AnchorPair(
        const Anchors&,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        bool adjacentInJourney);

    AnchorPair() {}

    AnchorPair(const AnchorPair& that) :
        anchorIdA(that.anchorIdA),
        anchorIdB(that.anchorIdB),
        orientedReadIds(that.orientedReadIds)
    {}

    // Copy from another AnchorPair, but excluding some OrientedReadIds.
    AnchorPair(
        const AnchorPair&,
        const vector<OrientedReadId>& excludedOrientedReadIds);

    // "Join" constructor from two AnchorPairs.
    // This constructs a new AnchorPair as follows:
    // - anchorIdA is the same as anchorPair0.anchorIdA.
    // - anchorIdB is the same as anchorPair1.anchorIdB.
    // - orientedReadIds are the intersection of
    //   anchorPair0.orientedReadIds and anchorPair1.orientedReadIds,
    //   with the additional requirement that the new journey offset
    //   is positive. That is, each OrientedReadId of the new AnchorPair
    //   visits anchorIdB after visiting the new anchorIdA.
    // This constructor is used in detangling.
    AnchorPair(
        const Anchors&,
        const AnchorPair& anchorPair0,
        const AnchorPair& anchorPair1);

    // This finds AnchorPairs as follows:
    // - anchorIdA is as specified.
    // - Coverage is at least minCoverage.
    // - All oriented reads have a journey offset equal to 1.
    static void createChildren(
        const Anchors&,
        const Journeys&,
        AnchorId anchorIdA,
        uint64_t minCoverage,
        vector<AnchorPair>&
        );


    AnchorId anchorIdA = invalid<AnchorId>;
    AnchorId anchorIdB = invalid<AnchorId>;

    vector<OrientedReadId> orientedReadIds;
    uint64_t size() const
    {
        return orientedReadIds.size();
    }

    uint32_t getAverageOffset(const Anchors&) const;
    void getOffsets(
        const Anchors&,
        uint32_t& averageBaseOffset,
        uint32_t& minBaseOffset,
        uint32_t& maxBaseOffset) const;

    // Get positions in journey, ordinals, base positions
    // for each of the two reads and for each of the two anchors.
    // The positions returned are the midpoint of the markers
    // corresponding to anchorIdA and anchorIdB.
    class Positions {
    public:
        uint32_t positionInJourney;
        uint32_t ordinal;
        uint32_t basePosition;
        Positions(
            uint32_t positionInJourney,
            uint32_t ordinal,
            uint32_t basePosition) :
            positionInJourney(positionInJourney),
            ordinal(ordinal),
            basePosition(basePosition)
        {}

        template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
        {
            ar & positionInJourney;
            ar & ordinal;
            ar & basePosition;
        }
    };
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions) const;

    // Same as the above, but also returns compute the sequences.
    void get(
        const Anchors&,
        vector< pair<Positions, Positions> >& positions,
        vector< vector<Base> >&) const;

    // Just return the ordinals.
    void getOrdinals(const Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Just return the positions in journeys.
    void getPositionsInJourneys(const Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Split the AnchorPair into one or more AnchorPairs with consistent offsets.
    // In the resulting AnchorPairs, if the position offsets are sorted in
    // increasing order, any two adjacent offsets D0 and D1
    // will satisfy D1 - D0 <= aDrift + bDrift * (D0 + D1) / 2.
    void splitByOffsets(
        const Anchors&,
        double aDrift,
        double bDrift,
        vector<AnchorPair>&
        ) const;

    // Split the AnchorPair using clustering of the oriented read journey portions
    // within this AnchorPair.
    // The new AnchorPairs are sorted by decreasing size.
    void splitByClustering(
        const Anchors&,
        const Journeys&,
        double clusteringMinJaccard,
        vector<AnchorPair>&
        ) const;

    // This returns true if a call to split with the same arguments would not split this Anchor.
    // The last two arguments are work vectors added as arguments for performance.
    bool isConsistent(
        const Anchors&,
        double aDrift,
        double bDrift,
        vector< pair<Positions, Positions> >&,
        vector<uint64_t>&) const;

    // Count OrientedReadIds in common with another AnchorPair.
    uint64_t countCommon(const AnchorPair&) const;

    // Remove from the AnchorPair OrientedReadIds that have negative offsets.
    void removeNegativeOffsets(const Anchors&);

    bool contains(OrientedReadId) const;



    // Various functions to get information about the AnchorIds present
    // in the portions between anchorIdA and anchorIdB of the journeys
    // of the OrientedReadIds in this AnchorPair.
    // They all return vectors sorted by AnchorId.
    // All = ancorIdA and anchorIdB are included.
    // Internal = ancorIdA and anchorIdB are excluded.
    // LocalCoverage = counting only the portions between anchorIdA and anchorIdB of the journeys
    // of the OrientedReadIds in this AnchorPair.
    void getAllAnchorIds(
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
        vector<AnchorId>&) const;
    void getInternalAnchorIds(
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
        vector<AnchorId>&) const;
    void getAllAnchorIdsAndLocalCoverage(
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
        vector<AnchorId>&,
        vector<uint64_t>& localCoverage) const;



    // Compute the clustering matrix.
    // This is a matrix with one row for each OrientedReadId
    // and a column for each internal AnchorId.
    // clusteringMatrix(i, j) is 1 if orientedReadId[i] visits internalAnchorIds[j]
    // in its journey portion between anchorIdA and anchorIdB.
    // Here, internalAnchorIds is as computed by getInternalAnchorIds
    // and is sorted by AnchorId.
    // The cluster matrix is stores as a column-major matrix
    // (Fortran compatible storage layout) so it can be later used
    // for a Singular Value Decomposition (SVD) for clustering.
    using Matrix = boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>;
    void computeClusteringMatrix(
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,  // As computed by getPositionsInJourneys.
        const vector<AnchorId>& internalAnchorIds,                      // As computed by getInternalAnchorIds.
        Matrix&
        ) const;

    // Singular value decomposition of the clustering matrix.
    void clusteringMatrixSvd(
        Matrix& clusteringMatrix,
        vector<double>& singularValues,
        Matrix& leftSingularVectors,
        Matrix& rightSingularVectors) const;

    // Use the scaled left singular values to compute a distance matrix
    // between oriented reads.
    void computeDistanceMatrix(
        uint64_t singularValueCount,    // Only use the first singular values
        const vector<double>& singularValues,
        const Matrix& leftSingularVectors,
        Matrix& distanceMatrix
        ) const;

    // Given the distance matrix, compute a similarity graph
    // between OrientedReadIds in which each vertex represents an OrientedReadId and
    // an edge is generated for OrientedReadId pairs
    // with distance below the given threshold.
    // Each vertex stores the id of the cluster it is assigned to.
    class OrientedReadIdSimilarityGraph :
        public boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, uint64_t> {
    public:
        OrientedReadIdSimilarityGraph(const Matrix& distanceMatrix, double maxDistance);
        vector< vector<vertex_descriptor> > clusters;
        void writeHtml(ostream&, const vector<OrientedReadId>&) const;
    private:
        void writeGraphviz(const string& fileName, const vector<OrientedReadId>&) const;
        void writeGraphviz(ostream&, const vector<OrientedReadId>&) const;
    };
    void writeClustersHtml(ostream&, const OrientedReadIdSimilarityGraph&) const;


    // A simple local anchor graph constructed using only the portions
    // between anchorIdA and anchorIdB of the journeys
    // of the OrientedReadIds in this AnchorPair.
    class SimpleLocalAnchorGraphVertex {
    public:
        AnchorId anchorId = invalid<AnchorId>;
        uint64_t localCoverage = invalid<uint64_t>;
        uint64_t color = 0;     // Only used by approximateTopologicalSort.
        uint64_t rank = 0;      // Only used by approximateTopologicalSort.
        SimpleLocalAnchorGraphVertex(AnchorId anchorId, uint64_t localCoverage) :
            anchorId(anchorId), localCoverage(localCoverage) {}
        SimpleLocalAnchorGraphVertex() {}
    };
    class SimpleLocalAnchorGraphEdge {
    public:
        uint64_t localCoverage = 0;
        bool isDagEdge = false; // Only used by approximateTopologicalSort.
    };
    class SimpleLocalAnchorGraph :
        public boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        SimpleLocalAnchorGraphVertex, SimpleLocalAnchorGraphEdge> {
    public:
        SimpleLocalAnchorGraph(const Anchors&, const Journeys&, const AnchorPair&);

        // The vertex descriptors, sorted in approximate topological order.
        void approximateTopologicalSort();
        vector<vertex_descriptor> approximateTopologicalOrder;

        void getInternalAnchorIdsInTopologicalOrder(vector<AnchorId>&) const;

        void writeGraphviz(const string& fileName) const;
        void writeGraphviz(ostream&) const;
    };



    // Return the url for the exploreAnchorPair1 page for this AnchorPair.
    string url() const;



    // Html output.
    void writeAllHtml(ostream&, const Anchors&, const Journeys&) const;
    void writeSummaryHtml(ostream&, const Anchors&) const;
    void writeOrientedReadIdsHtml(ostream&, const Anchors&) const;
    void writeJourneysHtml(
        ostream&,
        const Journeys&,
        const vector< pair<uint32_t, uint32_t> >& positionsInJourneys   // As computed by getPositionsInJourneys.
        ) const;
    void writeClusteringMatrix(
        ostream& html,
        // The internalAnchorIds as computed by getInternalAnchorIds.
        const vector<AnchorId>& internalAnchorIds,
        // The same AnchorIds, in the order in which the corresponding columns should be written out
        const vector<AnchorId>& internalAnchorIdsInOutputOrder,
        const Matrix&) const;
    void writeClusteringMatrixSvd(
        ostream& html,
        // The internalAnchorIds as computed by getInternalAnchorIds.
        const vector<AnchorId>& internalAnchorIds,
        // The same AnchorIds, in the order in which the corresponding columns should be written out
        const vector<AnchorId>& internalAnchorIdsInOutputOrder,
        uint64_t singularValueCount,    // Number of singular values/vectors to be writtten
        const vector<double>& singularValues,
        const Matrix& leftSingularVectors,
        const Matrix& rightSingularVectors) const;
    void writeDistanceMatrixHtml(
        ostream& html,
        const Matrix& distanceMatrix) const;
    void writeSimpleLocalAnchorGraphHtml(ostream&, const SimpleLocalAnchorGraph&) const;



    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIds;
    }
};
