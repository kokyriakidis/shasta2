#pragma once

// Shasta.
#include "Anchor.hpp"
#include "AnchorGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library;
#include <map>
#include "vector.hpp"

namespace shasta2 {
    class LocalAnchorGraph;
    class LocalAnchorGraphVertex;
    class LocalAnchorGraphEdge;

    using LocalAnchorGraphAssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalAnchorGraphVertex,
        LocalAnchorGraphEdge>;

    class LocalAnchorGraphDisplayOptions;

    class AnchorGraph;
    class AssemblyGraphPostprocessor;
}



class shasta2::LocalAnchorGraphDisplayOptions {
public:
    uint64_t sizePixels;
    string layoutMethod;

    // Vertices.
    double vertexSize;
    bool vertexSizeByCoverage;
    bool vertexLabels;
    string vertexColoring;
    string similarityMeasure;
    string referenceAnchorIdString;
    string assemblyStage;

    // Edges.
    string edgeColoring;
    double edgeThickness;
    double minimumEdgeLength;
    double additionalEdgeLengthPerKb;
    double arrowSize;
    bool edgeLabels;

    // Construct from an html request.
    LocalAnchorGraphDisplayOptions(const vector<string>& request);

    // Write the form.
    void writeForm(ostream& html) const;
};



class shasta2::LocalAnchorGraphVertex {
public:
    AnchorId anchorId;
    uint64_t distance;
    LocalAnchorGraphVertex(
        AnchorId anchorId,
        uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}
};


class shasta2::LocalAnchorGraphEdge {
public:
    // The edge of the global AnchorGraph that corresponds to this LocalAnchorGraphEdge.
    AnchorGraph::edge_descriptor eG;

};



class shasta2::LocalAnchorGraph : public LocalAnchorGraphAssemblyGraphBaseClass {
public:

    LocalAnchorGraph(
        const Anchors&,
        const AnchorGraph&,
        const vector<AnchorId>&,
        uint64_t maxDistance,
		uint64_t minCoverage,
        bool edgesMarkedForAssembly);

    const Anchors& anchors;
    const AnchorGraph* anchorGraphPointer = 0;
    uint64_t maxDistance;
    std::map<AnchorId, vertex_descriptor> vertexMap;

    void writeHtml(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*);
    void writeHtml1(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*) const;

    void writeGraphviz(
        const string& fileName,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*) const;
    void writeGraphviz(
        ostream&,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*) const;

private:

    // Html/svg output without using svg output created by Graphviz.
    void writeHtml2(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*);

    // The position of each vertex in the computed layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    void computeLayout(const LocalAnchorGraphDisplayOptions&);

    // The bounding box of the computed layout.
    class Box {
    public:
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        double xSize() {return xMax - xMin;}
        double ySize() {return yMax - yMin;}
        void makeSquare();
        void extend(double factor);
    };
    Box boundingBox;
    void computeLayoutBoundingBox();

    void writeVertices(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&,
        const AssemblyGraphPostprocessor*) const;

    void writeEdges(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&) const;

    void writeSvgControls(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&) const;
};
