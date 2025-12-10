#pragma once

// The read-anchor graph is an undirected bipartite graph
// in which each vertex represents an oriented read or an anchor.

// Shasta.
#include "Anchor.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include <map>
#include "vector.hpp"
#include "string.hpp"



namespace shasta2 {
    class LocalReadAnchorGraph;

    class LocalReadAnchorGraphVertex;

    using LocalReadAnchorGraphAssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::undirectedS,
        LocalReadAnchorGraphVertex>;

    class Journeys;
}



class shasta2::LocalReadAnchorGraphVertex
{
public:
    AnchorId anchorId = invalid<AnchorId>;
    OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);
    uint64_t distance = 0;

    // The position of this vertex in the computed layout.
    double x;
    double y;

    LocalReadAnchorGraphVertex(AnchorId anchorId, uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}

    LocalReadAnchorGraphVertex(OrientedReadId orientedReadId, uint64_t distance) :
        orientedReadId(orientedReadId),
        distance(distance)
    {}

    bool isAnchor() const
    {
        return anchorId != invalid<AnchorId>;
    }

    bool isOrientedRead() const
    {
        return not isAnchor();
    }


};



class shasta2::LocalReadAnchorGraph : public LocalReadAnchorGraphAssemblyGraphBaseClass{
public:

    // The constructor parses the request, creates the LocalReadAnchorGraph,
    // and displays it to html.
    LocalReadAnchorGraph(
        const Anchors&,
        const Journeys&,
        const vector<string>& request,
        ostream& html);

private:
    const Anchors& anchors;
    const Journeys& journeys;
    ostream& html;

    // Figure out if command "customLayout" is available.
    bool customLayoutIsAvailable = false;
    void findOutIfCustomLayoutIsAvailable();

    // The request options.
    string startVerticesString;
    uint64_t maxDistance = 3;
    uint64_t sizePixels = 600;
    string layoutMethod = "sfdp";
    void getRequestOptions(const vector<string>& request);

    // Write the page header.
    void writeHeader();

    // Write the form to enter the options.
    void writeForm();

    // The AnchorIds and OrientedReadIds for the starting vertices.
    vector<AnchorId> startAnchorIds;
    vector<OrientedReadId> startOrientedReadIds;
    void parseStartVertices();

    // Create the graph.
    void createVertices();
    void createEdges();
    std::map<AnchorId, vertex_descriptor> anchorVertexMap;
    std::map<OrientedReadId, vertex_descriptor> orientedReadVertexMap;

    // Compute the layout of the graph.
    // This stores a layout position in each vertex.
    void computeLayout();

    // The bounding box of the layout.
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
    void computeBoundingBox();

    // Display the graph to html.
    void write() const;
    void writeVertices() const;
    void writeEdges() const;
    void writeSvgControls() const;

};
