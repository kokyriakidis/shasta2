#pragma once

// Shasta.
// #include "invalid.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "vector.hpp"
#include "string.hpp"

namespace shasta2 {
    class LocalReadGraph;

    class LocalReadGraphVertex;

    using LocalReadGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::undirectedS,
        LocalReadGraphVertex>;

    class ReadGraph;

}



class shasta2::LocalReadGraphVertex
{
public:
    OrientedReadId orientedReadId;
    uint64_t distance = 0;

    // The position of this vertex in the computed layout.
    double x;
    double y;

    LocalReadGraphVertex(OrientedReadId orientedReadId, uint64_t distance) :
        orientedReadId(orientedReadId),
        distance(distance)
    {}

};



class shasta2::LocalReadGraph : public LocalReadGraphBaseClass{
public:

    // The constructor parses the request, creates the LocalReadGraph,
    // and displays it to html.
    LocalReadGraph(
        const ReadGraph&,
        const vector<string>& request,
        ostream& html);

private:
    const ReadGraph& readGraph;
    ostream& html;

    // Figure out if command "customLayout" is available.
    bool customLayoutIsAvailable = false;
    void findOutIfCustomLayoutIsAvailable();

    // The request options.
    string startOrientedReadIdsString;
    uint64_t maxDistance = 3;
    uint64_t sizePixels = 600;
    string layoutMethod = "sfdp";
    void getRequestOptions(const vector<string>& request);

    // Write the page header.
    void writeHeader();

    // Write the form to enter the options.
    void writeForm();

    // The OrientedReadIds for the starting vertices.
    vector<OrientedReadId> startOrientedReadIds;
    void parseStartOrientedReadIds();

    // Create the graph.
    void createVertices();
    void createEdges();
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

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

