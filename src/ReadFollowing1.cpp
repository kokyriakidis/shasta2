// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing1.hpp"
#include "Journeys.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>



ReadFollowing1::ReadFollowing1(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    findAppearances();
    countAppearances();
    createVertices();
    createEdges();

    // Enforce a minimum Jaccard when writing the graph.
    writeGraph(0.2);
}



void ReadFollowing1::findAppearances()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    initialAppearances.resize(orientedReadCount);
    finalAppearances.resize(orientedReadCount);

    // Loop over all Segments.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[segment];
        const uint64_t stepCount = edge.size();

        // Locate the initial representative region.
        const uint64_t initialBegin = 0;
        const uint64_t initialEnd = min(stepCount, representativeRegionLength);

        // Locate the final representative region.
        const uint64_t finalEnd = stepCount;
        const uint64_t finalBegin =
            ((stepCount >= representativeRegionLength) ? (stepCount - representativeRegionLength) : 0);

        // Appearances in the initial representative region of this edge.
        // For each OrientedReadId. store the last appearance in journey order.
        std::map<OrientedReadId, vector<AppearanceInfo> > initialAppearancesMap;
        for(uint64_t stepId=initialBegin; stepId!=initialEnd; stepId++) {

            // Compute the base offset to the end of the segment.
            uint64_t offset = 0;
            for(uint64_t i=stepId+1; i<initialEnd; i++) {
                const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
                if(assemblyGraphEdge.wasAssembled) {
                    offset += assemblyGraphEdge[i].sequence.size();
                } else {
                    offset += assemblyGraphEdge[i].offset;
                }
            }

            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdA;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                initialAppearancesMap[orientedReadId].push_back(AppearanceInfo(positionInJourney, stepId, offset));
            }
        }
        for(auto& p: initialAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<AppearanceInfo>& infos = p.second;
            sort(infos.begin(), infos.end());
            initialAppearances[orientedReadId.getValue()].push_back(
                Appearance(infos.back(), segment));
        }

        // Appearances in the final representative region of this edge.
        // For each OrientedReadId. store the first appearance in journey order.
        std::map<OrientedReadId, vector<AppearanceInfo> > finalAppearancesMap;
        for(uint64_t stepId=finalBegin; stepId!=finalEnd; stepId++) {

            // Compute the base offset from the beginning of the segment.
            uint64_t offset = 0;
            for(uint64_t i=finalBegin; i<stepId; i++) {
                const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
                if(assemblyGraphEdge.wasAssembled) {
                    offset += assemblyGraphEdge[i].sequence.size();
                } else {
                    offset += assemblyGraphEdge[i].offset;
                }
            }
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdB;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                finalAppearancesMap[orientedReadId].push_back(AppearanceInfo(positionInJourney, stepId, offset));
            }
        }
        for(auto& p: finalAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<AppearanceInfo>& infos = p.second;
            sort(infos.begin(), infos.end());
            finalAppearances[orientedReadId.getValue()].push_back(
                Appearance(infos.front(), segment));
        }
    }



    // For each OrientedReadId, sort the appearances by edge id.
    class AppearanceSorter {
    public:
        const AssemblyGraph& assemblyGraph;
        AppearanceSorter(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
        bool operator()(const Appearance& x, const Appearance& y) const
        {
            return assemblyGraph[x.segment].id < assemblyGraph[y.segment].id;
        }

    };
    for(vector<Appearance>& v: initialAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
    for(vector<Appearance>& v: finalAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
}



void ReadFollowing1::countAppearances()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();

    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Initial appearances.
        for(const Appearance appearance: initialAppearances[i]) {
            const Segment segment = appearance.segment;
            const auto it = initialAppearancesCount.find(segment);
            if(it == initialAppearancesCount.end()) {
                initialAppearancesCount.insert(make_pair(segment, 1));
            } else {
                ++(it->second);
            }
        }

        // Final appearances.
        for(const Appearance appearance: finalAppearances[i]) {
            const Segment segment = appearance.segment;
            const auto it = finalAppearancesCount.find(segment);
            if(it == finalAppearancesCount.end()) {
                finalAppearancesCount.insert(make_pair(segment, 1));
            } else {
                ++(it->second);
            }
        }

    }

}



uint64_t ReadFollowing1::getInitialAppearancesCount(Segment segment) const
{
    const auto it = initialAppearancesCount.find(segment);
    if(it == initialAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}



uint64_t ReadFollowing1::getFinalAppearancesCount(Segment segment) const
{
    const auto it = finalAppearancesCount.find(segment);
    if(it == finalAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}



// Create vertices of the ReadFollowing graph.
// Each vertex corresponds to a Segment, but not
// all Segments generate a vertex.
void ReadFollowing1::createVertices()
{
    Graph& graph = *this;

    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const ReadFollowing1Vertex vertex(assemblyGraph, segment);
        const vertex_descriptor v = add_vertex(ReadFollowing1Vertex(assemblyGraph, segment), graph);
        vertexMap.insert(make_pair(segment, v));
    }
}



ReadFollowing1Vertex::ReadFollowing1Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment) :
    segment(segment)
{
    const AssemblyGraphEdge& edge = assemblyGraph[segment];
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }
}



// Create edges of the ReadFollowing graph.
void ReadFollowing1::createEdges()
{
    Graph& graph = *this;
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();

    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Look over pairs (final appearance, initial appearance).
        for(const Appearance appearance0: finalAppearances[i]) {
            const Segment segment0 = appearance0.segment;
            const auto it0 = vertexMap.find(segment0);
            if(it0 == vertexMap.end()) {
                continue;
            }
            const vertex_descriptor v0 = it0->second;

            for(const Appearance appearance1: initialAppearances[i]) {
                const Segment segment1 = appearance1.segment;
                if(segment1 == segment0) {
                    continue;
                }
                const auto it1 = vertexMap.find(segment1);
                if(it1 == vertexMap.end()) {
                    continue;
                }
                const vertex_descriptor v1 = it1->second;

                if((appearance0.positionInJourney >= appearance1.positionInJourney)) {
                    continue;
                }

                // Increment the coverage of edge segment0->segment1,
                // creating the edge if necessary.
                edge_descriptor e;
                bool edgeExists;
                tie(e, edgeExists) = boost::edge(v0, v1, graph);
                if(not edgeExists) {
                    tie(e, edgeExists) = add_edge(v0, v1, graph);
                }
                ++graph[e].coverage;
            }
        }
    }


    // Fill in the Jaccard similarities of all the edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        graph[e].jaccard = jaccard(e);
    }

}


// Jaccard similarity for an EdgePairsGraph edge.
// Computed using the finalAppearancesCount of the source vertex
// and the initialAppearancesCount of the target vertex.
double ReadFollowing1::jaccard(edge_descriptor e) const
{
    const Graph& graph = *this;

    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);

    const Segment segment0 = graph[v0].segment;
    const Segment segment1 = graph[v1].segment;

    const uint64_t n0 = getFinalAppearancesCount(segment0);
    const uint64_t n1 = getInitialAppearancesCount(segment1);

    const uint64_t n01 = graph[e].coverage;

    const uint64_t intersectionSize = n01;
    const uint64_t unionSize = n0 + n1 - intersectionSize;

    return double(intersectionSize) / double(unionSize);

}



void ReadFollowing1::writeGraph(double minJaccard) const
{
    const Graph& graph = *this;

    ofstream dot("ReadFollowing1.dot");
    dot << "digraph ReadFollowing1 {\n";
    dot << std::fixed << std::setprecision(2);

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const ReadFollowing1Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id <<
            " [label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            getInitialAppearancesCount(segment) << "/" <<
            getFinalAppearancesCount(segment) <<
            "\"";
        dot <<
            "]"
            ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;
        const uint64_t coverage = graph[e].coverage;
        const double j = graph[e].jaccard;
        SHASTA_ASSERT(j <= 1.);

        if(j < minJaccard) {
            continue;
        }

        const double hue = j / 3.; // So 0=red, 1=green.

        dot <<
            assemblyGraph[segment0].id << "->" <<
            assemblyGraph[segment1].id <<
            " [penwidth=\"" << 0.4 * double(coverage) << "\""
            " tooltip=\"" << coverage << " " << j << "\""
            " color=\"" << hue << ",1,.9\""
            "]"
            ";\n";
    }

    dot << "}\n";
}



void ReadFollowing1::findPath(Segment segment, uint64_t direction, vector<vertex_descriptor>& path) const
{
    if(direction == 0) {
        findForwardPath(segment, path);
    } else if(direction == 1) {
        findBackwardPath(segment, path);
    } else {
        SHASTA_ASSERT(0);
    }
}



void ReadFollowing1::findForwardPath(Segment segment, vector<vertex_descriptor>& path) const
{
    const Graph& graph = *this;

    // Find the start vertex.
    const auto it = vertexMap.find(segment);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Each iteration adds one vertex to the path.
    path.clear();
    path.push_back(v);
    std::set<vertex_descriptor> pathVertices;
    pathVertices.insert(v);
    while(true) {

         // Find the best next vertex.
         vertex_descriptor vNext = null_vertex();
         double bestJaccardSum = 0.;
         BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
             vertex_descriptor v1 = target(e, graph);
             if(pathVertices.contains(v1)) {
                 continue;
             }

             // Compute the sum of Jaccard similarities between all
             // previous vertices in the path and v1.
             double jaccardSum = 0.;
             for(const vertex_descriptor vPrevious: path) {
                 edge_descriptor ePrevious;
                 bool edgeWasFound;
                 tie(ePrevious, edgeWasFound) = edge(vPrevious, v1, graph);
                 if(edgeWasFound) {
                    jaccardSum += graph[ePrevious].jaccard;
                 }
             }

             if(jaccardSum > bestJaccardSum) {
                 vNext = v1;
                 bestJaccardSum = jaccardSum;
             }

         }

         if(vNext == null_vertex()) {
             break;
         }

         v = vNext;

         path.push_back(v);
         pathVertices.insert(v);
     }

}



void ReadFollowing1::findBackwardPath(Segment, vector<vertex_descriptor>& /* path */) const
{
    SHASTA_ASSERT(0);
}



void ReadFollowing1::writePath(Segment segment, uint64_t direction) const
{
    const Graph& graph = *this;

    vector<vertex_descriptor> path;
    findPath(segment, direction, path);

    if(direction == 1) {
        std::ranges::reverse(path);
    }

    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << " ";
    }
    cout << endl;
}
