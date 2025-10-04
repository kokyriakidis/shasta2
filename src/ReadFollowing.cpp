// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing.hpp"
#include "deduplicate.hpp"
#include "findReachableVertices.hpp"
#include "longestPath.hpp"
#include "Journeys.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>



ReadFollowing::ReadFollowing(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    findAppearances();
    countAppearances();
    findEdgePairs();
    createGraph();
    writeGraph();
    findAssemblyPaths();
}



void ReadFollowing::findAppearances()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    initialAppearances.resize(orientedReadCount);
    finalAppearances.resize(orientedReadCount);



    // Loop over all edges of the AssemblyGraph.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
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
        std::map<OrientedReadId, vector<uint32_t> > initialAppearancesMap;
        for(uint64_t stepId=initialBegin; stepId!=initialEnd; stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdA;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                initialAppearancesMap[orientedReadId].push_back(positionInJourney);
            }
        }
        for(auto& p: initialAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<uint32_t>& positionsInJourney = p.second;
            std::ranges::sort(positionsInJourney);
            initialAppearances[orientedReadId.getValue()].push_back(
                Appearance(e, positionsInJourney.back()));
        }

        // Appearances in the final representative region of this edge.
        // For each OrientedReadId. store the first appearance in journey order.
        std::map<OrientedReadId, vector<uint32_t> > finalAppearancesMap;
        for(uint64_t stepId=finalBegin; stepId!=finalEnd; stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdB;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                finalAppearancesMap[orientedReadId].push_back(positionInJourney);
            }
        }
        for(auto& p: finalAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<uint32_t>& positionsInJourney = p.second;
            std::ranges::sort(positionsInJourney);
            finalAppearances[orientedReadId.getValue()].push_back(
                Appearance(e, positionsInJourney.front()));
        }
    }



    // For each OrientedReadId, sort the appearances by edge id.
    class AppearanceSorter {
    public:
        const AssemblyGraph& assemblyGraph;
        AppearanceSorter(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
        bool operator()(const Appearance& x, const Appearance& y) const
        {
            return assemblyGraph[x.e].id < assemblyGraph[y.e].id;
        }

    };
    for(vector<Appearance>& v: initialAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
    for(vector<Appearance>& v: finalAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
}



void ReadFollowing::countAppearances()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();


    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Initial appearances.
        for(const Appearance appearance: initialAppearances[i]) {
            const AEdge e = appearance.e;
            const auto it = initialAppearancesCount.find(e);
            if(it == initialAppearancesCount.end()) {
                initialAppearancesCount.insert(make_pair(e, 1));
            } else {
                ++(it->second);
            }
        }

        // Final appearances.
        for(const Appearance appearance: finalAppearances[i]) {
            const AEdge e = appearance.e;
            const auto it = finalAppearancesCount.find(e);
            if(it == finalAppearancesCount.end()) {
                finalAppearancesCount.insert(make_pair(e, 1));
            } else {
                ++(it->second);
            }
        }

    }

}



void ReadFollowing::findEdgePairs()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();

    vector<AssemblyGraphEdgePair> edgePairsVector;

    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Look over pairs (final appearance, initial appearance).
        for(const Appearance appearance0: finalAppearances[i]) {
            const AEdge e0 = appearance0.e;
            for(const Appearance appearance1: initialAppearances[i]) {
                const AEdge e1 = appearance1.e;
                if((e1 != e0) and (appearance0.positionInJourney < appearance1.positionInJourney)) {
                    edgePairsVector.push_back(make_pair(e0, e1));
                }
            }
        }
    }
    vector<uint64_t> coverage;
    deduplicateAndCount(edgePairsVector, coverage);
    SHASTA_ASSERT(coverage.size() == edgePairsVector.size());

    // Now store them in the edgePairs map.
    for(uint64_t i=0; i<edgePairsVector.size(); i++) {
        edgePairs.insert(make_pair(edgePairsVector[i], coverage[i]));
    }

    // Write out the edgePairs.
    {
        ofstream csv("ReadFollowing-EdgePairs.csv");
        csv << "Id0,Id1,Coverage,\n";
        for(const auto& p: edgePairs) {
            const auto& edgePair = p.first;
            const uint64_t coverage = p.second;
            const AEdge e0 = edgePair.first;
            const AEdge e1 = edgePair.second;
            csv << assemblyGraph[e0].id << ",";
            csv << assemblyGraph[e1].id << ",";
            csv << coverage << ",";
            csv << "\n";
        }
    }
}



void ReadFollowing::writeGraph() const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    ofstream dot("ReadFollowing.dot");
    dot << "digraph EdgePairsGraph {\n";
    dot << std::fixed << std::setprecision(2);

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const ReadFollowingVertex& vertex = graph[v];
        const AEdge ae = vertex.ae;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[ae];
        dot << assemblyGraphEdge.id <<
            " [label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            getInitialAppearancesCount(ae) << "/" <<
            getFinalAppearancesCount(ae) <<
            "\"";
        if(vertex.isLong) {
            dot << " style=filled fillcolor=CornflowerBlue";
        }
        dot <<
            "]"
            ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AEdge ae0 = graph[v0].ae;
        const AEdge ae1 = graph[v1].ae;
        const uint64_t coverage = graph[e].coverage;

        const double j = jaccard(e);
        SHASTA_ASSERT(j <= 1.);
        const double hue = j / 3.; // So 0=red, 1=green.

        dot <<
            assemblyGraph[ae0].id << "->" <<
            assemblyGraph[ae1].id <<
            " [penwidth=\"" << 0.4 * double(coverage) << "\""
            " tooltip=\"" << coverage << " " << j << "\""
            " color=\"" << hue << ",1,.9\""
            "]"
            ";\n";
    }

    dot << "}\n";
}



void ReadFollowing::createGraph()
{
    using Graph = ReadFollowing;
    Graph& graph = *this;
    const bool debug = false;

    // Create a vertex for each AssemblyGraph edge.
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        vertexMap.insert(make_pair(ae,
            add_vertex(ReadFollowingVertex(assemblyGraph, ae, longLengthThreshold), graph)));
    }

    // Group AssemblyGraphEdgePairs (e0, e1)
    // by (v0, v1), where v0 = target(e0), v1 = source(e1).
    std::set<AVertexPair> s;
    for(const auto& p: edgePairs) {
        const AEdgePair& aEdgePair = p.first;
        const AEdge e0 = aEdgePair.first;
        const AEdge e1 = aEdgePair.second;
        const AVertex v0 = target(e0, assemblyGraph);
        const AVertex v1 = source(e1, assemblyGraph);
        s.insert(make_pair(v0, v1));
    }

    // Loop over pairs of AssemblyGraph vertices (v0, v1) such that there
    // is at least one EdgePair (e0, e1) with v0 = target(e0, assemblyGraph)
    // and v1 = source(e1, assemblyGraph).
    for(const auto& p: s) {
        const AVertex v0 = p.first;
        const AVertex v1 = p.second;

        // Gather the in-edges of v0.
        vector<AEdge> in;
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            in.push_back(e);
        }
        sort(in.begin(), in.end(), assemblyGraph.orderById);
        const uint64_t nIn = in.size();

        // Gather the out-edges of v1.
        vector<AEdge> out;
        BGL_FORALL_OUTEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            out.push_back(e);
        }
        sort(out.begin(), out.end(), assemblyGraph.orderById);
        const uint64_t nOut = out.size();

        if(debug) {
            cout << "In:";
            for(const AEdge e: in) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;

            cout << "Out:";
            for(const AEdge e: out) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
        }

        // Create a tangle matrix using the edgePairs.
        vector< vector<uint64_t> > tangleMatrix(nIn, vector<uint64_t>(nOut));

        if(debug) {
            cout << "Tangle matrix:" << endl;
        }
        for(uint64_t i0=0; i0<nIn; i0++) {
            const AEdge e0 = in[i0];
            for(uint64_t i1=0; i1<nOut; i1++) {
                const AEdge e1 = out[i1];
                uint64_t coverage = 0;
                auto it = edgePairs.find(make_pair(e0, e1));
                if(it != edgePairs.end()) {
                    coverage = it->second;
                }
                if(debug) {
                    cout << coverage << ",";
                }
                tangleMatrix[i0][i1] = coverage;
            }
            if(debug) {
                cout << endl;
            }
        }

        // Compute the sum of tangle matrix values for each of the in-edges.
        vector<uint64_t> inSum(nIn, 0);
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                inSum[i0] += tangleMatrix[i0][i1];
            }
        }

        if(debug) {
            cout << "inSum: ";
            for(const uint64_t s: inSum) {
                cout << s << ",";
            }
            cout << endl;
        }

        // Compute the sum of tangle matrix values for each of the out-edges.
        vector<uint64_t> outSum(nOut, 0);
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                outSum[i1] += tangleMatrix[i0][i1];
            }
        }

        if(debug) {
            cout << "outSum: ";
            for(const uint64_t s: outSum) {
                cout << s << ",";
            }
            cout << endl;
        }



        // Each "strong" element of the tangle matrix generates a strong edge pair
        // and an edge of the EdgePairsGraph.
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                const AEdge e0 = in[i0];
                const AEdge e1 = out[i1];

                if(debug) {
                    cout << "Checking " << assemblyGraph[e0].id << " " << assemblyGraph[e1].id << endl;
                }

                const uint64_t coverage = tangleMatrix[i0][i1];
                if(coverage < minCoverage) {
                    if(debug) {
                        cout << "Discarded due to coverage." << endl;
                    }
                    continue;
                }
                if(double(coverage) / double(inSum[i0]) < minCoverageFraction) {
                    if(debug) {
                        cout << "Discarded due to coverage fraction on " << assemblyGraph[e0].id << endl;
                    }
                    continue;
                }
                if(double(coverage) / double(outSum[i1]) < minCoverageFraction) {
                    if(debug) {
                        cout << "Discarded due to coverage fraction on " << assemblyGraph[e1].id << endl;
                    }
                    continue;
                }

                if(graph[vertexMap[e0]].length < shortLengthThreshold) {
                    if(debug) {
                        cout << "Discarded due to length on " << assemblyGraph[e0].id << endl;
                    }
                    continue;
                }
                if(graph[vertexMap[e1]].length < shortLengthThreshold) {
                    continue;
                    if(debug) {
                        cout << "Discarded due to length on " << assemblyGraph[e1].id << endl;
                    }
                }

                if(initialAppearancesCount[e0] > maxAppearanceCount) {
                    if(debug) {
                        cout << "Discarded due to initial coverage on " << assemblyGraph[e0].id << endl;
                    }
                    continue;
                }
                if(finalAppearancesCount[e0] > maxAppearanceCount) {
                    if(debug) {
                        cout << "Discarded due to final coverage on " << assemblyGraph[e0].id << endl;
                    }
                    continue;
                }
                if(initialAppearancesCount[e1] > maxAppearanceCount) {
                    if(debug) {
                        cout << "Discarded due to initial coverage on " << assemblyGraph[e1].id << endl;
                    }
                    continue;
                }
                if(finalAppearancesCount[e1] > maxAppearanceCount) {
                    if(debug) {
                        cout << "Discarded due to final coverage " << assemblyGraph[e1].id << endl;
                    }
                    continue;
                }

                if(jaccard(e0, e1, coverage) < minJaccard) {
                    if(debug) {
                        cout << "Discarded due to Jaccard." << endl;
                    }
                    continue;
                }

                if(debug) {
                    cout << "Adding edge " << assemblyGraph[e0].id << " " << assemblyGraph[e1].id << endl;
                }
                add_edge(vertexMap[e0], vertexMap[e1], coverage, graph);
            }
        }
    }
}



// In the graph, find a path that starts at a given AEdge
// and moves forward/backward. At each step we choose the child vertex
// corresponding to the longest AEdge.
void ReadFollowing::findPath(AEdge ae, uint64_t direction, vector<vertex_descriptor>& path) const
{
    if(direction == 0) {
        findForwardPath(ae, path);
    } else if(direction == 1) {
        findBackwardPath(ae, path);
    } else {
        SHASTA_ASSERT(0);
    }
}



// Test version.
void ReadFollowing::findForwardPath(AEdge ae, vector<vertex_descriptor>& path) const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    // Find the start vertex.
    const auto it = vertexMap.find(ae);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Each iteration adds one vertex to the path.
    path.clear();
    path.push_back(v);
    while(true) {

         // Find the best next vertex.
         vertex_descriptor vNext = null_vertex();
         double bestJaccardSum = 0.;
         BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
             vertex_descriptor v1 = target(e, graph);

             // Compute the sum of Jaccard similarities between all
             // previous vertices in the path and v1.
             double jaccardSum = 0.;
             for(const vertex_descriptor vPrevious: path) {
                 edge_descriptor ePrevious;
                 bool edgeWasFound;
                 tie(ePrevious, edgeWasFound) = edge(vPrevious, v1, graph);
                 if(edgeWasFound) {
                    jaccardSum += jaccard(ePrevious);
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
         if(graph[v].isLong) {
             break;
         }
     }

    /*
     cout << "Path:" << endl;
     for(const vertex_descriptor v: path) {
         cout << assemblyGraph[graph[v].ae].id << ",";
     }
     cout << endl;
     */
}



void ReadFollowing::findBackwardPath(AEdge ae, vector<vertex_descriptor>& path) const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    // Find the start vertex.
    const auto it = vertexMap.find(ae);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    // Each iteration adds one vertex to the path.
    path.clear();
    path.push_back(v);
    while(true) {

         // Find the best next vertex.
         vertex_descriptor vNext = null_vertex();
         double bestJaccardSum = 0.;
         BGL_FORALL_INEDGES(v, e, graph, Graph) {
             vertex_descriptor v1 = source(e, graph);

             // Compute the sum of Jaccard similarities between all
             // previous vertices in the path and v1.
             double jaccardSum = 0.;
             for(const vertex_descriptor vPrevious: path) {
                 edge_descriptor ePrevious;
                 bool edgeWasFound;
                 tie(ePrevious, edgeWasFound) = edge(v1, vPrevious, graph);
                 if(edgeWasFound) {
                    jaccardSum += jaccard(ePrevious);
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
         if(graph[v].isLong) {
             break;
         }
     }
     std::ranges::reverse(path);

     /*
     cout << "Path:" << endl;
     for(const vertex_descriptor v: path) {
         cout << assemblyGraph[graph[v].ae].id << ",";
     }
     cout << endl;
     */
}



// Old version.
#if 0
void ReadFollowing::findForwardPath(AEdge ae) const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    const auto it = vertexMap.find(ae);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    vector<AEdge> path;
    while(true) {
        path.push_back(graph[v].ae);

        vertex_descriptor vNext = null_vertex();
        // uint64_t bestLength = 0;
        double bestJaccard = 0.;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            vertex_descriptor v1 = target(e, graph);
            // const uint64_t length = graph[v1].length;
            const double j = jaccard(e);
            if(j > bestJaccard /*length > bestLength */) {
                vNext = v1;
                // bestLength = length;
                bestJaccard = j;
            }
        }

        if(vNext == null_vertex()) {
            break;
        }

        v = vNext;
    }

    cout << "Path:" << endl;
    for(const AEdge ae: path) {
        cout << assemblyGraph[ae].id << ",";
    }
    cout << endl;
}



void ReadFollowing::findBackwardPath(AEdge ae) const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    const auto it = vertexMap.find(ae);
    SHASTA_ASSERT(it != vertexMap.end());
    vertex_descriptor v = it->second;

    vector<AEdge> path;
    while(true) {
        path.push_back(graph[v].ae);

        vertex_descriptor vNext = null_vertex();
        // uint64_t bestLength = 0;
        double bestJaccard = 0.;
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            vertex_descriptor v1 = source(e, graph);
            // const uint64_t length = graph[v1].length;
            const double j = jaccard(e);
            if(j > bestJaccard /*length > bestLength */) {
                vNext = v1;
                // bestLength = length;
                bestJaccard = j;
            }
        }

        if(vNext == null_vertex()) {
            break;
        }

        v = vNext;
    }
    std::ranges::reverse(path);

    cout << "Path:" << endl;
    for(const AEdge ae: path) {
        cout << assemblyGraph[ae].id << ",";
    }
    cout << endl;
}
#endif



ReadFollowingVertex::ReadFollowingVertex(
    const AssemblyGraph& assemblyGraph,
    AssemblyGraph::edge_descriptor ae,
    uint64_t longLengthThreshold) :
    ae(ae)
{
    const AssemblyGraphEdge& edge = assemblyGraph[ae];
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }

    isLong = (length >= longLengthThreshold);
}



uint64_t ReadFollowing::getInitialAppearancesCount(AssemblyGraph::edge_descriptor ae) const
{
    const auto it = initialAppearancesCount.find(ae);
    if(it == initialAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}



uint64_t ReadFollowing::getFinalAppearancesCount(AEdge ae) const
{
    const auto it = finalAppearancesCount.find(ae);
    if(it == finalAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}



// Jaccard similarity for an EdgePairsGraph.
// Computed using the finalAppearancesCount of the surce vertex
// and the initialAppearancesCount of the target vertex.
double ReadFollowing::jaccard(edge_descriptor e) const
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;

    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);

    const AEdge ae0 = graph[v0].ae;
    const AEdge ae1 = graph[v1].ae;

    const uint64_t n0 = getFinalAppearancesCount(ae0);
    const uint64_t n1 = getInitialAppearancesCount(ae1);

    const uint64_t n01 = graph[e].coverage;

    const uint64_t intersectionSize = n01;
    const uint64_t unionSize = n0 + n1 - intersectionSize;

    return double(intersectionSize) / double(unionSize);
}



double ReadFollowing::jaccard(AEdge ae0, AEdge ae1, uint64_t coverage) const
{
    const uint64_t n0 = getFinalAppearancesCount(ae0);
    const uint64_t n1 = getInitialAppearancesCount(ae1);

    const uint64_t n01 = coverage;

    const uint64_t intersectionSize = n01;
    const uint64_t unionSize = n0 + n1 - intersectionSize;

    return double(intersectionSize) / double(unionSize);

}




void ReadFollowing::findAssemblyPaths()
{
    using Graph = ReadFollowing;
    const Graph& graph = *this;
    const bool debug = true;

    // Gather the vertices flagged as long.
    vector<vertex_descriptor> longVertices;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(graph[v].isLong) {
            longVertices.push_back(v);
        }
    }
    if(debug) {
        cout << "Found the following " << longVertices.size() <<
            " long assembly graph edges:";
        for(const vertex_descriptor v: longVertices) {
            cout << " " << assemblyGraph[graph[v].ae].id;
        }
        cout << endl;
    }



    // Call findPath in both directions for each of the longVertices.
    // Store all paths found, keyed by the first and last vertex.
    using Path = vector<vertex_descriptor>;
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<Path> > pathTable;
    vector<vertex_descriptor> path;
    for(const vertex_descriptor v0: longVertices) {

        // Forward.
        findPath(graph[v0].ae, 0, path);
        if(path.size() > 1) {
            const vertex_descriptor v1 = path.back();
            if(graph[v1].isLong) {
                pathTable[make_pair(v0, v1)].push_back(path);
            }
        }

        // Backward.
        findPath(graph[v0].ae, 1, path);
        if(path.size() > 1) {
            const vertex_descriptor v1 = path.front();
            if(graph[v1].isLong) {
                pathTable[make_pair(v1, v0)].push_back(path);
            }
        }
    }



    // The ones that were found twice (that is, in both directions) will
    // be used to create assembly paths.
    assemblyPaths.clear();
    uint64_t debugOutputIndex = 0;
    for(const auto& p: pathTable) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        const vector<Path>& paths = p.second;
        if(debug) {
            cout << "Found " << paths.size() << " paths between " << id(v0) << " and " << id(v1) << endl;
        }

        if(paths.size() != 2) {
            continue;
        }

        if(debug) {
            cout << "Looking for an optimal path between " << id(v0) << " and " << id(v1) << endl;
            for(const Path& path: paths) {
                cout << "One-directional path of length " << path.size() << endl;
                for(const vertex_descriptor v: path) {
                    cout << id(v) << ",";
                }
                cout << endl;
            }
        }



        // Create a filtered graph containing only the vertices found in these two paths.
        class Filter {
        public:
            bool operator()(const vertex_descriptor& v) const
            {
                return vertices.contains(v);
            }
            bool operator()(const edge_descriptor& e) const
            {
                const vertex_descriptor v0 = source(e, *graph);
                const vertex_descriptor v1 = target(e, *graph);
                return vertices.contains(v0) and vertices.contains(v1);
            }
            Filter(const Graph* graph = 0) : graph(graph) {}
            const Graph* graph;
            std::set<vertex_descriptor> vertices;
        };
        Filter filter(&graph);
        for(const Path& path: paths) {
            for(const vertex_descriptor v: path) {
                filter.vertices.insert(v);
            }
        }
        if(debug) {
            cout << "Using " << filter.vertices.size() <<
                " vertices to compute an optimal path." << endl;
        }
        using FilteredGraph = boost::filtered_graph<Graph, Filter>;
        FilteredGraph filteredGraph(graph, filter);

        if(debug) {
            ofstream dot("ReadFollowing-FilteredGraph-" + to_string(debugOutputIndex) + ".dot");
            dot << "digraph FilteredGraph {\n";
            BGL_FORALL_EDGES(e, filteredGraph, FilteredGraph) {
                const vertex_descriptor v0 = source(e, filteredGraph);
                const vertex_descriptor v1 = target(e, filteredGraph);
                dot << id(v0) << "->" << id(v1) << ";\n";
            }
            dot << "}\n";
        }

        // Find the longest path in this filtered graph.
        vector<edge_descriptor> edgePath;
        longestPath(filteredGraph, edgePath);
        vector<vertex_descriptor> optimalPath;
        optimalPath.push_back(source(edgePath.front(), graph));
        for(const edge_descriptor e: edgePath) {
            optimalPath.push_back(target(e, graph));
        }

        // Turn it into a sequence of AEdges and store it.
        assemblyPaths.emplace_back();
        auto& assemblyPath = assemblyPaths.back();
        for(const vertex_descriptor v: optimalPath) {
            assemblyPath.push_back(graph[v].ae);
        }

        if(debug) {
            cout << "The assembly path has " << edgePath.size() + 1 << " vertices:" << endl;
            for(const AEdge ae: assemblyPath) {
                cout << assemblyGraph[ae].id << ",";
            }
            cout << endl;

            // Also write a fasta file numbering edges along this path.
            ofstream fasta("ReadFollowing-OptimalPath-" + to_string(debugOutputIndex) + ".fasta");
            for(uint64_t i=0; i<assemblyPath.size(); i++) {
                const AEdge ae = assemblyPath[i];
                const AssemblyGraphEdge& aEdge = assemblyGraph[ae];
                vector<shasta::Base> sequence;
                aEdge.getSequence(sequence);
                fasta << ">" << i << " " << assemblyGraph[ae].id << " " << sequence.size() << endl;
                std::ranges::copy(sequence, ostream_iterator<shasta::Base>(fasta));
                fasta << endl;
            }

        }


        ++debugOutputIndex;
    }
}




