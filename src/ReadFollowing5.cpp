// Shasta2.
#include "ReadFollowing5.hpp"
#include "deduplicate.hpp"
#include "Options.hpp"
#include "SegmentStepSupport.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "fstream.hpp"



Graph::Graph(
    const AssemblyGraph& assemblyGraph,
    uint64_t tangleId,
    const Tangle& tangle) :
    assemblyGraph(assemblyGraph),
    tangle(tangle),
    tangleId(tangleId),
    isSelfComplementaryTangle(tangle.isSelfComplementary())
{
    Graph& graph = *this;
    const bool debug = true;

    createVertices();
    createEdges();

    if(debug) {
        cout << "The read following graph for this tangle has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
        writeGraphviz("ReadFollowingGraph-Tangle-" + to_string(tangleId) + ".dot");
    }

    findShortestPaths();
}



void Graph::createVertices()
{
    Graph& graph = *this;

    for(const Segment segment: tangle.tangleEdges) {
        add_vertex(Vertex(assemblyGraph, segment, false, false), graph);
    }
    for(const Segment segment: tangle.entrances) {
        entranceVertices.push_back(add_vertex(Vertex(assemblyGraph, segment, true, false), graph));
    }
    for(const Segment segment: tangle.exits) {
        exitVertices.push_back(add_vertex(Vertex(assemblyGraph, segment, false, true), graph));
    }
}



Vertex::Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment,
    bool isEntrance,
    bool isExit) :
    segment(segment),
    isEntrance(isEntrance),
    isExit(isExit)
{
    const uint32_t representativeRegionStepCount = uint32_t(assemblyGraph.options.representativeRegionStepCount);
    vector<SegmentStepSupport> support;

    // Fill in initial support.
    // Don't include ReadIds that appear on both strands.
    if(not isEntrance) {
        SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, support);
        initialSupport.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            initialSupport.push_back(orientedReadId);
        }
    }

    // Fill in final support.
    // Don't include ReadIds that appear on both strands.
    if(not isExit) {
        SegmentStepSupport::getFinalLast(assemblyGraph, segment, representativeRegionStepCount, support);
        finalSupport.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            finalSupport.push_back(orientedReadId);
        }
    }
}



void Graph::createEdges()
{
    Graph& graph = *this;
    const bool debug = false;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    const double a = 3.;                // dB
    const double b = 15.;               // dB
    // const double logPThreshold = 10.;   // dB

    // No html output from analyzeSegmentPair.
    ostream html(0);



    // For each OrientedReadId, gather the vertices that have that OrientedReadId
    // in their initialSupport/finalSupport.
    std::map<OrientedReadId, vector<vertex_descriptor> > initialSupportMap;
    std::map<OrientedReadId, vector<vertex_descriptor> > finalSupportMap;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        for(const OrientedReadId orientedReadId: vertex.initialSupport) {
            initialSupportMap[orientedReadId].push_back(v);
        }
        for(const OrientedReadId orientedReadId: vertex.finalSupport) {
            finalSupportMap[orientedReadId].push_back(v);
        }
    }

    // Gather vertex pairs (v0, v1) that have common OrientedReadIds
    // between the final support of v0 and the intial support of v1.
    vector< pair<vertex_descriptor, vertex_descriptor> > vertexPairs;
    for(const auto&[orientedReadId, finalSupportVertices]: finalSupportMap) {
        const auto& initialSupportVertices = initialSupportMap[orientedReadId];
        for(const vertex_descriptor v0: finalSupportVertices) {
            for(const vertex_descriptor v1: initialSupportVertices) {
                vertexPairs.push_back(make_pair(v0, v1));
            }
        }
    }

    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(vertexPairs, count, minCommonCount);
    SHASTA2_ASSERT(vertexPairs.size() == count.size());



    // Each vertex pair can generate an edge, but we have to do it in a strand-symmetric way.
    vector<OrientedReadId> commonOrientedReadIds;
    for(const auto&[v0, v1]: vertexPairs) {
        if(v0 == v1) {
            continue;
        }
        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];
        const Segment segment0 = vertex0.segment;
        const Segment segment1 = vertex1.segment;

        // If this is a self-complementary tangle, generate edges in pairs.
        // See if this a vertex pair that we want to use (to ensure strand symmetry).
        if(isSelfComplementaryTangle) {
            const uint64_t id0 = assemblyGraph[segment0].id;
            const uint64_t id1 = assemblyGraph[segment1].id;
            const Segment segment0Rc = assemblyGraph[segment0].eRc;
            const Segment segment1Rc = assemblyGraph[segment1].eRc;
            SHASTA2_ASSERT(segment0Rc != assemblyGraphNullEdge);
            SHASTA2_ASSERT(segment1Rc != assemblyGraphNullEdge);
            SHASTA2_ASSERT(segment0Rc != segment0);
            SHASTA2_ASSERT(segment0Rc != segment1);
            SHASTA2_ASSERT(segment1Rc != segment0);
            SHASTA2_ASSERT(segment1Rc != segment1);
            const uint64_t id0Rc = assemblyGraph[segment0Rc].id;
            const uint64_t id1Rc = assemblyGraph[segment1Rc].id;
            if(min(id0, id1) >= min(id0Rc, id1Rc)) {
                continue;
            }
        }

        const SegmentPairInformation segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
            html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

        if(debug) {
            cout << "Edge candidate " << assemblyGraph.id(segment0) << " " <<
                assemblyGraph.id(segment1) << ", offset " << segmentPairInformation.segmentOffset <<
                ", common count " << segmentPairInformation.commonCount <<
                ", missing " << segmentPairInformation.missing() << endl;
        }

        // If the offset is negative, discard it.
        if(segmentPairInformation.segmentOffset < 0) {
            if(debug) {
                cout << "Discarded due to negative offset." << endl;
            }
            continue;
        }

        // If it does not satisfy our requirements, discard it.
        if(segmentPairInformation.commonCount < minCommonCount) {
            if(debug) {
                cout << "Discarded due low common count." << endl;
            }
            continue;
        }
        const double logP = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing());
        /*
        if(logP < logPThreshold) {
            if(debug) {
                cout << "Discarded due low logP." << endl;
            }
            continue;
        }
        */
        if(not assemblyGraph.canConnect(segment0, segment1)) {
            if(debug) {
                cout << "Discarded because cannot connect." << endl;
            }
            continue;
        }


        if(debug) {
            cout << "Keeping " << assemblyGraph.id(segment0) << " " << assemblyGraph.id(segment1) << " " << logP << endl;
        }
        add_edge(v0, v1, Edge(segmentPairInformation.commonCount, logP), graph);
        if(isSelfComplementaryTangle) {
            const Segment segment0Rc = assemblyGraph[segment0].eRc;
            const Segment segment1Rc = assemblyGraph[segment1].eRc;
            if(debug) {
                cout << "Keeping " << assemblyGraph.id(segment1Rc) << " " << assemblyGraph.id(segment0Rc) << " " << logP << endl;
            }
            SHASTA2_ASSERT(0);
            // Must add the reverse complement edge.
        }
    }
}



void Graph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void Graph::writeGraphviz(ostream& dot) const
{
    const Graph& graph = *this;
    dot << "digraph ReadFollowingGraph {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        dot << assemblyGraph.id(graph[v].segment);
        if(graph[v].isEntrance) {
            dot << "[style=filled fillcolor=pink]";
        }
        if(graph[v].isExit) {
            dot << "[style=filled fillcolor=green]";
        }
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << assemblyGraph.id(graph[v0].segment) << "->" <<
            assemblyGraph.id(graph[v1].segment) <<
            "[ label=\"" << int64_t(std::round(graph[e].logP)) << "\"]"
            "\n";
    }

    dot << "}\n";
}


void Graph::findShortestPaths() const
{
    const Graph& graph = *this;
    using boost::make_assoc_property_map;

    // Create a vertex index map, needed below.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Loop over entrances.
    for(uint64_t iEntrance=0; iEntrance<entranceVertices.size(); iEntrance++) {
        const vertex_descriptor vEntrance = entranceVertices[iEntrance];

        // Create a shortest path tree rooted at this entrance.
        std::map<vertex_descriptor, double> distanceMap;
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
        dijkstra_shortest_paths(graph, vEntrance,
           weight_map(boost::get(&Edge::weight, graph)).
           vertex_index_map(make_assoc_property_map(vertexIndexMap)).
           distance_map(make_assoc_property_map(distanceMap)).
           predecessor_map(make_assoc_property_map(predecessorMap)));

        // Write out the distances of the exits from this entrance.
        for(uint64_t iExit=0; iExit<exitVertices.size(); iExit++) {
            const vertex_descriptor vExit = exitVertices[iExit];
            const double distance = distanceMap.at(vExit);
            cout << iEntrance << " " << iExit << " " << distance << endl;
        }

    }
}
