#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "ReadFollowing.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>


// Constructor from another AssemblyGraph and assembly paths
// like the ones computed by ReadFollowing.
// (Those are not actually paths in the AssemblyGraph but simply
// sequences of edges).
AssemblyGraph::AssemblyGraph(
    const Assembler& assembler,
    const Options& options,
    const AssemblyGraph& oldAssemblyGraph,
    const vector< vector<edge_descriptor> >& assemblyPaths) :
    MappedMemoryOwner(assembler),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(assembler.anchors()),
    journeys(assembler.journeys()),
    options(options),
    orderById(*this)
{
    AssemblyGraph& newAssemblyGraph = *this;

    // The terminal edges of the paths are special in that
    // only one copy of them is carried in the new AssemblyGraph.

    // Gather the terminal edges.
    vector<edge_descriptor> terminalEdges;
    for(const vector<edge_descriptor>& path: assemblyPaths) {
        terminalEdges.push_back(path.front());
        terminalEdges.push_back(path.back());
    }
    deduplicate(terminalEdges);
    sort(terminalEdges.begin(), terminalEdges.end(), oldAssemblyGraph.orderById);



    // Generate an edge in the newAssemblyGraph for each terminal
    // edge in the old assembly graph.
    // Keep track of them, keyed by edge_descriptor in the oldAssemblyGraph.
    std::map<edge_descriptor, edge_descriptor> terminalEdgesMap;
    for(const edge_descriptor eOld: terminalEdges) {
        const AssemblyGraphEdge& oldEdge = oldAssemblyGraph[eOld];
        const AnchorId anchorId0 = oldEdge.front().anchorPair.anchorIdA;
        const AnchorId anchorId1 = oldEdge.back().anchorPair.anchorIdB;
        const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
        const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(oldEdge.id), newAssemblyGraph);
        AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
        newEdge = oldEdge;
        newEdge.id = newAssemblyGraph.nextEdgeId++;
        terminalEdgesMap.insert(make_pair(eOld, eNew));
    }



    // Now generate the remaining edges for each path.
    // Note that non-terminal edges can get more than one copy
    // in the newAssemblyGraph.
    vector< vector<edge_descriptor> > newAssemblyPaths;
    for(const auto& oldAssemblyPath: assemblyPaths) {
        newAssemblyPaths.emplace_back();
        auto& newAssemblyPath = newAssemblyPaths.back();

        for(const edge_descriptor eOld: oldAssemblyPath) {
            const auto it = terminalEdgesMap.find(eOld);
            if(it != terminalEdgesMap.end()) {
                newAssemblyPath.push_back(it->second);
            } else {

                // This is not a terminal edge and we have to create.
                const AssemblyGraphEdge& oldEdge = oldAssemblyGraph[eOld];
                const AnchorId anchorId0 = oldEdge.front().anchorPair.anchorIdA;
                const AnchorId anchorId1 = oldEdge.back().anchorPair.anchorIdB;
                const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
                const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
                edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(oldEdge.id), newAssemblyGraph);
                AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
                newEdge = oldEdge;
                newEdge.id = newAssemblyGraph.nextEdgeId++;
                newAssemblyPath.push_back(eNew);
            }
        }
    }



    // For each pair of consecutive edges in a path (e0, e1),
    // generate a new edge in-between to bridge between them.
    // The code is similar to Tangle1::addConnectPair and Tangle1::detangle,
    // but simpler.
    for(const auto& newAssemblyPath: newAssemblyPaths) {

        for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const edge_descriptor e0 = newAssemblyPath[i0];
            const edge_descriptor e1 = newAssemblyPath[i1];

            const vertex_descriptor v0 = target(e0, newAssemblyGraph);
            const vertex_descriptor v1 = source(e1, newAssemblyGraph);

            const AnchorId anchorId0 = newAssemblyGraph[v0].anchorId;
            const AnchorId anchorId1 = newAssemblyGraph[v1].anchorId;

            // Create the new edge.
            // If the two anchors are the same, leave it empty without any steps.
            // Otherwise use the same process in Tangle1::addConnectPair.
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(newAssemblyGraph.nextEdgeId++), newAssemblyGraph);
            AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
            if(anchorId0 != anchorId1) {

                // Create the RestrictedAnchorGraph, then:
                // - Remove vertices not accessible from anchorId0 and anchorId1.
                // - Remove cycles.
                // - Find the longest path.
                // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

                ostream html(0);
                const TangleMatrix1 tangleMatrix(
                    newAssemblyGraph,
                    vector<edge_descriptor>(1, e0),
                    vector<edge_descriptor>(1, e1),
                    html);

                RestrictedAnchorGraph restrictedAnchorGraph(
                    newAssemblyGraph.anchors, newAssemblyGraph.journeys, tangleMatrix, 0, 0, html);
                restrictedAnchorGraph.removeLowCoverageEdges(anchorId0, anchorId1);
                restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
                restrictedAnchorGraph.removeCycles();
                restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
                vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                // restrictedAnchorGraph.findLongestPath(longestPath);
                restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                    const auto& rEdge = restrictedAnchorGraph[re];
                    if(rEdge.anchorPair.size() < newAssemblyGraph.options.detangleMinCoverage) {
                        newEdge.clear();
                        SHASTA_ASSERT(0);
                    }
                    newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
                }
            }
        }
    }

    newAssemblyGraph.removeEmptyEdges();
    newAssemblyGraph.compress();
}



// Note assemlbyPaths are not necessarily paths in the AssemblyGraph.
// There may be jumps, which are bridged using local assemblies.
void AssemblyGraph::findAssemblyPaths(vector< vector<edge_descriptor> >& assemblyPaths) const
{
	ReadFollowing::Graph readFollowingGraph(*this);
	readFollowingGraph.findPaths(assemblyPaths);
}



// Note assemlbyPaths are not necessarily paths in the AssemblyGraph.
// There may be jumps, which are bridged using local assemblies.
// The assembly paths must satisfy the following rules:
// - They must consist of at least 2 edge_descriptors.
// - If an edge_descriptor appears at the beginning or end of an assembly path,
//   it cannot appear elsewhere, in other paths or in the same path.
// - Edges that appear inside assembly paths can appear multiple times
//   in the same assembly paths or different assembly paths.
// These rules are satisfied when the assembly paths are computed using
// findAssemblyPaths.
void AssemblyGraph::connectAssemblyPaths(const vector< vector<edge_descriptor> >&  assemblyPaths)
{
	AssemblyGraph& assemblyGraph = *this;

	const bool debug = true;
	if(debug) {
		cout << "Connecting " << assemblyPaths.size() << " assembly paths:" << endl;
		for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
			SHASTA_ASSERT(assemblyPath.size() > 1);
			const edge_descriptor e0 = assemblyPath.front();
			const edge_descriptor e1 = assemblyPath.back();
			cout << "Assembly path with " << assemblyPath.size() << " segments beginning at " <<
				assemblyGraph[e0].id << " and ending at " <<
				assemblyGraph[e1].id << endl;
		}
	}

	// Gather edge_descriptors that appear at the beginning/end of assembly paths.
	std::set<edge_descriptor> pathInitialSegments;
	std::set<edge_descriptor> pathFinalSegments;
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
		SHASTA_ASSERT(assemblyPath.size() > 1);
		const edge_descriptor e0 = assemblyPath.front();
		const edge_descriptor e1 = assemblyPath.back();
		SHASTA_ASSERT(e0 != e1);
		SHASTA_ASSERT(not pathInitialSegments.contains(e0));
		SHASTA_ASSERT(not pathFinalSegments.contains(e0));
		SHASTA_ASSERT(not pathInitialSegments.contains(e1));
		SHASTA_ASSERT(not pathFinalSegments.contains(e1));
		pathInitialSegments.insert(e0);
		pathFinalSegments.insert(e1);
	}

	// Gather edge_descriptors that appear internally to assembly paths.
	// Count how many times each of them appear.
	// The ones that appear only once will keep their id.
	std::map<edge_descriptor, uint64_t> pathInternalSegments;
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
		for(uint64_t i=1; i<assemblyPath.size()-1; i++) {
			const edge_descriptor e = assemblyPath[i];
			SHASTA_ASSERT(not pathInitialSegments.contains(e));
			SHASTA_ASSERT(not pathFinalSegments.contains(e));
			const auto it = pathInternalSegments.find(e);
			if(it == pathInternalSegments.end()) {
				pathInternalSegments.insert({e, 1});
			} else {
				++it->second;
			}
		}
	}

	if(debug) {
		for(const auto& p: pathInternalSegments) {
			if(p.second > 1) {
				const edge_descriptor e = p.first;
				cout << "Segment " << assemblyGraph[e].id <<
					" appears more than once in assembly paths." << endl;
			}
		}
	}



	// Each assembly path generates a linear sequence of new edges that
	// can be later collapsed into a single edge by calling compress.
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {

		// Generate the new edges.
		vector<edge_descriptor> newAssemblyPath;
		for(uint64_t i=0; i<assemblyPath.size(); i++) {
			const edge_descriptor e = assemblyPath[i];
			const AssemblyGraphEdge& edge = assemblyGraph[e];
			const vertex_descriptor v0 = source(e, assemblyGraph);
			const vertex_descriptor v1 = target(e, assemblyGraph);
			const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
			const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

			// Define the source and target vertices of the new edge, and its id.
			vertex_descriptor v0New = null_vertex();
			vertex_descriptor v1New = null_vertex();
			uint64_t idNew = invalid<uint64_t>;
			if(i == 0) {
				// Initial segment.
				v0New = v0;
				v1New = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
				idNew = edge.id;
			} else if(i == assemblyPath.size() - 1) {
				// Final segment.
				v0New = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
				v1New = v1;
				idNew = edge.id;
			} else {
				// Internal segment.
				v0New = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
				v1New = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
				if(pathInternalSegments[e] == 1) {
					idNew = edge.id;
				} else {
					idNew = nextEdgeId++;
				}
			}

			// Create the new edge.
	        edge_descriptor eNew;
	        tie(eNew, ignore) = add_edge(v0New, v1New, AssemblyGraphEdge(idNew), assemblyGraph);
	        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
	        edgeNew = edge;
	        edgeNew.id = idNew;
	        newAssemblyPath.push_back(eNew);
		}



	    // For each pair of consecutive edges in this path,
	    // generate a new edge in-between to bridge between them.
	    // The code is similar to Tangle1::addConnectPair and Tangle1::detangle,
	    // but simpler.
		for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
			const uint64_t i0 = i1 - 1;
			const edge_descriptor e0 = newAssemblyPath[i0];
			const edge_descriptor e1 = newAssemblyPath[i1];

			const vertex_descriptor v0 = target(e0, assemblyGraph);
			const vertex_descriptor v1 = source(e1, assemblyGraph);

			const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
			const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

			// Create the new edge.
			// If the two anchors are the same, leave it empty without any steps.
			// Otherwise use the same process in Tangle1::addConnectPair.
			edge_descriptor eNew;
			tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
			AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
			if(anchorId0 != anchorId1) {

				// Create the RestrictedAnchorGraph, then:
				// - Remove vertices not accessible from anchorId0 and anchorId1.
				// - Remove cycles.
				// - Find the longest path.
				// - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

				ostream html(0);
				const TangleMatrix1 tangleMatrix(
					assemblyGraph,
					vector<edge_descriptor>(1, e0),
					vector<edge_descriptor>(1, e1),
					html);

				RestrictedAnchorGraph restrictedAnchorGraph(anchors, journeys, tangleMatrix, 0, 0, html);
				restrictedAnchorGraph.removeLowCoverageEdges(anchorId0, anchorId1);
				restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
				restrictedAnchorGraph.removeCycles();
				restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
				vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
				// restrictedAnchorGraph.findLongestPath(longestPath);
				restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

				for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
					const auto& rEdge = restrictedAnchorGraph[re];
					if(rEdge.anchorPair.size() < options.detangleMinCoverage) {
						newEdge.clear();
						SHASTA_ASSERT(0);
					}
					newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
				}
			}
	    }
	}



	// Remove all edges that appear in one or more assembly paths.
	vector<edge_descriptor> edgesToBeRemoved;
	BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
		if(
			pathInitialSegments.contains(e) or
			pathFinalSegments.contains(e) or
			pathInternalSegments.contains(e)) {
			edgesToBeRemoved.push_back(e);
		}
	}
	for(const edge_descriptor e: edgesToBeRemoved) {
		boost::remove_edge(e, assemblyGraph);
	}

	// Compress the lienar chains we created.
	compress();
}



void AssemblyGraph::findAndConnectAssemblyPaths()
{
	vector< vector<edge_descriptor> > assemblyPaths;
	findAssemblyPaths(assemblyPaths);
	connectAssemblyPaths(assemblyPaths);
}

