// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "ExactDetangler.hpp"
#include "findLinearChains.hpp"
#include "findConvergingVertex.hpp"
#include "color.hpp"
#include "inducedSubgraphIsomorphisms.hpp"
#include "Journeys.hpp"
#include "LikelihoodRatioDetangler.hpp"
#include "LocalAssembly2.hpp"
#include "MurmurHash2.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "rle.hpp"
#include "SearchGraph.hpp"
#include "SimpleDetangler.hpp"
#include "Superbubble.hpp"
#include "SuperbubbleChain.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
#include "TrivialDetangler.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph>;



// Initial construction from the AnchorGraph.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const AnchorGraph& anchorGraph,
    const Options& options) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options)
{
    AssemblyGraph& assemblyGraph = *this;

    // Create a filtered AnchorGraph containing only the edges marked as "useForAssembly".
    class EdgePredicate {
    public:
        bool operator()(const AnchorGraph::edge_descriptor& e) const
        {
            return (*anchorGraph)[e].useForAssembly;
        }
        EdgePredicate(const AnchorGraph& anchorGraph) : anchorGraph(&anchorGraph) {}
        EdgePredicate() : anchorGraph(0) {}
        const AnchorGraph* anchorGraph;
    };
    using FilteredAnchorGraph = boost::filtered_graph<AnchorGraph, EdgePredicate>;
    FilteredAnchorGraph filteredAnchorGraph(anchorGraph, EdgePredicate(anchorGraph));



    // Find linear chains of edges in the FilteredAnchorGraph.
    vector< std::list<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(filteredAnchorGraph, 1, chains);

    // Generate vertices.
    // At this stage there is a vertex for each AnchorGraph vertex
    // that is at the beginning or end of a linear chain,
    // so there is only one vertex for a given AnchorId.
    // However, after detangling there can be more than one vertex
    // with a given AnchorId. So the vertexMap is only used in this constructor.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        if(not vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(not vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    SHASTA_ASSERT(vertexMap.size() == num_vertices(assemblyGraph));



    // Generate the edges. There is an edge for each linear chain.
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Each AnchorGraph edge in the chain contributes a step to this AssemblyGraph edge.
        for(const AnchorGraph::edge_descriptor eA: chain) {
            const AnchorGraphEdge& edgeA = anchorGraph[eA];
            edge.emplace_back(edgeA.anchorPair, edgeA.offset);
        }
    }


    check();
}



// Deserialize constructor.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const Options& options,
    const string& stage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options)
{
    load(stage);
}



// Detangle, phase, assemble sequence, output.
void AssemblyGraph::run(uint64_t threadCount)
{
    // AssemblyGraph& assemblyGraph = *this;
    // const uint64_t detangleMaxIterationCount = 10;

    // Initial output.
    write("A");

    // Prune.
    prune();
    compress();

    // Simplify Superbubbles and remove or simplify bubbles likely caused by errors.
    simplifySuperbubbles();
    write("B");
    bubbleCleanup(threadCount);
    compress();
    write("C");

    // Phase SuperbubbleChains.
    phaseSuperbubbleChains();
    write("D");

#if 0
    // Detangling.
    LikelihoodRatioDetangler detangler(
        options.detangleMinCommonCoverage,
        options.detangleEpsilon,
        options.detangleMaxLogP,
        options.detangleMinLogPDelta);
    detangle(detangleMaxIterationCount, std::numeric_limits<uint64_t>::max(), detangler);
    write("E");
#endif

    // Sequence assembly.
    assembleAll(threadCount);
    write("Z");
    writeFasta("Z");
}



void AssemblyGraph::check() const
{
    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());



        // Check that the first/last AnchorIds of this edge are consistent
        // with the ones in the source/target vertices.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

        SHASTA_ASSERT(edge.front().anchorPair.anchorIdA == anchorId0);
        SHASTA_ASSERT(edge.back().anchorPair.anchorIdB == anchorId1);



        // Check that AnchorPairs in this edge are adjacent to each other.
        for(uint64_t i1=1; i1<edge.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA_ASSERT(edge[i0].anchorPair.anchorIdB == edge[i1].anchorPair.anchorIdA);
        }    }

}



uint64_t AssemblyGraphEdge::offset() const
{
    uint64_t sum = 0;
    for(const auto& step: *this) {
        sum += step.offset;
    }
    return sum;
}



void AssemblyGraph::write(const string& stage)
{
    cout << "Assembly graph at stage " << stage << " has " <<
        num_vertices(*this) << " vertices and " <<
        num_edges(*this) << " edges (segments)." << endl;

    save(stage);
    writeGfa("AssemblyGraph-" + stage + ".gfa");
    writeGraphviz("AssemblyGraph-" + stage + ".dot");
}



void AssemblyGraph::writeFasta(const string& stage) const
{
    const AssemblyGraph& assemblyGraph = *this;;

    ofstream fasta("AssemblyGraph-" + stage + ".fasta");

    vector<shasta::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.getSequence(sequence);

        fasta << ">" << edge.id << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(fasta));
        fasta << "\n";
    }
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}



void AssemblyGraph::writeGfa(ostream& gfa) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Each edge generates a gfa segment.
    vector<shasta::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(gfa));
            gfa << "\t";
            gfa << "LN:i:" << sequence.size() << "\n";

        } else {
            gfa << "*\t";
            gfa << "LN:i:" << edge.offset() << "\n";
        }
    }



    // For each vertex, generate a link between each pair of
    // incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;

                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }


}



void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AssemblyGraph::writeGraphviz(ostream& dot) const
{
    const AssemblyGraph& assemblyGraph = *this;

    dot << "digraph AssemblyGraph {\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
    	const AssemblyGraphVertex& vertex = assemblyGraph[v];
    	dot <<
    		vertex.id <<
    		" [label=\"" << anchorIdToString(vertex.anchorId) << "\\n" << vertex.id << "\"]"
    	    ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
    	const vertex_descriptor v0 = source(e, assemblyGraph);
    	const vertex_descriptor v1 = target(e, assemblyGraph);
    	const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
    	const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];
    	dot <<
    	    vertex0.id << "->" <<
    	    vertex1.id <<
    	    " [label=\"" << edge.id << "\\n" <<
    	    (edge.wasAssembled ? edge.sequenceLength() : edge.offset()) <<
    	    "\\n" << edge.size() <<
    	    "\"]"
    	    ";\n";
    }

    dot << "}\n";
}



// Assemble sequence for all edges.
void AssemblyGraph::assembleAll(uint64_t threadCount)
{
    performanceLog << timestamp << "Sequence assembly begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const AssemblyGraph& assemblyGraph = *this;

    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble(threadCount);
    edgesToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



// Assemble sequence for the specified edge.
void AssemblyGraph::assemble(edge_descriptor e, uint64_t threadCount)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble(threadCount);
}



// Assemble sequence for step i of the specified edge.
// This is the lowest level sequence assembly function and is not multithreaded.
// It runs a LocalAssembly2 on the AnchorPair for that step.
void AssemblyGraph::assembleStep(edge_descriptor e, uint64_t i)
{
    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphEdge& edge = assemblyGraph[e];
    AssemblyGraphEdgeStep& step = edge[i];

    if(step.anchorPair.anchorIdA == step.anchorPair.anchorIdB) {
        step.sequence.clear();
        return;
    }

    // Run the LocalAssembly2.
    ofstream html;  // Not open, so no html output takes place.
    LocalAssembly2 localAssembly(
        anchors, html, false,
        options.aDrift,
        options.bDrift,
        step.anchorPair);
    if(localAssembly.coverage() == 0) {
        throw runtime_error("No coverage for local assembly at assembly graph edge " +
            to_string(edge.id) + " step " + to_string(i));
    }
    localAssembly.run(false, options.localAssemblyOptions.maxAbpoaLength);
    localAssembly.getSequence(step.sequence);
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the stepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph::assemble(uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    stepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            stepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(stepsToBeAssembled.size(), batchCount);
    runThreads(&AssemblyGraph::assembleThreadFunction, threadCount);

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }

    edgesToBeAssembled.clear();
    stepsToBeAssembled.clear();
}



void AssemblyGraph::assembleThreadFunction(uint64_t /* threadId */)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {
            if((j % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << j << "/" << stepsToBeAssembled.size() << endl;
            }

            const auto& p = stepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            AssemblyGraphEdge& edge = assemblyGraph[e];
            SHASTA_ASSERT(i < edge.size());
            assembleStep(e, i);
        }
    }
}



// Clear sequence from all steps of all edges.
void AssemblyGraph::clearSequence()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.wasAssembled= false;
        for(AssemblyGraphEdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}



void AssemblyGraphEdge::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(const auto& step: *this) {
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}



uint64_t AssemblyGraphEdge::sequenceLength() const
{
    SHASTA_ASSERT(wasAssembled);

    uint64_t length = 0;
    for(const auto& step: *this) {
        length += step.sequence.size();
    }
    return length;
}



void AssemblyGraph::findBubbles(vector<Bubble>& bubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;
    bubbles.clear();

    // Look at bubbles with source v0.
    std::map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }

        for(const auto& p: m) {
            const vertex_descriptor v1 = p.first;
            const vector<edge_descriptor>& edges = p.second;
            if(edges.size() > 1) {
                Bubble bubble;
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                bubbles.push_back(bubble);
            }
        }
    }

    cout << "Found " << bubbles.size() << " bubbles." << endl;

    // Count the bubbles for each ploidy.
    vector<uint64_t> histogram;
    for(const Bubble& bubble: bubbles) {
        const uint64_t ploidy = bubble.edges.size();
        if(ploidy >= histogram.size()) {
            histogram.resize(ploidy+1, 0);
        }
        ++histogram[ploidy];
    }
    for(uint64_t ploidy=0; ploidy<histogram.size(); ploidy++) {
        const uint64_t frequency = histogram[ploidy];
        if(frequency) {
            cout << "Ploidy " << ploidy << ": " << frequency << " bubbles." << endl;
        }
    }

}



void AssemblyGraph::bubbleCleanup(uint64_t threadCount)
{
    while(bubbleCleanupIteration(threadCount) > 0);
}



uint64_t AssemblyGraph::bubbleCleanupIteration(uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    const uint64_t maxRepeatPeriod = 4;


    // Find all bubbles.
    vector<Bubble> allBubbles;
    findBubbles(allBubbles);



    // For each Bubble, compute the bridge AnchorPair we would use if
    // we were to remove the bubble. If the bridge AnchorPair
    // has low coverage, we can't remove the bubble.
    vector< pair<Bubble, AnchorPair> > candidateBubbles;
    for(const Bubble& bubble: allBubbles) {
        const vertex_descriptor v0 = bubble.v0;
        const vertex_descriptor v1 = bubble.v1;

        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
        const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];

        const AnchorId anchorId0 = vertex0.anchorId;
        const AnchorId anchorId1 = vertex1.anchorId;

        const AnchorPair bridgeAnchorPair(anchors, anchorId0, anchorId1, false);

        // If coverage of the bridgeAnchorPair is sufficient, add this bubble to our list of candidates.
        if(bridgeAnchorPair.orientedReadIds.size() >= options.bubbleCleanupMinCommonCount) {
            candidateBubbles.push_back(make_pair(bubble, bridgeAnchorPair));
        }


#if 0
        // This is the old code that assumes that the bubble is preceded/followed by a
        // single segment.
        SHASTA_ASSERT(in_degree(v0, assemblyGraph) == 1);
        SHASTA_ASSERT(out_degree(v1, assemblyGraph) == 1);

        in_edge_iterator it0;
        tie(it0, ignore) = in_edges(v0, assemblyGraph);
        const edge_descriptor e0 = *it0;
        const AnchorPair& anchorPair0 = assemblyGraph[e0].back().anchorPair;

        out_edge_iterator it1;
        tie(it1, ignore) = out_edges(v1, assemblyGraph);
        const edge_descriptor e1 = *it1;
        const AnchorPair& anchorPair1 = assemblyGraph[e1].front().anchorPair;

        // Construct the bridge AnchorPair.
        const AnchorPair bridgeAnchorPair = anchors.bridge(
            anchorPair0, anchorPair1,
            options.aDrift,
            options.bDrift);

        // If coverage of the bridgeAnchorPair is sufficient, add this bubble to our list of candidates.
        if(bridgeAnchorPair.orientedReadIds.size() >= options.assemblyGraphOptions.bubbleCleanupMinCommonCount) {
            candidateBubbles.push_back(make_pair(bubble, bridgeAnchorPair));
        }
#endif
    }
    cout << candidateBubbles.size() << " bubbles are candidate for removal." << endl;


    // Assemble sequence for all the edges of these bubbles.
    edgesToBeAssembled.clear();
    for(const auto& p: candidateBubbles) {
        const Bubble& bubble = p.first;
        for(const edge_descriptor e: bubble.edges) {
            if(not assemblyGraph[e].wasAssembled) {
                edgesToBeAssembled.push_back(e);
            }
        }
    }
    assemble(threadCount);



    // Now that we have sequence for the candidate bubbles, we can decide
    // which ones should be removed.
    uint64_t removedCount = 0;
    for(const auto& p: candidateBubbles) {
        const Bubble& bubble = p.first;

        // Gather the sequences of all the sides of this bubble
        vector< vector<shasta::Base> > sequences;
        for(const edge_descriptor e: bubble.edges) {
            sequences.emplace_back();
            SHASTA_ASSERT(assemblyGraph[e].wasAssembled);
            assemblyGraph[e].getSequence(sequences.back());
        }

        // This bubble can be removed if all the raw sequences are identical.
        bool allRawSequenceAreEqual = true;
        for(uint64_t i=1; i<sequences.size(); i++) {
            if(sequences[i] != sequences[0]) {
                allRawSequenceAreEqual = false;
                break;
            }
        }


        // If all raw sequence are equal, we can remove the bubble without checking the RLE sequences.
        // Otherwise we have to also check the RLE sequences.
        bool removeBubble = allRawSequenceAreEqual;
        if(not allRawSequenceAreEqual) {

            // Compute the RLE sequences.
            vector< vector<shasta::Base> > rleSequences;
            for(const vector<shasta::Base>& sequence: sequences) {
                rleSequences.emplace_back();
                rle(sequence, rleSequences.back());
            }

            // Check if they are all the same.
            bool allRleSequenceAreEqual = true;
            for(uint64_t i=1; i<sequences.size(); i++) {
                if(rleSequences[i] != rleSequences[0]) {
                    allRleSequenceAreEqual = false;
                    break;
                }
            }

            removeBubble = allRleSequenceAreEqual;
        }

        // For a diploid bubble, also remove it if the two sides
        // differ just by a copy number in a repeat with short period.
        if(sequences.size() == 2) {
            if(isCopyNumberDifference(sequences[0], sequences[1], maxRepeatPeriod)) {
                removeBubble = true;
            }
        }


        // Remove the bubble, if we decided that we can do that.
        if(removeBubble) {
            ++removedCount;

            // Remove the edges of the bubble.
            for(const edge_descriptor e: bubble.edges) {
                boost::remove_edge(e, assemblyGraph);
            }

            // Add a new edge with a single step to replace the bubble.
            edge_descriptor e;
            bool edgeWasAdded;
            tie(e, edgeWasAdded) = add_edge(bubble.v0, bubble.v1,
                AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            AssemblyGraphEdge& edge = assemblyGraph[e];

            const AnchorPair& anchorPair = p.second;
            const uint64_t offset = anchorPair.getAverageOffset(anchors);
            edge.emplace_back(anchorPair, offset);
        }
    }

    cout << "Out of " << allBubbles.size() << " bubbles, " <<
        candidateBubbles.size() << " were candidate for removal and " <<
        removedCount << " were actually removed." << endl;

    return removedCount;
}



// Compress linear chains of edges into a single edge.
void AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Find linear chains of 2 or more edges.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    for(const auto& chain: chains) {
        SHASTA_ASSERT(chain.size() > 1);

        // Get the first and last edge of this chain.
        const edge_descriptor e0 = chain.front();
        const edge_descriptor e1 = chain.back();

        // Get the first and last edge of this chain.
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

        // Concatenate the steps of all the edges in the chain.
        for(const edge_descriptor e: chain) {
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            copy(edge.begin(), edge.end(), back_inserter(edgeNew));
        }

        if(debug) {
            cout << "Compress";
            for(const edge_descriptor e: chain) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << " into " << edgeNew.id << endl;
        }

        // Now we can remove the edges of the chain and its internal vertices.
        bool isFirst = true;
        for(const edge_descriptor e: chain) {
            if(isFirst) {
                isFirst = false;
            } else {
                const vertex_descriptor v = source(e, assemblyGraph);
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }
        }

    }

}



void AssemblyGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AssemblyGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AssemblyGraph::save(const string& stage) const
{
    // If not using persistent binary data, do nothing.
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    // First save to a string.
    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    // cout << "Serialization of AssemblyGraph-" + stage << " needs " << dataString.size() << " bytes." << endl;

    // Now save the string to binary data.
    const string name = largeDataName("AssemblyGraph-" + stage);
    MemoryMapped::Vector<char> data;
    data.createNew(name, largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AssemblyGraph::load(const string& assemblyStage)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        const string name = largeDataName("AssemblyGraph-" + assemblyStage);
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error("Assembly graph at stage " + assemblyStage +
            " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error("Error reading assembly graph at stage " + assemblyStage +
            ": " + e.what());
    }
}



uint64_t AssemblyGraph::detangleVertices(uint64_t maxIterationCount, Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<maxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleVerticesIteration(detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            compress();
            cout << "Detangle vertices iteration " << iteration << ": " << changeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleVerticesIteration(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather vertices on which we will attempt detangling.
    // Each generates a tangle with just one vertex.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {

        if(
            (in_degree(v, assemblyGraph) > 1) and
            (out_degree(v, assemblyGraph) > 1)
             ) {
            detanglingCandidates.emplace_back(vector<vertex_descriptor>({v}));
        }

    }

    // Do the detangling.
    return detangle(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleEdges(
    uint64_t maxIterationCount,
    uint64_t maxEdgeLength,
    Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<maxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleEdgesIteration(maxEdgeLength, detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            compress();
            cout << "Detangle edges iteration " << iteration << ": " << iterationChangeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleEdgesIteration(
    uint64_t maxEdgeLength,
    Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    // cout << "Edge detangling begins." << endl;

    // Gather edges on which we will attempt detangling.
    // Each generates a tangle with just two vertices.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const bool isTangleEdge =
            (in_degree(v0, assemblyGraph) > 1) and
            (out_degree(v0, assemblyGraph) == 1) and
            (in_degree(v1, assemblyGraph) == 1) and
            (out_degree(v1, assemblyGraph) > 1);

        const bool isShort = assemblyGraph[e].offset() <= maxEdgeLength;

        if(isTangleEdge and isShort) {
            detanglingCandidates.emplace_back(vector<vertex_descriptor>({v0, v1}));
        }
    }

    // Do the detangling.
    return detangle(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleTemplates(uint64_t maxIterationCount, Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t templateId=0; templateId<tangleTemplates.size(); templateId++) {
        const uint64_t thisTemplateChangeCount = detangleTemplate(templateId, maxIterationCount, detangler);
        changeCount += thisTemplateChangeCount;
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleTemplate(
    uint64_t templateId,
    uint64_t maxIterationCount,
    Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<maxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleTemplateIteration(templateId, detangler);
        if(iterationChangeCount > 0) {
            compress();
            changeCount += iterationChangeCount;
            cout << "Detangle tangle template " << templateId << " iteration " << iteration << ": " << changeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleTemplateIteration(uint64_t templateId, Detangler& detangler)
{
    const TangleTemplate& tangleTemplate = tangleTemplates[templateId];

    AssemblyGraph& assemblyGraph = *this;
    vector< vector<vertex_descriptor> > detanglingCandidates;
    inducedSubgraphIsomorphisms(assemblyGraph, tangleTemplate, detanglingCandidates);
    // cout << "Found " << detanglingCandidates.size() << " instances of this tangle template." << endl;

    return detangle(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangle(
    const vector< vector<vertex_descriptor> >& detanglingCandidates,
    Detangler& detangler)
{
    const bool debug = false;
    // detangler.debug = true;

    AssemblyGraph& assemblyGraph = *this;

    std::set<vertex_descriptor> removedVertices;
    uint64_t attemptCount = 0;
    uint64_t successCount = 0;
    for(const vector<vertex_descriptor>& tangleVertices: detanglingCandidates) {

        // If any of the vertices in this tangle have been removed, by previous
        // detangling operations, skip it.
        bool skip = false;
        for(const vertex_descriptor v: tangleVertices) {
            if(removedVertices.contains(v)) {
                skip = true;
                break;
            }
        }
        if(skip) {
            continue;
        }



        // Attempt detangling for the tangle defined by these vertices.
        ++attemptCount;
        Tangle tangle(assemblyGraph, tangleVertices,
            options.aDrift,
            options.bDrift);
        if(debug) {
            const TangleMatrix& tangleMatrix = *(tangle.tangleMatrix);
            cout << "Tangle with " << tangleMatrix.entrances.size() << " entrances and " <<
                tangleMatrix.exits.size() << " exits." << endl;

            cout << "Entrances:";
            for(const auto& entrance: tangleMatrix.entrances) {
                cout << " " << assemblyGraph[entrance.e].id;
            }
            cout << endl;

            cout << "Exits:";
            for(const auto& exit: tangleMatrix.exits) {
                cout << " " << assemblyGraph[exit.e].id;
            }
            cout << endl;

            for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
                for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                    cout << tangleMatrix.tangleMatrix[iEntrance][iExit].orientedReadIds.size() << " ";
                }
            }
            cout << endl;
        }

        const bool success = detangler(tangle);
        if(success) {
            if(debug) {
                cout << "Detangle was successful." << endl;
            }
            for(vertex_descriptor v: tangle.removedVertices) {
                removedVertices.insert(v);
            }
            ++successCount;
        } else {
            if(debug) {
                cout << "Detangle failed." << endl;
            }
        }


    }
    // cout << "Attempted detangling for " << attemptCount << " tangles." << endl;
    // cout << "Detangling was successful for " << successCount << " tangles." << endl;



    return successCount;
}


#if 0
// One iteration of all usable detangling functions using the given Detangler.
bool AssemblyGraph::detangleIteration(Detangler& detangler)
{
    bool success = false;

    const bool verticesSuccess = detangleVertices(detangler);
    success = success or verticesSuccess;
    compress();

    const bool edgesSuccess = detangleEdges(detangler);
    success = success or edgesSuccess;
    compress();

    for(uint64_t tangleTemplateId=0; tangleTemplateId<tangleTemplates.size(); tangleTemplateId++) {
        const TangleTemplate& tangleTemplate = tangleTemplates[tangleTemplateId];
        /*
        cout << "Working on tangle template " << tangleTemplateId <<
            " with " << num_vertices(tangleTemplate) <<
            " vertices and " << num_edges(tangleTemplate) << " edges." << endl;
        */
        const bool templateSuccess = detangle(tangleTemplate, detangler);
        success = success or templateSuccess;
        compress();
    }

    return success;
}
#endif



uint64_t AssemblyGraph::detangle(
    uint64_t maxIterationCount,
    uint64_t maxEdgeLength,
    Detangler& detangler)
{

    const uint64_t verticesChangeCount = detangleVertices(maxIterationCount, detangler);
    const uint64_t edgesChangeCount = detangleEdges(maxIterationCount, maxEdgeLength, detangler);
    // const uint64_t templateChangeCount = detangleTemplates(maxIterationCount, detangler);

    const uint64_t changeCount = verticesChangeCount + edgesChangeCount; // + templateChangeCount;

    return changeCount;
}



void AssemblyGraph::prune()
{
    AssemblyGraph& assemblyGraph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    vector<edge_descriptor> edgesToBeRemoved;
    for(uint64_t iteration=0; iteration<options.pruneIterationCount; iteration++) {

        // Edge pruning for this iteration.
        edgesToBeRemoved.clear();
        uint64_t pruneCount = 0;
        uint64_t prunedLength = 0;
        BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
            const uint64_t offset = assemblyGraph[e].offset();

            // If long enough, don't prune it.
            if(offset > options.pruneLength) {
                continue;
            }

            // See if it can be pruned.
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if((in_degree(v0, assemblyGraph) == 0) or (out_degree(v1, assemblyGraph) == 0)) {
                edgesToBeRemoved.push_back(e);
                ++pruneCount;
                prunedLength += offset;
            }
        }

        // Remove the edges.
        for(const edge_descriptor e: edgesToBeRemoved) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Now remove any vertices that are left isolated.
        verticesToBeRemoved.clear();
        BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
            if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
                verticesToBeRemoved.push_back(v);
            }
        }
        for(const vertex_descriptor v: verticesToBeRemoved) {
            boost::remove_vertex(v, assemblyGraph);
        }

        cout << "Prune iteration removed " << pruneCount <<
            " edges with total estimated length " << prunedLength << endl;

        if(pruneCount == 0) {
            break;
        }
    }
}



// Tangle teplates are processed in the order defined here.
// They should entered in order of increasing complexity (number of edges).
void AssemblyGraph::createTangleTemplates()
{
    // Templates with 3 edges.

    tangleTemplates.emplace_back(3);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        tangleTemplates.push_back(reverse(g));
    }



    // Templates with 4 edges.

    tangleTemplates.emplace_back(4);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
    }

    tangleTemplates.emplace_back(4);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 3, g);
        tangleTemplates.push_back(reverse(g));
    }



    // Templates with 6 edges.

    tangleTemplates.emplace_back(5);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(1, 4, g);
        add_edge(2, 3, g);
        add_edge(2, 4, g);
        tangleTemplates.push_back(reverse(g));
    }



    // Templates with 7 edges.

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        tangleTemplates.push_back(reverse(g));
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        tangleTemplates.push_back(reverse(g));
    }

    tangleTemplates.emplace_back(6);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 3, g);
        add_edge(2, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
    }



    // Templates with 10 edges.

    tangleTemplates.emplace_back(8);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
        add_edge(5, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
    }

    tangleTemplates.emplace_back(8);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        add_edge(4, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
    }



    // Templates with 13 edges.

    tangleTemplates.emplace_back(10);
    {
        TangleTemplate& g = tangleTemplates.back();
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(1, 2, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(3, 4, g);
        add_edge(4, 5, g);
        add_edge(5, 6, g);
        add_edge(5, 6, g);
        add_edge(6, 7, g);
        add_edge(7, 8, g);
        add_edge(7, 8, g);
        add_edge(8, 9, g);
    }



    for(uint64_t i=0; i<tangleTemplates.size(); i++) {
        const TangleTemplate& tangleTemplate = tangleTemplates[i];
        const string dotFileName = "TangleTemplate-" + to_string(i) + ".dot";
        ofstream dot(dotFileName);
        writeGraphviz(dot, tangleTemplate);
        dot.close();
        std::system(("dot -O -T svg -Nshape=rectangle " + dotFileName).c_str());
    }
}




// This reverses all of the edges of the TangleTemplate.
// We can't use boost::reverse_graph because that creates
// a graph of a different type.
AssemblyGraph::TangleTemplate AssemblyGraph::reverse(const TangleTemplate& x)
{
    TangleTemplate y(num_vertices(x));

    BGL_FORALL_EDGES(e, x, TangleTemplate) {
        add_edge(target(e, x), source(e, x), y);
    }
    return y;
}



void AssemblyGraph::writeGraphviz(ostream& s, const TangleTemplate& g)
{
    s << "digraph TangleTemplate {\n";

    BGL_FORALL_VERTICES(v, g, TangleTemplate) {
        s << v << ";\n";
    }

    BGL_FORALL_EDGES(e, g, TangleTemplate) {
        const auto v0 = source(e, g);
        const auto v1 = target(e, g);
        s << v0 << "->" << v1 << ";\n";
    }

    s << "}\n";
}



// This finds all Superbubbles seen using the specified maxDistance.
// Some pairs of Superbubble can intersect (that is, they can have common edges).
void AssemblyGraph::findSuperbubbles(
    vector<Superbubble>& superbubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(assemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                forwardPairs.emplace_back(vA, vB);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    // cout << "Found " << forwardPairs.size() << " forward pairs." << endl;

    const boost::reverse_graph<AssemblyGraph> reverseAssemblyGraph(assemblyGraph);
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(reverseAssemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                backwardPairs.emplace_back(vB, vA);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(backwardPairs.begin(), backwardPairs.end());
    // cout << "Found " << backwardPairs.size() << " backward pairs." << endl;

    // Find pairs that appeared in both directions.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    // Each of these pairs generates a Superbubble.
    superbubbles.clear();
    for(const auto& p: bidirectionalPairs) {
        superbubbles.emplace_back(assemblyGraph, p.first, p.second);

    }

}



// Remove Superbubbles that are entirely contained in a larger superbubble.
void AssemblyGraph::removeContainedSuperbubbles(vector<Superbubble>& superbubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Find pairs of intersecting superbubbles.
    // Two superbubbles intersect if they have one or more internal edges in common.
    std::map<edge_descriptor, vector<uint64_t> > m;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles[superbubbleId];
        for(const edge_descriptor e: superbubble.internalEdges) {
            m[e].push_back(superbubbleId);
        }
    }
    std::set< pair<uint64_t, uint64_t> > intersectingPairs;
    for(const auto& p: m) {
        const vector<uint64_t>& edgeSuperbubbles = p.second;
        for(uint64_t i0=0; i0<edgeSuperbubbles.size()-1; i0++) {
            const uint64_t superbubbleId0 = edgeSuperbubbles[i0];
            for(uint64_t i1=i0+1; i1<edgeSuperbubbles.size(); i1++) {
                const uint64_t superbubbleId1 = edgeSuperbubbles[i1];
                intersectingPairs.insert({
                    min(superbubbleId0, superbubbleId1),
                    max(superbubbleId0, superbubbleId1)});
            }
        }
    }
    // cout << "Found " << intersectingPairs.size() << " intersecting superbubble pairs:" << endl;

    vector<edge_descriptor> commonEdges;
    vector<uint64_t> superbubblesToBeRemoved;
    for(const auto& p: intersectingPairs) {
        const uint64_t superbubbleId0 = p.first;
        const uint64_t superbubbleId1 = p.second;
        const Superbubble& superbubble0 = superbubbles[superbubbleId0];
        const Superbubble& superbubble1 = superbubbles[superbubbleId1];
        const auto& internalEdges0 = superbubble0.internalEdges;
        const auto& internalEdges1 = superbubble1.internalEdges;

        // Find the common edges.
        commonEdges.clear();
        std::set_intersection(
            internalEdges0.begin(), internalEdges0.end(),
            internalEdges1.begin(), internalEdges1.end(),
            back_inserter(commonEdges),
            AssemblyGraph::OrderById(assemblyGraph)
            );

        if(commonEdges.size() == internalEdges0.size()) {
            // cout << "Superbubble " << superbubbleId0 << " is contained in superbubble " << superbubbleId1 << endl;
            superbubblesToBeRemoved.push_back(superbubbleId0);
        } else if(commonEdges.size() == internalEdges1.size()) {
            // cout << "Superbubble " << superbubbleId1 << " is contained in superbubble " << superbubbleId0 << endl;
            superbubblesToBeRemoved.push_back(superbubbleId1);
        } else {
            cout << "Superbubbles " << superbubbleId0 << " and " << superbubbleId1 << " intersect." << endl;
            cout << "Superbubbles " << superbubbleId0 << " has " << internalEdges0.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges0) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Superbubbles " << superbubbleId1 << " has " << internalEdges1.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges1) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Found " << commonEdges.size() << " common edges." << endl;
            // This should never happen, but just in case we remove both of them.
            SHASTA_ASSERT(0);
            // superbubblesToBeRemoved.push_back(superbubbleId0);
            // superbubblesToBeRemoved.push_back(superbubbleId1);
        }
    }
    deduplicate(superbubblesToBeRemoved);
    // cout << "Removing " << superbubblesToBeRemoved.size() << " superbubbles that "
    //     "are contained in another superbubble." << endl;

    vector<Superbubble> newSuperbubbles;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(not binary_search(superbubblesToBeRemoved.begin(), superbubblesToBeRemoved.end(), superbubbleId)) {
            newSuperbubbles.push_back(superbubbles[superbubbleId]);
        }
    }
    newSuperbubbles.swap(superbubbles);
}



// This creates a csv file with one line of information for each superbubble.
void AssemblyGraph::writeSuperbubbles(
    const vector<Superbubble>& superbubbles,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Superbubble id,Type,Source ploidy,Target ploidy,Internal edges\n";

    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Superbubble& superbubble = superbubbles[id];

#if 0
        for(const auto& e: superbubble.sourceEdges) {
            cout << "Source edge " << assemblyGraph[e].id << endl;
        }
        for(const auto& e: superbubble.targetEdges) {
            cout << "Target edge " << assemblyGraph[e].id << endl;
        }
        for(const auto& e: superbubble.internalEdges) {
            cout << "Internal edge " << assemblyGraph[e].id << endl;
        }
#endif

        csv << id << ",";

        if(superbubble.isBubble()) {
            const uint64_t ploidy = superbubble.ploidy();
            if(ploidy == 1) {
                csv << "Edge,";
            } else {
                csv << "Bubble,";
            }
        } else {
            csv << "Superbubble,";
        }

        csv << superbubble.sourcePloidy() << ",";
        csv << superbubble.targetPloidy() << ",";

        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << ",";
        }
        csv << "\n";
    }
}



// This creates a csv file that can be loaded in Bandage to see the Superbubbles.
void AssemblyGraph::writeSuperbubblesForBandage(
    const vector<Superbubble>& superbubbles,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Segment,Color\n";

    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Superbubble& superbubble = superbubbles[id];
        const string color = randomHslColor(id, 0.75, 0.5);

        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << "," << color << "\n";
        }
    }

}



void AssemblyGraph::writeSuperbubbleChains(
    const vector<SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "ChainId,Position,Type,Source ploidy,Target ploidy,"
        "Internal vertices count,Internal edges count,Internal edges\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const  SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Superbubble& superbubble = superbubbleChain[position];
            csv << chainId << ",";
            csv << position << ",";

            if(superbubble.isBubble()) {
                const uint64_t ploidy = superbubble.ploidy();
                if(ploidy == 1) {
                    csv << "Edge,";
                } else {
                    csv << "Bubble,";
                }
            } else {
                csv << "Superbubble,";
            }

            csv << superbubble.sourcePloidy() << ",";
            csv << superbubble.targetPloidy() << ",";
            csv << superbubble.internalVertices.size() << ",";
            csv << superbubble.internalEdges.size() << ",";
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

}



void AssemblyGraph::writeSuperbubbleChainsForBandage(
    const vector<SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Segment,Color,Chain,Position\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        const string color = randomHslColor(chainId, 0.75, 0.5);

        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Superbubble& superbubble = superbubbleChain[position];
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv <<
                    assemblyGraph[e].id << "," <<
                    color << "," <<
                    chainId << "," <<
                    position << "\n";
            }
        }
    }
}



void AssemblyGraph::phaseSuperbubbleChains()
{

    // Find superbubbles.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    writeSuperbubbles(superbubbles, "Superbubbles-WithOverlaps.csv");
    removeContainedSuperbubbles(superbubbles);
    cout << "Found " << superbubbles.size() << " non-overlapping superbubbles." << endl;
    writeSuperbubbles(superbubbles, "Superbubbles.csv");
    writeSuperbubblesForBandage(superbubbles, "Superbubbles-Bandage.csv");

    // Find Superbubble chains.
    vector<SuperbubbleChain> superbubbleChains;
    findSuperbubbleChains(superbubbles, superbubbleChains);
    cout << "Found " << superbubbleChains.size() << " superbubble chains." << endl;
    writeSuperbubbleChains(superbubbleChains, "SuperbubbleChains.csv");
    writeSuperbubbleChainsForBandage(superbubbleChains, "SuperbubbleChains-Bandage.csv");

    // Phase them.
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
        superbubbleChain.phase(*this, superbubbleChainId);
    }

    compress();

}


// Simplify Superbubbles by turning them into bubbles via clustering
// of oriented read journeys.
void AssemblyGraph::simplifySuperbubbles()
{

    // Find the superbubbles, then remove superbubbles that are entirely
    // contained in another superbubble.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    removeContainedSuperbubbles(superbubbles);

    // Count the number of true Superbubbles, excluding bubbles.
    uint64_t count = 0;
    for(const Superbubble& superbubble: superbubbles) {
        if(not superbubble.isBubble()) {
            ++count;
        }
    }
    cout << "Found " << superbubbles.size() << " superbubbles of which " <<
        count << " are not simple bubbles." << endl;

    for(const Superbubble& superbubble: superbubbles) {
        if(not superbubble.isBubble()) {
            simplifySuperbubble(superbubble, options.simplifySuperbubbleMinCoverage);
        }
    }

}



void AssemblyGraph::simplifySuperbubble(
    const Superbubble& superbubble,
    uint64_t minCoverage)
{
    AssemblyGraph& assemblyGraph = *this;

    const AnchorId anchorIdA = assemblyGraph[superbubble.sourceVertex].anchorId;
    const AnchorId anchorIdB = assemblyGraph[superbubble.targetVertex].anchorId;

    const bool debug = false;

    if(debug) {
        cout << "Working on a superbubble consisting of the following " <<
            superbubble.internalEdges.size() << " edges:" << endl;
        for(const edge_descriptor e: superbubble.internalEdges) {
            cout << assemblyGraph[e].id << " ";
        }
        cout << endl;
    }

    // Create an AnchorPair between anchorIdA and anchorIdB
    // using all oriented reads in common between anchorIdA and anchorIdB.
    AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);

    if(debug) {
        cout << "The initial anchor pair has " << anchorPair.orientedReadIds.size() <<
            " oriented reads." << endl;
    }

    // We only want to use OrientedReadIds that appear at least once in the internal
    // edges of the superbubble.
    vector<OrientedReadId> allowedOrientedReadIds;
    for(const edge_descriptor e: superbubble.internalEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AssemblyGraphEdgeStep& step: edge) {
            const AnchorPair& stepAnchorPair = step.anchorPair;
            copy(stepAnchorPair.orientedReadIds.begin(), stepAnchorPair.orientedReadIds.end(),
                back_inserter(allowedOrientedReadIds));
        }
    }
    deduplicate(allowedOrientedReadIds);
    if(debug) {
        cout << "The internal edges of the superbubble contain " <<
            allowedOrientedReadIds.size() << " oriented reads." << endl;
    }

    // Only keep OrientedReadIds in the allowed set.
    vector<OrientedReadId> newOrientedReadIds;
    std::set_intersection(
        anchorPair.orientedReadIds.begin(), anchorPair.orientedReadIds.end(),
        allowedOrientedReadIds.begin(), allowedOrientedReadIds.end(),
        back_inserter(newOrientedReadIds));
    anchorPair.orientedReadIds.swap(newOrientedReadIds);

    if(debug) {
        cout << "The final anchor pair has " << anchorPair.orientedReadIds.size() <<
            " oriented reads." << endl;
    }

    // Cluster the oriented reads in the AnchorPair.
    vector<AnchorPair> newAnchorPairs;
    anchorPair.splitByClustering(anchors, journeys, options.clusteringMinJaccard, newAnchorPairs);



    if(debug) {
        cout << "Found " << newAnchorPairs.size() << " split AnchorPairs with sizes:";
        for(const AnchorPair& newAnchorPair: newAnchorPairs) {
            cout << " " << newAnchorPair.orientedReadIds.size();
        }
        cout << endl;

#if 0
        // Also write out assembled sequences.
        for(uint64_t i=0; i<newAnchorPairs.size(); i++) {
            const AnchorPair& newAnchorPair = newAnchorPairs[i];
            ostream html(0);
            LocalAssembly2 localAssembly(
                anchors, html, false,
                options.aDrift,
                options.bDrift,
                newAnchorPair);
            localAssembly.run(false, options.localAssemblyOptions.maxAbpoaLength);
            vector<shasta::Base> sequence;
            localAssembly.getSequence(sequence);
            cout << ">" << i << endl;
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(cout));
            cout << endl;
        }
#endif
    }

    // Only keep the ones with coverage at least minCoverage.
    for(uint64_t i=0; i<newAnchorPairs.size(); i++) {
        if(newAnchorPairs[i].orientedReadIds.size() < minCoverage) {
            newAnchorPairs.resize(i);
            break;
        }
    }

    if(debug) {
        cout << "Kept " << newAnchorPairs.size() << " split AnchorPairs with sizes:";
        for(const AnchorPair& newAnchorPair: newAnchorPairs) {
            cout << " " << newAnchorPair.orientedReadIds.size();
        }
        cout << endl;
    }

    // If there are no AnchorPairs with sufficient coverage, we can't simplify this Superbubble.
    if(newAnchorPairs.empty()) {
        if(debug) {
            cout << "This superbubble cannot be simplified because there are no usable anchor pairs." << endl;
        }
        return;
    }

    // We replace the Superbubble with a bubble created using these new AnchorPairs.
    // Each of the new AnchorPairs we kept will generate a branch of the new bubble.
    for(const AnchorPair& newAnchorPair: newAnchorPairs) {
        const uint64_t offset = newAnchorPair.getAverageOffset(anchors);
        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(
            superbubble.sourceVertex, superbubble.targetVertex,
            AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.emplace_back(newAnchorPair, offset);
        if(debug) {
            cout << "Created new edge " << edge.id << endl;
        }
    }

    // Now we can remove all the internal edges of the Superbubble.
    for(const edge_descriptor e: superbubble.internalEdges) {
        if(debug) {
            cout << "Removed edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // Also remove the internal vertices of the Superbubble.
    for(const vertex_descriptor v: superbubble.internalVertices) {
        SHASTA_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA_ASSERT(out_degree(v, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
    }
}



void AssemblyGraph::findSuperbubbleChains(
    const vector<Superbubble>& superbubbles,
    vector<SuperbubbleChain>& superbubbleChains
    ) const
{
    // Index the superbubbles by their source and target vertex.
    // FOR REPRODUCIBILITY WE SHOULD CHANGE THIS TO INDEX BY VERTEX ID INSTEAD.
    std::map<vertex_descriptor, vector<uint64_t> > mapBySource;
    std::map<vertex_descriptor, vector<uint64_t> > mapByTarget;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles[superbubbleId];
        mapBySource[superbubble.sourceVertex].push_back(superbubbleId);
        mapByTarget[superbubble.targetVertex].push_back(superbubbleId);
    }

    // Sanity check: a vertex can only be a source or target of a single Superbubble.
    for(const auto& p: mapBySource) {
        SHASTA_ASSERT(p.second.size() == 1);
    }
    for(const auto& p: mapByTarget) {
        SHASTA_ASSERT(p.second.size() == 1);
    }

    // A vector to keep track of the Superbubbles we already added to a chain.
    vector<bool> wasUsed(superbubbles.size(), false);

    // Work vectors used below to construct chains.
    vector<uint64_t> forward;
    vector<uint64_t> backward;

    // Generate Superbubble chains.
    superbubbleChains.clear();
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(wasUsed[superbubbleId]) {
            continue;
        }
        // cout << "Starting a new chain at Superbubble " << superbubbleId << endl;

        // This Superbubble has not yet been added to any chain.
        // We will start a new chain here.

        // Create the forward portion of this chain.
        forward.clear();
        vertex_descriptor v = superbubbles[superbubbleId].targetVertex;
        while(true) {
            const auto it = mapBySource.find(v);
            if(it == mapBySource.end()) {
                break;
            }
            const vector<uint64_t>& nextVector = it->second;
            SHASTA_ASSERT(nextVector.size() == 1);
            const uint64_t nextSuperbubbleId = nextVector.front();
            forward.push_back(nextSuperbubbleId);
            // cout << "Forward: " << nextSuperbubbleId << endl;

            v = superbubbles[nextSuperbubbleId].targetVertex;
        }

        // Create the backward portion of this chain.
        backward.clear();
        v = superbubbles[superbubbleId].sourceVertex;
        while(true) {
            const auto it = mapByTarget.find(v);
            if(it == mapByTarget.end()) {
                break;
            }
            const vector<uint64_t>& previousVector = it->second;
            SHASTA_ASSERT(previousVector.size() == 1);
            const uint64_t previousSuperbubbleId = previousVector.front();
            backward.push_back(previousSuperbubbleId);
            // cout << "Backward: " << previousSuperbubbleId << endl;

            v = superbubbles[previousSuperbubbleId].sourceVertex;
        }

        // Now we can create the new SuperbubbleChain.
        superbubbleChains.emplace_back();
        SuperbubbleChain& superbubbleChain = superbubbleChains.back();
        std::reverse(backward.begin(), backward.end());
        for(const uint64_t id: backward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
        superbubbleChain.push_back(superbubbles[superbubbleId]);
        wasUsed[superbubbleId] = true;
        for(const uint64_t id: forward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
    }
}



// Find the non-trivial strongly connected components.
// Each component is stored with vertices sorted to permit binary searches.
void AssemblyGraph::findStrongComponents(
    vector< vector<vertex_descriptor> >& strongComponents) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        assemblyGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    strongComponents.clear();
    for(const auto& p: componentVertices) {
        const vector<vertex_descriptor>& component = p.second;
        if(component.size() > 1) {
            strongComponents.push_back(component);
            sort(strongComponents.back().begin(), strongComponents.back().end());
        }
    }


}



// This creates a csv file that can be loaded in bandage to see
// the strongly connected components.
void AssemblyGraph::colorStrongComponents() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< vector<vertex_descriptor> > strongComponents;
    findStrongComponents(strongComponents);

    ofstream csv("StrongComponents.csv");
    csv << "Id,Color,Component\n";

    // Loop over the non-trivial strongly connected components.
    for(uint64_t i=0; i<strongComponents.size(); i++) {
        const vector<vertex_descriptor>& strongComponent = strongComponents[i];

        for(const vertex_descriptor v0: strongComponent) {
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                if(std::binary_search(strongComponent.begin(), strongComponent.end(), v1)) {
                    csv << assemblyGraph[e].id << ",Green," << i << "\n";
                }
            }
        }

    }


}




double AssemblyGraphEdge::averageCoverage() const
{
    uint64_t sum = 0;
    for(const AssemblyGraphEdgeStep& step: *this) {
        sum += step.anchorPair.orientedReadIds.size();
    }

    return double(sum) / double(size());
}



// Compute oriented read journeys in the AssemblyGraph.
void AssemblyGraph::computeJourneys() const
{
    const uint64_t orientedReadCount = journeys.size();

    const AssemblyGraph& assemblyGraph = *this;

    class AssemblyGraphJourneyEntry {
    public:
        edge_descriptor e;
        uint64_t stepId;

        // These are positions in the standard Journeys on the Anchors.
        uint64_t positionInJourneyA;
        uint64_t positionInJourneyB;

        bool operator<(const AssemblyGraphJourneyEntry& that) const {
            return positionInJourneyA < that.positionInJourneyA;
        }
    };
    vector< vector<AssemblyGraphJourneyEntry> > assemblyGraphJourneys(orientedReadCount);


    // Loop over AssemblyGrapg edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Loop over steps of this edge.
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];

            // Locate the anchors for this step.
            const AnchorPair& anchorPair = step.anchorPair;
            const AnchorId anchorIdA = anchorPair.anchorIdA;
            const AnchorId anchorIdB = anchorPair.anchorIdB;
            const Anchor anchorA = anchors[anchorIdA];
            const Anchor anchorB = anchors[anchorIdB];

            // Loop over OrientedReadIds of this step.
            auto itA = anchorA.begin();
            auto itB = anchorB.begin();
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {

                // Locate this OrientedReadId in the two anchors.
                for(; (itA != anchorA.end()) and (itA->orientedReadId != orientedReadId); ++itA) {}
                SHASTA_ASSERT(itA != anchorA.end());
                SHASTA_ASSERT(itA->orientedReadId == orientedReadId);
                const AnchorMarkerInfo& infoA = *itA;
                for(; (itB != anchorB.end()) and (itB->orientedReadId != orientedReadId); ++itB) {}
                SHASTA_ASSERT(itB != anchorB.end());
                const AnchorMarkerInfo& infoB = *itB;
                SHASTA_ASSERT(itB->orientedReadId == orientedReadId);

                const uint32_t positionInJourneyA = infoA.positionInJourney;
                const uint32_t positionInJourneyB = infoB.positionInJourney;

                AssemblyGraphJourneyEntry entry;
                entry.e = e;
                entry.stepId = stepId;
                entry.positionInJourneyA = positionInJourneyA;
                entry.positionInJourneyB = positionInJourneyB;

                assemblyGraphJourneys[orientedReadId.getValue()].push_back(entry);

            }

        }
    }

    for(vector<AssemblyGraphJourneyEntry>& v: assemblyGraphJourneys) {
        sort(v.begin(), v.end());
    }

    ofstream csv("AssemblyGraphJourneys.csv");
    csv << "OrientedReadId,Segment,Step,PositionInJourneyA,PositionInJourneyB\n";
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
        const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
        for(const AssemblyGraphJourneyEntry& entry: assemblyGraphJourney) {
            csv << orientedReadId << ",";
            csv << assemblyGraph[entry.e].id << ",";
            csv << entry.stepId << ",";
            csv << entry.positionInJourneyA << ",";
            csv << entry.positionInJourneyB << "\n";
        }
    }
}



// Search starting at a given edge (segment) and moving in the specified direction.
void AssemblyGraph::search(
    edge_descriptor eStart,
    uint64_t direction) const
{
    const AssemblyGraph& assemblyGraph = *this;

    const bool debug = true;

    // Create the SearchGraph
    SearchGraph searchGraph;
    const SearchGraph::vertex_descriptor svStart = add_vertex(SearchGraphVertex(eStart), searchGraph);
    std::map<AssemblyGraph::edge_descriptor, SearchGraph::vertex_descriptor> vertexMap;
    vertexMap.insert({eStart, svStart});

    // Work vectors.
    vector<AssemblyGraph::edge_descriptor> eStartVector(1, eStart);
    vector<AssemblyGraph::edge_descriptor> adjacentEdges;

    // Initialize a BFS.
    std::queue<SearchGraph::vertex_descriptor> q;
    q.push(svStart);

    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a SearchGraph vertex, which corresponds to an AssemblyGraph edge.
        const SearchGraph::vertex_descriptor sv0 = q.front();
        q.pop();
        const AssemblyGraph::edge_descriptor e0 = searchGraph[sv0].e;
        if(debug) {
            cout << "Dequeued " << assemblyGraph[e0].id << endl;
        }

        // Find adjacent edges in the specified direction.
        const vertex_descriptor v0 = (direction == 0) ? target(e0, assemblyGraph) : source(e0, assemblyGraph);
        adjacentEdges.clear();
        if(direction == 0) {
            BGL_FORALL_OUTEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
                adjacentEdges.push_back(e1);
            }
        } else {
            BGL_FORALL_INEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
                adjacentEdges.push_back(e1);
            }
        }

        // Compute a TangleMatrix between the start edge and these edges.
        TangleMatrix tangleMatrix(
            assemblyGraph,
            eStartVector,
            adjacentEdges,
            options.aDrift,
            options.bDrift);


        if(debug) {
            cout << "Tangle matrix:" << endl;
            for(uint64_t i=0; i<adjacentEdges.size(); i++) {
                const AssemblyGraph::edge_descriptor e1 = adjacentEdges[i];
                const uint64_t tangleMatrixCoverage = tangleMatrix.tangleMatrix[0][i].orientedReadIds.size();
                cout << assemblyGraph[e1].id << " " << tangleMatrixCoverage << endl;
            }
        }


        // Enqueue the segments with non-zero tangle matrix.
        for(uint64_t i=0; i<adjacentEdges.size(); i++) {
            const uint64_t tangleMatrixCoverage = tangleMatrix.tangleMatrix[0][i].orientedReadIds.size();
            if(tangleMatrixCoverage == 0) {
                continue;
            }
            const AssemblyGraph::edge_descriptor e1 = adjacentEdges[i];

            // Get the SearchGraph vertex corresponding to this AssemblyGraph edge, creating it
            // if necessary.
            SearchGraph::vertex_descriptor sv1;
            auto it = vertexMap.find(e1);
            if(it == vertexMap.end()) {
                sv1 = add_vertex(SearchGraphVertex(e1), searchGraph);
                vertexMap.insert({e1, sv1});
                q.push(sv1);
            } else {
                sv1 = it->second;
            }
            add_edge(sv0, sv1, SearchGraphEdge(tangleMatrixCoverage), searchGraph);
        }
    }
}



// More systematic search functionality that uses indexes.
void AssemblyGraph::findEdgePairs(uint64_t minCoverage)
{
    const AssemblyGraph& assemblyGraph = *this;

    const ReadId readCount = anchors.reads.readCount();
    const ReadId orientedReadCount = 2 * readCount;
    const OrderById orderById(*this);

    // Create an index that will be used for searches.
    // index[orientedReadId.getValue()] will contain the edges that
    // have orientedReadId in the AnchorPair of their first step.
    // Edge descriptors in the index are sorted using OrderById.
    vector< vector<edge_descriptor> > searchIndex(orientedReadCount);

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const AssemblyGraphEdgeStep& firstStep = edge.front();
        const AnchorPair& anchorPair = firstStep.anchorPair;
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            searchIndex[orientedReadId.getValue()].push_back(e);
        }
    }
    for(auto& v: searchIndex) {
        sort(v.begin(), v.end(), orderById);
    }



    // Now do a search from the last step of each edge.
    vector<edge_descriptor> currentEdges;
    vector<uint64_t> count;
    edgePairsBySource.clear();
    edgePairsByTarget.clear();
    BGL_FORALL_EDGES(e0, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
        const AssemblyGraphEdgeStep& lastStep = edge0.back();

        currentEdges.clear();
        const AnchorPair& anchorPair = lastStep.anchorPair;
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            const vector<edge_descriptor>& v = searchIndex[orientedReadId.getValue()];
            for(const edge_descriptor e1: v) {
                currentEdges.push_back(e1);
            }
        }
        deduplicateAndCountWithThreshold(currentEdges, count, minCoverage);

        for(const edge_descriptor e1: currentEdges) {
            edgePairsBySource[e0].push_back(e1);
            edgePairsByTarget[e1].push_back(e0);
        }
    }

    // Sort the pairs we found for each edge.
    for(auto& p:edgePairsBySource) {
        auto& v = p.second;
        sort(v.begin(), v.end(), orderById);
    }
    for(auto& p:edgePairsByTarget) {
        auto& v = p.second;
        sort(v.begin(), v.end(), orderById);
    }

}



void AssemblyGraph::testSearch(uint64_t edgeId0, uint64_t direction, uint64_t minCoverage) const
{
    const AssemblyGraph& assemblyGraph = *this;
    OrderById orderById(*this);


    edge_descriptor e0;
    bool found = false;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].id == edgeId0) {
            e0 = e;
            found = true;
            break;
        }
    }
    if(not found) {
        cout << "Edge with id " << edgeId0 << " does not exist." << endl;
        return;
    }

    std::map<edge_descriptor, vector<edge_descriptor> >::const_iterator it;
    if(direction == 0) {
        it = edgePairsBySource.find(e0);
        SHASTA_ASSERT(it != edgePairsBySource.end());
    } else {
        it = edgePairsByTarget.find(e0);
        SHASTA_ASSERT(it != edgePairsByTarget.end());
    }
    const vector<edge_descriptor> e1s = it->second;
    cout << "Found " << e1s.size() << " edges." << endl;




    // Compute the AnchorPairs between e0 and each of these e1.
    vector<edge_descriptor> e1sGood;
    for(const edge_descriptor e1: e1s) {
        edge_descriptor eA = e0;
        edge_descriptor eB = e1;
        if(direction == 1) {
            std::swap(eA, eB);
        }

        const AnchorPair bridgeAnchorPair = anchors.bridge(
            assemblyGraph[eA].back().anchorPair,
            assemblyGraph[eB].front().anchorPair,
            options.aDrift,
            options.bDrift);
        if(bridgeAnchorPair.orientedReadIds.size() >= minCoverage) {
            e1sGood.push_back(e1);
        }
    }
    cout << "Found " << e1sGood.size() << " good edges." << endl;



    // Write a csv file that can be loaded in Bandage.
    ofstream csv("Search.csv");
    csv << "Id,Color\n";
    BGL_FORALL_EDGES(e1, assemblyGraph, AssemblyGraph) {
        csv << assemblyGraph[e1].id << ",";
        if(e1 == e0) {
            csv << "Red";
        } else {
            const auto it = std::lower_bound(e1sGood.begin(), e1sGood.end(), e1, orderById);
            if((it != e1s.end()) and (*it == e1)) {
                csv << "Green";
            } else {
                csv << "Grey";
            }
        }
        csv << endl;
    }

}



// Local search that continues as long as we have only one way to move.
void AssemblyGraph::forwardLocalSearch(
    edge_descriptor eStart,
    uint64_t minCoverage,
    vector<edge_descriptor>& edges) const
{
    const bool debug = false;
    const AssemblyGraph& assemblyGraph = *this;

    edges.clear();
    edge_descriptor e0 = eStart;

    vector<edge_descriptor> entrances(1, eStart);
    vector<edge_descriptor> exits;

    // Main iteration.
    while(true) {
        const vertex_descriptor v0 = target(e0, assemblyGraph);
        exits.clear();
        BGL_FORALL_OUTEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
            exits.push_back(e1);
        }

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                cout << assemblyGraph[exits[iExit]].id << " " <<
                    tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() << endl;
            }
        }

        // Count the exits with coverage at least minCoverage.
        uint64_t goodExitCount = 0;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() >= minCoverage) {
                ++goodExitCount;
            }
        }

        // Only keep going if we have exactly one good exit.
        if(goodExitCount != 1) {
            return;
        }

        // Find the one and only good exit.
        edge_descriptor e1;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() >= minCoverage) {
                e1 = exits[iExit];
                break;
            }
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::backwardLocalSearch(
    edge_descriptor eStart,
    uint64_t minCoverage,
    vector<edge_descriptor>& edges) const
{
    const bool debug = false;
    const AssemblyGraph& assemblyGraph = *this;

    edges.clear();
    edge_descriptor e0 = eStart;

    vector<edge_descriptor> entrances;
    vector<edge_descriptor> exits(1, eStart);

    // Main iteration.
    while(true) {
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        entrances.clear();
        BGL_FORALL_INEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
            entrances.push_back(e1);
        }

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                cout << assemblyGraph[entrances[iEntrance]].id << " " <<
                    tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() << endl;
            }
        }

        // Count the enrances with coverage at least minCoverage.
        uint64_t goodEntranceCount = 0;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() >= minCoverage) {
                ++goodEntranceCount;
            }
        }

        // Only keep going if we have exactly one good entrance.
        if(goodEntranceCount != 1) {
            return;
        }

        // Find the one and only good entrance.
        edge_descriptor e1;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() >= minCoverage) {
                e1 = entrances[iEntrance];
                break;
            }
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::testLocalSearch(
    uint64_t edgeIdStart,
    uint64_t direction,
    uint64_t minCoverage) const
{
    const AssemblyGraph& assemblyGraph = *this;
    OrderById orderById(*this);

    edge_descriptor eStart;
    bool found = false;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].id == edgeIdStart) {
            eStart = e;
            found = true;
            break;
        }
    }
    if(not found) {
        cout << "Edge with id " << edgeIdStart << " does not exist." << endl;
        return;
    }

    vector<edge_descriptor> edges;
    if(direction == 0) {
        forwardLocalSearch(eStart, minCoverage, edges);
    } else {
        backwardLocalSearch(eStart, minCoverage, edges);
    }

    cout << "Found " << edges.size() << " edges:" << endl;
    for(const edge_descriptor e: edges) {
        cout << assemblyGraph[e].id << " ";
    }
    cout << endl;


}



void AssemblyGraph::createSearchGraph(uint64_t minCoverage) const
{
    shasta::SearchGraph searchGraph(*this, minCoverage);
    searchGraph.writeGraphviz("SearchGraph.dot");
}
