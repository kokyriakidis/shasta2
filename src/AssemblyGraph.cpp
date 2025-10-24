// Must include this first due to some Boost include file problem.
#include <boost/pending/disjoint_sets.hpp>

// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "areSimilarSequences.hpp"
#include "color.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "findLinearChains.hpp"
#include "findConvergingVertex.hpp"
#include "inducedSubgraphIsomorphisms.hpp"
#include "Journeys.hpp"
#include "LikelihoodRatioDetangler.hpp"
#include "LocalAssembly3.hpp"
#include "MurmurHash2.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "rle.hpp"
#include "SearchGraph.hpp"
#include "Superbubble.hpp"
#include "SuperbubbleChain.hpp"
#include "Tangle.hpp"
#include "Tangle1.hpp"
#include "TangleMatrix.hpp"
#include "TangleMatrix1.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/container/small_vector.hpp>
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
    options(options),
    orderById(*this)
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
    options(options),
    orderById(*this)
{
    load(stage);
}



// Detangle, phase, assemble sequence, output.
void AssemblyGraph::simplifyAndAssemble()
{
    performanceLog << timestamp << "Assembly graph::simplifyAndAssemble begins." << endl;

    LikelihoodRatioDetangler detangler(
        options.detangleEpsilon,
        options.detangleMaxLogP,
        options.detangleMinLogPDelta,
        options.detangleMinCoverage,
        false,  // Consider all hypotheses
        false,  // Consider all hypotheses
        true,   // Require top hypothesis to be injective
        true);  // Require top hypothesis to be a permutation.

    // Initial output.
    write("A");

    for(uint64_t iteration=0; iteration<options.simplifyMaxIterationCount; iteration++) {

        uint64_t changeCount = 0;

        // Prune.
        changeCount += prune();
        changeCount += compress();

        // Simplify Superbubbles.
        changeCount += simplifySuperbubbles(detangler);
        write("B" + to_string(iteration));

        // Remove or simplify bubbles likely caused by errors.
        changeCount += bubbleCleanup();
        changeCount += compress();
        write("C" + to_string(iteration));

        // Phase SuperbubbleChains, considering all hypotheses.
        changeCount += phaseSuperbubbleChains(true, true);
        write("D" + to_string(iteration));

        // Detangling.
        changeCount += detangleHighLevel(detangler);
        write("E" + to_string(iteration));

        cout << "Total change count at iteration " << iteration << " was " << changeCount << endl;
    }


    // Sequence assembly.
    assembleAll();
    write("Z");
    writeFasta("Z");

    performanceLog << timestamp << "Assembly graph::simplifyAndAssemble ends." << endl;
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
    cout << "Stage " << stage << ": " <<
        num_vertices(*this) << " vertices, " <<
        num_edges(*this) << " edges. Next edge id is " << nextEdgeId << "." << endl;

    save(stage);
    writeGfa("AssemblyGraph-" + stage + ".gfa");
    writeGraphviz("AssemblyGraph-" + stage + ".dot");
    writeCsv("AssemblyGraph-" + stage + ".csv");
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
        const double coverage = edge.averageCoverage();

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(gfa));
            const uint64_t length = sequence.size();
            gfa << "\tLN:i:" << length;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(length)));
            gfa << "\n";

        } else {
            const uint64_t offset = edge.offset();
            gfa << "*\tLN:i:" << offset;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
            gfa << "\n";
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
void AssemblyGraph::assembleAll()
{
    performanceLog << timestamp << "Sequence assembly begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;

    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble();
    edgesToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



// Assemble sequence for the specified edge.
void AssemblyGraph::assemble(edge_descriptor e)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble();
}



// Assemble sequence for step i of the specified edge.
// This is the lowest level sequence assembly function and is not multithreaded.
// It runs a LocalAssembly2 on the AnchorPair for that step.
void AssemblyGraph::assembleStep(edge_descriptor e, uint64_t i)
{
    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphEdge& edge = assemblyGraph[e];
    AssemblyGraphEdgeStep& step = edge[i];

    // cout << "Begin local assembly for edge " << edge.id << " " << " step " << i << endl;

    if(step.anchorPair.anchorIdA == step.anchorPair.anchorIdB) {
        step.sequence.clear();
        return;
    }

#if 0
    // Run the LocalAssembly.
    ofstream html;  // Not open, so no html output takes place.
    LocalAssembly localAssembly(
        anchors, html, false,
        options.aDrift,
        options.bDrift,
        step.anchorPair);
    if(localAssembly.coverage() == 0) {
        throw runtime_error("No coverage for local assembly at assembly graph edge " +
            to_string(edge.id) + " step " + to_string(i));
    }
    localAssembly.run(false, options.maxAbpoaLength);
    localAssembly.getSequence(step.sequence);
#endif


    // Let LocalAssembly3 use OrientedReadIds from the previous and next step.
    vector<OrientedReadId> additionalOrientedReadIds;
    if(i > 0) {
        std::ranges::copy(edge[i - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    if(i < edge.size() - 1) {
        std::ranges::copy(edge[i + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    deduplicate(additionalOrientedReadIds);

    ostream html(0);
    LocalAssembly3 localAssembly(
        anchors,
        html,
        false,
        edge[i].anchorPair,
        additionalOrientedReadIds);
    step.sequence = localAssembly.sequence;
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the stepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph::assemble()
{
    performanceLog << timestamp << "Sequence assembly begins for " << edgesToBeAssembled.size() <<
        " assembly graph edges." << endl;
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
    runThreads(&AssemblyGraph::assembleThreadFunction, options.threadCount);

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }

    edgesToBeAssembled.clear();
    stepsToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



void AssemblyGraph::assembleThreadFunction(uint64_t /* threadId */)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {

            /*
            if((j % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << j << "/" << stepsToBeAssembled.size() << endl;
            }
            */

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
                std::ranges::sort(bubble.edges, orderById);
                bubbles.push_back(bubble);
            }
        }
    }

#if 0
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
#endif
}



uint64_t AssemblyGraph::bubbleCleanup()
{
    uint64_t modifiedCount = 0;

    vector< pair<vertex_descriptor, vertex_descriptor> > excludeList;
    while(true) {
        const uint64_t modifiedCountThisIteration = bubbleCleanupIteration(excludeList);
        if(modifiedCountThisIteration == 0) {
            break;
        }
        modifiedCount += modifiedCountThisIteration;
    }

    return modifiedCount;
}



uint64_t AssemblyGraph::bubbleCleanupIteration(
    vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList)
{
    performanceLog << timestamp << "Bubble cleanup iteration begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find all bubbles.
    vector<Bubble> allBubbles;
    findBubbles(allBubbles);
    // cout << "Found " << allBubbles.size() << " bubbles." << endl;

    // Find candidate bubbles.
    // These are the ones that are not in the exclude list and
    // in which no branch has offset greater than
    // options.bubbleCleanupMaxBubbleLength.
    vector<Bubble> candidateBubbles;
    for(const Bubble& bubble: allBubbles) {

        if(std::ranges::contains(excludeList, make_pair(bubble.v0, bubble.v1))) {
            continue;
        }

        bool hasLongBranch = false;
        for(const edge_descriptor e: bubble.edges) {
            if(assemblyGraph[e].offset() > options.bubbleCleanupMaxBubbleLength) {
                hasLongBranch = true;
                break;
            }
        }

        if(not hasLongBranch) {
            candidateBubbles.push_back(bubble);
        }
    }
    // cout << candidateBubbles.size() << " bubbles are candidate for clean up." << endl;

    // Assemble sequence for all the edges of these bubbles.
    edgesToBeAssembled.clear();
    for(const Bubble& bubble: candidateBubbles) {
        for(const edge_descriptor e: bubble.edges) {
            if(not assemblyGraph[e].wasAssembled) {
                edgesToBeAssembled.push_back(e);
            }
        }
    }
    assemble();


    uint64_t modifiedCount = 0;
    for(const Bubble& bubble: candidateBubbles) {
        if(bubbleCleanup(bubble)) {
            ++modifiedCount;
        }
    }
    // cout << "Bubble cleanup modified " << modifiedCount << " bubbles." << endl;

    // Update the excludeList.
    for(const Bubble& bubble: candidateBubbles) {
        excludeList.push_back(make_pair(bubble.v0, bubble.v1));
    }
    std::ranges::sort(excludeList);

    performanceLog << timestamp << "Bubble cleanup iteration ends." << endl;
    return modifiedCount;
}



bool AssemblyGraph::bubbleCleanup(const Bubble& bubble)
{
    // EXPOSE WHEN CODE STABILIZES.
    const vector<uint64_t> minRepeatCount = {0, 2, 2, 2, 2, 2, 2};

    AssemblyGraph& assemblyGraph = *this;
    const uint64_t ploidy = bubble.edges.size();

    const bool debug = false; // (assemblyGraph[bubble.edges.front()].id)== 130581;
    if(debug) {
        cout << "Attempting cleanup for bubble";
        for(uint64_t i=0; i<ploidy; i++) {
            const edge_descriptor e = bubble.edges[i];
            cout << " (" << i << " " << assemblyGraph[e].id << ")";
        }
        cout << endl;
    }

    // Find which pairs of branches in the bubble have similar sequence
    // as defined by the minRepeatCount vector.
    // See analyzeBubble for more information.
    vector< pair<uint64_t, uint64_t> > similarPairs;
    analyzeBubble(bubble, minRepeatCount, similarPairs);

    // If no similar pairs were found, leave this bubble alone.
    if(similarPairs.empty()) {
        if(debug) {
            cout << "No similar branch sequences were found. Nothing done for this bubble." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "These pairs of branches have similar sequence:" << endl;
        for(const auto& p: similarPairs) {
            const uint64_t i0 = p.first;
            const uint64_t i1 = p.second;
            cout << i0 << " " << i1 << " " << assemblyGraph[bubble.edges[i0]].id <<
                " " << assemblyGraph[bubble.edges[i1]].id << endl;
        }
    }

    // Find the OrientedReadIds that appear in all the steps of each branch of the bubble.
    vector< vector<OrientedReadId> > allOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        const AssemblyGraphEdge& edge = assemblyGraph[bubble.edges[i]];
        for(const auto& step: edge) {
            std::ranges::copy(step.anchorPair.orientedReadIds, back_inserter(allOrientedReadIds[i]));
        }
        deduplicate(allOrientedReadIds[i]);
        if(debug) {
            cout << "Branch " << i << " has " << allOrientedReadIds[i].size() <<
                " total oriented read ids." << endl;
        }
    }

#if 0
    // Gather all OrientedReadIds that appear in exactly one of the branches.
    vector<OrientedReadId> unambiguousOrientedReadIdsAllBranches;
    for(const auto& v: allOrientedReadIds) {
        std::ranges::copy(v, back_inserter(unambiguousOrientedReadIdsAllBranches));
    }
    deduplicateAndCountAndKeepUnique(unambiguousOrientedReadIdsAllBranches);

    // The unambiguous OrientedReadIds for each branch are the intersection
    // of unambiguousOrientedReadIdsAllBranches with allOrientedReadIds for that branch.
    vector< vector<OrientedReadId> > unambiguousOrientedReadIds(bubble.edges.size());
    for(uint64_t i=0; i<ploidy; i++) {
        std::ranges::set_intersection(
            unambiguousOrientedReadIdsAllBranches,
            allOrientedReadIds[i],
            back_inserter(unambiguousOrientedReadIds[i]));
        if(debug) {
            cout << "Branch " << i << " has " << unambiguousOrientedReadIds[i].size() <<
                " unambiguous oriented read ids." << endl;
        }
    }
#endif

    // Create an AnchorPair between the source and target of this bubble.
    const AnchorId anchorId0 = assemblyGraph[bubble.v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[bubble.v1].anchorId;
    const AnchorPair bubbleAnchorPair(anchors, anchorId0, anchorId1, false);

    // The usable OrientedReadIds for each branch are the intersection of
    // allOrientedReads for the branch with the OrientedReadIds in the bubbleAnchorPair.
    vector< vector<OrientedReadId> > usableOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        std::ranges::set_intersection(
            bubbleAnchorPair.orientedReadIds,
            allOrientedReadIds[i],
            back_inserter(usableOrientedReadIds[i]));
        if(debug) {
            cout << "Branch " << i << " has " << usableOrientedReadIds[i].size() <<
                " usable oriented read ids." << endl;
        }
    }



    // Compute groups of similar branches. Each group will generate a new branch.
    // This uses a disjoint sets data structure created from the similarPairs vector,
    // except for some common special cases.
    vector< vector<uint64_t> > branchGroups;
    if(ploidy == 2) {

        SHASTA_ASSERT(similarPairs.size() == 1);    // We already checked that is is not empty.
        // Create a single branch group that includes both branches.
        branchGroups.push_back({0, 1});

    } else if(similarPairs.size() == (ploidy * (ploidy - 1)) / 2) {

        // Create a single branch group that includes all branches.
        branchGroups.resize(1);
        for(uint64_t i=0; i<ploidy; i++) {
            branchGroups.front().push_back(i);
        }

    } else {

        // Initialize the disjoint sets data structure.
        vector<uint64_t> rank(ploidy);
        vector<uint64_t> parent(ploidy);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<ploidy; i++) {
            disjointSets.make_set(i);
        }

        // Connect the similar pairs.
        for(const auto& p: similarPairs) {
            disjointSets.union_set(p.first, p.second);
        }

        // Gather the branches in each group.
        vector< vector<uint64_t> > groups(ploidy);
        for(uint64_t i=0; i<ploidy; i++) {
            groups[disjointSets.find_set(i)].push_back(i);
        }

        // Each of the non-empty ones generates a branch group.
        for(const auto& group: groups) {
            if(not group.empty()) {
                branchGroups.push_back(group);
            }
        }
    }

    if(debug) {
        cout << "Found the following similar groups:" << endl;
        for(const auto& branchGroup: branchGroups) {
            for(const uint64_t i: branchGroup) {
                cout << i << " ";
            }
            cout << endl;
        }
    }



    // Each branch group generates a new branch, if it has enough coverage.
    for(const auto& branchGroup: branchGroups) {

        // If this branch group consists of a single branch, don't do anything.
        if(branchGroup.size() == 1) {
            continue;
        }

        if(debug) {
            cout << "Working on branch group";
            for(const uint64_t i: branchGroup) {
                cout << " " << i;
            }
            cout << endl;
        }

        // Create an AnchorPair using all of the usableOrientedReadIds
        // for the branches in this group.
        AnchorPair newAnchorPair;
        newAnchorPair.anchorIdA = anchorId0;
        newAnchorPair.anchorIdB = anchorId1;
        for(const uint64_t i: branchGroup) {
            std::ranges::copy(usableOrientedReadIds[i], back_inserter(newAnchorPair.orientedReadIds));
        }
        deduplicate(newAnchorPair.orientedReadIds);
        if(debug) {
            cout << "The new branch has coverage " << newAnchorPair.size() << endl;
        }

        if(newAnchorPair.size() < options.bubbleCleanupMinCommonCount) {
            if(debug) {
                cout << "Coverage for this branch group is too low, no new branch generated." << endl;
            }
            continue;
        }

        edge_descriptor e;
        tie(e, ignore) = add_edge(bubble.v0, bubble.v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.emplace_back(newAnchorPair, newAnchorPair.getAverageOffset(anchors));
        if(debug) {
            cout << "Generated a new branch " << edge.id << " by combining the branches in this group." << endl;
        }

        // Now we can delete the old branches in this group.
        for(const uint64_t i: branchGroup) {
            if(debug) {
                cout << "Removing " << assemblyGraph[bubble.edges[i]].id << endl;
            }
            boost::remove_edge(bubble.edges[i], assemblyGraph);
        }
    }

    SHASTA_ASSERT(out_degree(bubble.v0, assemblyGraph) > 0);
    SHASTA_ASSERT(in_degree(bubble.v1, assemblyGraph) > 0);

    return true;
}



// Analyze a Bubble and finds pairs of "similar" branches.
// See analyzeSimilarSequences.h for more information.
bool AssemblyGraph::analyzeBubble(
    const Bubble& bubble,
    const vector<uint64_t> minRepeatCount,
    vector< pair<uint64_t, uint64_t> >& similarPairs
    ) const
{
    const AssemblyGraph& assemblyGraph = *this;
    using shasta::Base;

    const bool debug = false; // (assemblyGraph[bubble.edges.front()].id == 130581);
    if(debug) {
        cout << "Analyzing bubble";
        for (const edge_descriptor e: bubble.edges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;
    }

    SHASTA_ASSERT(bubble.edges.size() > 1);

    // Gather the sequences of all the sides of this bubble
    vector< vector<Base> > sequences;
    for(const edge_descriptor e: bubble.edges) {
        sequences.emplace_back();
        SHASTA_ASSERT(assemblyGraph[e].wasAssembled);
        assemblyGraph[e].getSequence(sequences.back());
    }

    if(debug) {
        cout << "Bubble sequences:" << endl;
        for(uint64_t i=0; i<bubble.edges.size(); i++) {
            const edge_descriptor e = bubble.edges[i];
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            cout << ">" << edge.id <<  " " << sequences[i].size() << endl;
            copy(sequences[i], ostream_iterator<Base>(cout));
            cout << endl;
        }
    }

    // Loop over pairs of bubble edges.
    ostream html(0);
    similarPairs.clear();
    for(uint64_t i0=0; i0<bubble.edges.size()-1; i0++) {
        const vector<Base>& sequence0 = sequences[i0];
        for(uint64_t i1=i0+1; i1<bubble.edges.size(); i1++) {
            const vector<Base>& sequence1 = sequences[i1];
             if(debug) {
                cout << "Checking " << assemblyGraph[bubble.edges[i0]].id <<
                    " against " << assemblyGraph[bubble.edges[i1]].id << endl;
            }

            if(areSimilarSequences(sequence0, sequence1, minRepeatCount, html)) {
                similarPairs.emplace_back(i0, i1);
            }
        }

    }

    return false;
}



// Compress linear chains of edges into a single edge.
uint64_t AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;
    uint64_t compressCount = 0;

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



        // Compact debug output.
        if(compressDebugLevel >= 1) {
            cout << "Compress";
            for(const edge_descriptor e: chain) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << " into " << edgeNew.id << endl;
        }



        // Detailed debug output.
        if(compressDebugLevel >= 2) {
            uint64_t stepCount = 0;
            for(const edge_descriptor e: chain) {
                const AssemblyGraphEdge& edge = assemblyGraph[e];
                const uint64_t stepBegin = stepCount;
                const uint64_t stepEnd = stepBegin + edge.size();
                cout << edge.id << " " << edge.size() << " steps become " <<
                    edgeNew.id << " steps " << stepBegin << "-" << stepEnd << endl;
                stepCount = stepEnd;
            }
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

        ++compressCount;
    }

    return compressCount;
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



uint64_t AssemblyGraph::detangleVertices(Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<options.detangleMaxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleVerticesIteration(detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            changeCount += compress();
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
    return detangleLowLevel(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleEdges(Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<options.detangleMaxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleEdgesIteration(detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            changeCount += compress();
            cout << "Detangle edges iteration " << iteration << ": " << iterationChangeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleEdgesIteration(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    // cout << "Edge detangling begins." << endl;

    // Gather edges on which we will attempt detangling.
    // Each generates a tangle with just two vertices.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const uint64_t in0 = in_degree(v0, assemblyGraph);
        const uint64_t out0 = out_degree(v0, assemblyGraph);
        const uint64_t in1 = in_degree(v1, assemblyGraph);
        const uint64_t out1 = out_degree(v1, assemblyGraph);

        const bool isTangleEdge = (in0 > 1) and (out0 == 1) and (in1 == 1) and (out1 > 1);

        // For cross edges, don't use in-degree and out-degree.
        // This is needed to make sure that bubble edges are not classified as cross-edges.
        const bool isCrossEdge =
            (countDistinctTargetVertices(v0) > 1) and
            (countDistinctSourceVertices(v1) > 1);

        const bool isShort = assemblyGraph[e].offset() <= options.detangleMaxCrossEdgeLength;

        // Tangle edges without length limitations.
        // Cross edges limited by length.
        if(isTangleEdge or(isShort and isCrossEdge)) {
            detanglingCandidates.emplace_back(vector<vertex_descriptor>({v0, v1}));
        }
    }

    // Do the detangling.
    return detangleLowLevel(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleLowLevel(
    const vector< vector<vertex_descriptor> >& detanglingCandidates,
    Detangler& detangler)
{
    const bool debug = false;

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
        Tangle1 tangle(assemblyGraph, tangleVertices);
        if(debug) {
            const TangleMatrix1& tangleMatrix = tangle.tangleMatrix();
            cout << "Tangle with " << tangleMatrix.entrances.size() << " entrances and " <<
                tangleMatrix.exits.size() << " exits." << endl;

            cout << "Entrances:";
            for(const edge_descriptor e: tangleMatrix.entrances) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;

            cout << "Exits:";
            for(const edge_descriptor e: tangleMatrix.exits) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;

            for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
                for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                    cout << tangleMatrix.tangleMatrix[iEntrance][iExit] << " ";
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

    // The detangling process can generate empty edges. Remove them.
    removeEmptyEdges();

    // cout << "Attempted detangling for " << attemptCount << " tangles." << endl;
    // cout << "Detangling was successful for " << successCount << " tangles." << endl;
    return successCount;
}



uint64_t AssemblyGraph::detangleHighLevel(Detangler& detangler)
{
    performanceLog << timestamp << "AssemblyGraph::detangle begins." << endl;

    const uint64_t verticesChangeCount = detangleVertices(detangler);
    const uint64_t edgesChangeCount = detangleEdges(detangler);
    const uint64_t shortTanglesChangeCount = detangleShortTangles(detangler);

    const uint64_t changeCount = verticesChangeCount + edgesChangeCount + shortTanglesChangeCount;

    performanceLog << timestamp << "AssemblyGraph::detangle ends." << endl;
    return changeCount;
}



uint64_t AssemblyGraph::prune()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    uint64_t totalPruneCount = 0;

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
            if(debug) {
                cout << "Pruned " << assemblyGraph[e].id << ", offset " << assemblyGraph[e].offset() << endl;
            }
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

        totalPruneCount += pruneCount;
    }

    return totalPruneCount;
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeCsv(csv);
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(ostream& csv) const
{
    const AssemblyGraph& assemblyGraph = *this;

    csv << "Segment,Number of steps,Average coverage,Estimated length,Actual length\n";
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        csv <<
            edge.id << "," <<
            edge.size() << "," <<
            coverage << "," <<
            edge.offset() << ",";
        if(edge.wasAssembled) {
            csv << edge.sequenceLength();
        }
        csv << ",\n";
    }
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
            orderById);

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



uint64_t AssemblyGraph::phaseSuperbubbleChains(
    bool onlyConsiderInjective,
    bool onlyConsiderPermutation)
{
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains begins." << endl;

    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    data.onlyConsiderInjective = onlyConsiderInjective;
    data.onlyConsiderPermutation = onlyConsiderPermutation;
    data.superbubbleChains = make_shared< vector<SuperbubbleChain> >();
    vector<SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);
    data.totalChangeCount = 0;

    // Find superbubbles.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    writeSuperbubbles(superbubbles, "Superbubbles-WithOverlaps.csv");
    removeContainedSuperbubbles(superbubbles);
    cout << "Found " << superbubbles.size() << " non-overlapping superbubbles." << endl;
    writeSuperbubbles(superbubbles, "Superbubbles.csv");
    writeSuperbubblesForBandage(superbubbles, "Superbubbles-Bandage.csv");

    // Find Superbubble chains.
    findSuperbubbleChains(superbubbles, superbubbleChains);
    cout << "Found " << superbubbleChains.size() << " superbubble chains." << endl;
    writeSuperbubbleChains(superbubbleChains, "SuperbubbleChains.csv");
    writeSuperbubbleChainsForBandage(superbubbleChains, "SuperbubbleChains-Bandage.csv");

    // Phase them.
    findOrientedReadEdgeInformation();
#if 0
    uint64_t changeCount = 0;
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
        changeCount += superbubbleChain.phase1(
            *this,
            superbubbleChainId,
            options.detangleMinCoverage,
            onlyConsiderInjective,
            onlyConsiderPermutation);
    }
#endif
    setupLoadBalancing(superbubbleChains.size(), 1);
    runThreads(&AssemblyGraph::phaseSuperbubbleChainsThreadFunction, options.threadCount);
    data.superbubbleChains = 0;
    uint64_t changeCount = data.totalChangeCount;
    clearOrientedReadEdgeInformation();

    changeCount += compress();
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains ends." << endl;

    return changeCount;
}



void AssemblyGraph::phaseSuperbubbleChainsThreadFunction(uint64_t /* threadId */)
{
    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    const bool onlyConsiderInjective = data.onlyConsiderInjective;
    const bool onlyConsiderPermutation = data.onlyConsiderPermutation;
    vector<SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all superbubble chains assigned to this batch.
        for(uint64_t superbubbleChainId=begin; superbubbleChainId<end; superbubbleChainId++) {
            SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
            const uint64_t changeCount = superbubbleChain.phase1(
                *this,
                superbubbleChainId,
                options.detangleMinCoverage,
                onlyConsiderInjective,
                onlyConsiderPermutation);
            __sync_fetch_and_add(&data.totalChangeCount, changeCount);
        }
    }
}



// Simplify Superbubbles by turning them into bubbles via clustering
// of oriented read journeys.
uint64_t AssemblyGraph::simplifySuperbubbles(Detangler& detangler)
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

    uint64_t simplifiedCount = 0;
    for(const Superbubble& superbubble: superbubbles) {
        if(not superbubble.isBubble()) {

            // First try simplify by detangling.
            bool success = simplifySuperbubbleByDetangling(superbubble, detangler);

            // If that did not work, try simplify by clustering.
            if(not success) {
                simplifySuperbubbleByClustering(superbubble,
                    options.simplifySuperbubbleMinCoverage,
                    options.simplifySuperbubbleMaxOffset);
            }

            // If one of the above was successful, increment our counter.
            if(success) {
                ++simplifiedCount;
            }
        }
    }

    // The detangling process can generate empty edges. Remove them.
    removeEmptyEdges();

    return simplifiedCount;
}



bool AssemblyGraph::simplifySuperbubbleByDetangling(
    const Superbubble& superbubble,
    Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    const bool debug = false;

    Tangle1 tangle(assemblyGraph, superbubble.internalVertices);

    if(debug) {
        cout << "simplifySuperbubbleByDetangling begins." << endl;

        const TangleMatrix1& tangleMatrix = tangle.tangleMatrix();
        cout << "Tangle with " << tangleMatrix.entrances.size() << " entrances and " <<
            tangleMatrix.exits.size() << " exits." << endl;

        cout << "Entrances:";
        for(const edge_descriptor e: tangleMatrix.entrances) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        cout << "Exits:";
        for(const edge_descriptor e: tangleMatrix.exits) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
            for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                cout << tangleMatrix.tangleMatrix[iEntrance][iExit] << " ";
            }
        }
        cout << endl;
    }

    const bool success = detangler(tangle);
    if(success) {
        if(debug) {
            cout << "Detangle was successful." << endl;
        }
    } else {
        if(debug) {
            cout << "Detangle failed." << endl;
        }
    }

    return success;

}



bool AssemblyGraph::simplifySuperbubbleByClustering(
    const Superbubble& superbubble,
    uint64_t minCoverage,
    uint64_t maxOffset)
{
    AssemblyGraph& assemblyGraph = *this;

    const AnchorId anchorIdA = assemblyGraph[superbubble.sourceVertex].anchorId;
    const AnchorId anchorIdB = assemblyGraph[superbubble.targetVertex].anchorId;

    const bool debug = false;



    if(debug) {
        cout << "Working on a superbubble with source vertex " <<
            assemblyGraph[superbubble.sourceVertex].id <<
            " and target vertex " << assemblyGraph[superbubble.targetVertex].id << endl;

        cout << "Source edges:";
        for(const edge_descriptor e: superbubble.sourceEdges) {
            cout << " " << assemblyGraph[e].id << " ";
        }
        cout << endl;

        cout << "Target edges:";
        for(const edge_descriptor e: superbubble.targetEdges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        cout << "This superbubble consists of the following " <<
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
        cout << "The initial anchor pair " <<
            anchorIdToString(anchorIdA) << " " <<
            anchorIdToString(anchorIdB) <<
            " has " << anchorPair.orientedReadIds.size() <<
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
        cout << "Found " << newAnchorPairs.size() << " split AnchorPairs:" << endl;
        for(const AnchorPair& newAnchorPair: newAnchorPairs) {
            cout << "AnchorPair with coverage " << newAnchorPair.orientedReadIds.size() <<
                ", offset " << newAnchorPair.getAverageOffset(anchors) << endl;
        }

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
            localAssembly.run(false, options.maxAbpoaLength);
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
        return false;
    }

    // If any of these AnchorPairs have a long offset, don't do it.
    // This only works well for small superbubbles, and for long superbubbles
    // it can destroy correct sequence.
    for(const AnchorPair& newAnchorPair: newAnchorPairs) {
        if(newAnchorPair.getAverageOffset(anchors) > maxOffset) {
            if(debug) {
                cout << "Skipping this superbubble due to large offset." << endl;
            }
            return false;
        }
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

    return true;
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
void AssemblyGraph::computeJourneys()
{
    const uint64_t orientedReadCount = journeys.size();

    AssemblyGraph& assemblyGraph = *this;


    // Compute uncompressed journeys for all oriented reads.
    class AssemblyGraphJourneyEntry {
    public:
        edge_descriptor e;
        uint64_t stepId;

        // These are positions in the standard Journeys on the Anchors.
        uint64_t positionInJourneyA;
        uint64_t positionInJourneyB;

        bool operator<(const AssemblyGraphJourneyEntry& that) const {
            return positionInJourneyA + positionInJourneyB <
                that.positionInJourneyA + that.positionInJourneyB;
        }
    };
    vector< vector<AssemblyGraphJourneyEntry> > assemblyGraphJourneys(orientedReadCount);

    // Loop over AssemblyGraph edges.
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


    // Sort the journeys.
    for(vector<AssemblyGraphJourneyEntry>& v: assemblyGraphJourneys) {
        sort(v.begin(), v.end());
    }



    // Create the compressed journeys.
    // Here, we only consider transitions between AssemblyGraph edges.
    // Compressed journeys consisting of only one edge are considered empty.
    compressedJourneys.clear();
    compressedJourneys.resize(orientedReadCount);
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
        vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];

        for(uint64_t i1=0; i1<assemblyGraphJourney.size(); i1++) {
            const AssemblyGraphJourneyEntry& entry1 = assemblyGraphJourney[i1];
            const edge_descriptor e1 = entry1.e;

            if(i1 == 0) {
                compressedJourney.push_back(e1);
            } else {
                const uint64_t i0 = i1 - 1;
                const AssemblyGraphJourneyEntry& entry0 = assemblyGraphJourney[i0];
                const edge_descriptor e0 = entry0.e;
                if(e1 != e0) {
                    compressedJourney.push_back(e1);
                }
            }
        }
        if(compressedJourney.size() == 1) {
            compressedJourney.clear();
        }
    }



    // Write out the compressed journeys.
    {
        ofstream csv("AssemblyGraphCompressedJourneys.csv");
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];
            csv << orientedReadId << ",";
            for(const edge_descriptor e: compressedJourney) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

    {
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
}




// Find appearances of each OrientedReadId in the AssemblyGraphSteps
// of each AssemblyGraphEdge.
// This also stores in each AssemblyGraphEdge::transitioningOrientedReadIds
// the OrientedReadIds that visit the edge and at least one other edge.
void AssemblyGraph::findOrientedReadEdgeInformation()
{
    AssemblyGraph& assemblyGraph = *this;
    const uint64_t orientedReadCount = journeys.size();

    // Gather the edges that each OrientedReadId appears in, allowing duplicates.
    vector< vector<edge_descriptor> > edgesByOrientedReadId(orientedReadCount);
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AssemblyGraphEdgeStep& step: edge) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                edgesByOrientedReadId[orientedReadId.getValue()].push_back(e);
            }
        }
    }

    // Deduplicate and count for each OrientedReadId, then store ordered by segmentId.
    orientedReadEdgeInformation.clear();
    orientedReadEdgeInformation.resize(orientedReadCount);
    OrientedReadEdgeInformationOrderById orientedReadEdgeInformationOrderById(assemblyGraph);
    vector<uint64_t> count;
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        vector<edge_descriptor>& x = edgesByOrientedReadId[orientedReadIdValue];
        deduplicateAndCount(x, count);
        vector<OrientedReadEdgeInformation>& y = orientedReadEdgeInformation[orientedReadIdValue];
        for(uint64_t i=0; i<x.size(); i++) {
            y.emplace_back(x[i], count[i]);
        }
        std::ranges::sort(y, orientedReadEdgeInformationOrderById);
    }



    // Store in each AssemblyGraphEdge the OrientedReadIds that visit
    // the edge and at least one other edge.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        assemblyGraph[e].transitioningOrientedReadIds.clear();
    }
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
        const vector<OrientedReadEdgeInformation>& v = orientedReadEdgeInformation[orientedReadIdValue];
        if(v.size() > 1) {
            for(const OrientedReadEdgeInformation& info: v) {
                assemblyGraph[info.e].transitioningOrientedReadIds.push_back(make_pair(orientedReadId, info.stepCount));
            }
        }
    }
}



void AssemblyGraph::clearOrientedReadEdgeInformation()
{
    AssemblyGraph& assemblyGraph = *this;

    orientedReadEdgeInformation.clear();

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.transitioningOrientedReadIds.clear();
    }
}



void AssemblyGraph::writeOrientedReadEdgeInformation()
{
    const AssemblyGraph& assemblyGraph = *this;
    const uint64_t orientedReadCount = journeys.size();

    {
        ofstream csv("OrientedReadEdgeInformaton.csv");
        csv << "OrientedReadId,SegmentId,StepCount,\n";

        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<OrientedReadEdgeInformation>& v = orientedReadEdgeInformation[orientedReadIdValue];
            for(const OrientedReadEdgeInformation x: v) {
                csv << orientedReadId << ",";
                csv << assemblyGraph[x.e].id << ",";
                csv << x.stepCount << ",\n";
            }
        }
    }



    {
        ofstream csv("TransitioningOrientedReads.csv");
        csv << "SegmentId,OrientedReadId,StepCount,\n";

        BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            for(const auto& p: edge.transitioningOrientedReadIds) {
                csv << edge.id << "," << p.first << "," << p.second << ",\n";
            }
        }

    }
}



void AssemblyGraph::computeExtendedTangleMatrix(
    vector<edge_descriptor>& entrances,
    vector<edge_descriptor>& exits,
    vector< vector<double> >& tangleMatrix,
    ostream& html
    ) const
{
    const AssemblyGraph& assemblyGraph = *this;
    SHASTA_ASSERT(not orientedReadEdgeInformation.empty());
    const bool debug = false;

    SHASTA_ASSERT(std::ranges::is_sorted(entrances, orderById));
    SHASTA_ASSERT(std::ranges::is_sorted(exits, orderById));



    if(debug) {
        cout << "Extended tangle matrix computation begins." << endl;

        cout << entrances.size() << " entrances:" << endl;
        for(uint64_t i=0; i<entrances.size(); i++) {
            cout << assemblyGraph[entrances[i]].id << " ";
        }
        cout << endl;

        cout << exits.size() << " exits:" << endl;
        for(uint64_t i=0; i<exits.size(); i++) {
            cout << assemblyGraph[exits[i]].id << " ";
        }
        cout << endl;

        for(uint64_t i=0; i<entrances.size(); i++) {
            cout << "Entrance " << assemblyGraph[entrances[i]].id << " has " << assemblyGraph[entrances[i]].transitioningOrientedReadIds.size() <<
                " oriented reads:" << endl;
            for(const auto& p: assemblyGraph[entrances[i]].transitioningOrientedReadIds) {
                cout << p.first << " " << p.second << " steps" << endl;
            }
        }
        for(uint64_t i=0; i<exits.size(); i++) {
            cout << "Exit " << assemblyGraph[exits[i]].id << " has " << assemblyGraph[exits[i]].transitioningOrientedReadIds.size() <<
                " oriented reads:" << endl;
            for(const auto& p: assemblyGraph[exits[i]].transitioningOrientedReadIds) {
                cout << p.first << " " << p.second << " steps" << endl;
            }
        }
    }


    // Find orientedReads that appear in one or more entrances.
    vector<OrientedReadId> entranceOrientedReadIds;
    for(uint64_t i=0; i<entrances.size(); i++) {
        for(const auto& p: assemblyGraph[entrances[i]].transitioningOrientedReadIds) {
            entranceOrientedReadIds.push_back(p.first);
        }
    }
    deduplicate(entranceOrientedReadIds);

    // Find orientedReads that appear in one or more exits.
    vector<OrientedReadId> exitOrientedReadIds;
    for(uint64_t i=0; i<exits.size(); i++) {
        for(const auto& p: assemblyGraph[exits[i]].transitioningOrientedReadIds) {
            exitOrientedReadIds.push_back(p.first);
        }
    }
    deduplicate(exitOrientedReadIds);


    // Find OrientedReadIds that appear in at least one entrance
    // and at least one exit.
    vector<OrientedReadId> orientedReadIds;
    std::set_intersection(
        entranceOrientedReadIds.begin(), entranceOrientedReadIds.end(),
        exitOrientedReadIds.begin(), exitOrientedReadIds.end(),
        back_inserter(orientedReadIds));

    if(debug) {
        cout << "Found " << orientedReadIds.size() <<
            " oriented reads that appear in at least one entrance and one exit." << endl;
    }


    // Gather the number of appearances in entrances and exits for each of these oriented reads.
    class OrientedReadInfo {
    public:
        vector<uint64_t> entranceStepCount;
        vector<uint64_t> exitStepCount;
    };
    std::map<OrientedReadId, OrientedReadInfo> m;
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        OrientedReadInfo& info = m[orientedReadId];
        info.entranceStepCount.resize(entrances.size(), 0);
        info.exitStepCount.resize(exits.size(), 0);
    }
    for(uint64_t i=0; i<entrances.size(); i++) {
        for(const auto& p: assemblyGraph[entrances[i]].transitioningOrientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            auto it = m.find(orientedReadId);
            if(it != m.end()) {
                const uint64_t count = p.second;
                OrientedReadInfo& info = it->second;
                info.entranceStepCount[i] = count;
            }
        }
    }
    for(uint64_t i=0; i<exits.size(); i++) {
        for(const auto& p: assemblyGraph[exits[i]].transitioningOrientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            auto it = m.find(orientedReadId);
            if(it != m.end()) {
                const uint64_t count = p.second;
                OrientedReadInfo& info = it->second;
                info.exitStepCount[i] = count;
            }
        }
    }

    if(debug) {
        cout << "OrientedReadInfo table:" << endl;
        for(const auto& p: m) {
            const OrientedReadId orientedReadId = p.first;
            const auto& info = p.second;
            cout << orientedReadId;
            for(uint64_t i=0; i<entrances.size(); i++) {
                cout << " " << info.entranceStepCount[i];
            }
            for(uint64_t i=0; i<exits.size(); i++) {
                cout << " " << info.exitStepCount[i];
            }
            cout << endl;
        }
    }



    if(html) {
        html << std::setprecision(2) << "<h3>Extended tangle matrix</h3>"
            "<h5>Oriented read contributions</h5>"
            "<table>"
            "<tr><th>Oriented<br>read<br>id";
        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<th>Entrance<br>" << assemblyGraph[entrances[i]].id;
        }
        for(uint64_t i=0; i<exits.size(); i++) {
            html << "<th>Exit<br>" << assemblyGraph[exits[i]].id;
        }
        html << "<th>Tangle<br>matrix";

        // Loop over oriented reads that contribute to this tangle matrix.
        for(const auto& p: m) {
            const OrientedReadId orientedReadId = p.first;
            const auto& info = p.second;
            const double entranceSum = double(std::accumulate(info.entranceStepCount.begin(), info.entranceStepCount.end(), 0UL));
            const double exitSum = double(std::accumulate(info.exitStepCount.begin(), info.exitStepCount.end(), 0UL));

            html << "<tr><th>" << orientedReadId;
            for(uint64_t i=0; i<entrances.size(); i++) {
                html << "<td class=centered>" << info.entranceStepCount[i];
            }
            for(uint64_t i=0; i<exits.size(); i++) {
                html << "<td class=centered>" << info.exitStepCount[i];
            }

            // Write the contribution of this oriented read to the total tangle matrix.
            html << "<td class=centered><table style='margin: 0 auto;'>";

            for(uint64_t i=0; i<entrances.size(); i++) {
                html << "<tr>";
                for(uint64_t j=0; j<exits.size(); j++) {
                    html << "<td class=centered>" << double(info.entranceStepCount[i] * info.exitStepCount[j]) / (entranceSum * exitSum);
                }
            }
            html << "</table>";
        }
        html << "</table>";
    }



    // Compute the extended tangle matrix.
    tangleMatrix.clear();
    tangleMatrix.resize(entrances.size(), vector<double>(exits.size(), 0));
    for(const auto& p: m) {
        const auto& info = p.second;
        const double entranceSum = double(std::accumulate(info.entranceStepCount.begin(), info.entranceStepCount.end(), 0UL));
        const double exitSum = double(std::accumulate(info.exitStepCount.begin(), info.exitStepCount.end(), 0UL));
        if(debug) {
            cout << "Tangle matrix contribution of " << p.first << ":";
        }
        for(uint64_t i=0; i<entrances.size(); i++) {
            for(uint64_t j=0; j<exits.size(); j++) {
                tangleMatrix[i][j] += double(info.entranceStepCount[i] * info.exitStepCount[j]) / (entranceSum * exitSum);
                if(debug) {
                    cout << " " << double(info.entranceStepCount[i] * info.exitStepCount[j]) / (entranceSum * exitSum);
                }
            }
        }
        if(debug) {
            cout << endl;
        }
    }

    if(debug) {
        cout << "Extended tangle matrix:" << endl;
        for(uint64_t i=0; i<entrances.size(); i++) {
            for(uint64_t j=0; j<exits.size(); j++) {
                cout << tangleMatrix[i][j] << " ";
            }
            cout << endl;
        }
    }



    if(html) {
        html << "<h5>Total tangle matrix</h5><table><tr><th>";
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<th>" << assemblyGraph[exits[j]].id;
        }

        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<tr><th>" << assemblyGraph[entrances[i]].id;
            for(uint64_t j=0; j<exits.size(); j++) {
                html << "<td class=centered>" << tangleMatrix[i][j];
            }
        }

        html << "</table>";
    }
}



// Search starting at a given edge (segment) and moving in the specified direction.
void AssemblyGraph::search(
    edge_descriptor eStart,
    uint64_t direction) const
{
    const AssemblyGraph& assemblyGraph = *this;

    const bool debug = false;

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
            0,
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
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
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

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, 0, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                cout << assemblyGraph[exits[iExit]].id << " " <<
                    tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() << endl;
            }
        }

        // Counts exits by type.
        // Significant: coverage >= highCoverageThreshold
        // Insignificant: coverage <= lowCoverageThreshold
        // Ambiguous: lowCoverageThreshold< coverage < highCoverageThreshold.
        uint64_t significantCount = 0;
        uint64_t insignificantCount = 0;
        uint64_t ambiguousCount = 0;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size();
            if(coverage <= lowCoverageThreshold) {
                ++insignificantCount;
            } else if(coverage >= highCoverageThreshold) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }

        // Only keep going if we have exactly one significant exit and no ambiguous exits.
        if(significantCount != 1) {
            return;
        }
        if(ambiguousCount > 0) {
            return;
        }

        // Find the one and only significant exit.
        edge_descriptor e1;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() >= highCoverageThreshold) {
                e1 = exits[iExit];
                break;
            }
        }

        // Check for loops.
        if(e1 == eStart) {
            break;
        }
        bool loopDetected = false;
        for(const edge_descriptor e: edges) {
            if(e == e1) {
                loopDetected = true;
                break;
            }
        }
        if(loopDetected) {
            break;
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::backwardLocalSearch(
    edge_descriptor eStart,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
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

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, 0, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                cout << assemblyGraph[entrances[iEntrance]].id << " " <<
                    tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() << endl;
            }
        }

        // Counts exits by type.
        // Significant: coverage >= highCoverageThreshold
        // Insignificant: coverage <= lowCoverageThreshold
        // Ambiguous: lowCoverageThreshold< coverage < highCoverageThreshold.
        uint64_t significantCount = 0;
        uint64_t insignificantCount = 0;
        uint64_t ambiguousCount = 0;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size();
            if(coverage <= lowCoverageThreshold) {
                ++insignificantCount;
            } else if(coverage >= highCoverageThreshold) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }

        // Only keep going if we have exactly one significant exit and no ambiguous exits.
        if(significantCount != 1) {
            return;
        }
        if(ambiguousCount > 0) {
            return;
        }

        // Find the one and only good entrance.
        edge_descriptor e1;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() >= highCoverageThreshold) {
                e1 = entrances[iEntrance];
                break;
            }
        }

        // Check for loops.
        if(e1 == eStart) {
            break;
        }
        bool loopDetected = false;
        for(const edge_descriptor e: edges) {
            if(e == e1) {
                loopDetected = true;
                break;
            }
        }
        if(loopDetected) {
            break;
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::testLocalSearch(
    uint64_t edgeIdStart,
    uint64_t direction,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold) const
{
    const AssemblyGraph& assemblyGraph = *this;

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
        forwardLocalSearch(eStart, lowCoverageThreshold, highCoverageThreshold, edges);
    } else {
        backwardLocalSearch(eStart, lowCoverageThreshold, highCoverageThreshold, edges);
    }

    cout << "Found " << edges.size() << " edges:" << endl;
    for(const edge_descriptor e: edges) {
        cout << assemblyGraph[e].id << " ";
    }
    cout << endl;


}



void AssemblyGraph::createSearchGraph(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    using shasta::SearchGraph;
    AssemblyGraph& assemblyGraph = *this;

    // Create the SearchGraph.
    SearchGraph searchGraph(*this, lowCoverageThreshold, highCoverageThreshold);

    // Compute connected components.
    vector<SearchGraph> components;
    searchGraph.computeConnectedComponents(components);
    // cout << "Found " << components.size() << " non-trivial connected components of the SearchGraph." << endl;



    // Process each connected component separately.
    ofstream csv("SearchGraph-Chains.csv");
    ofstream csvBandage("SearchGraph-Bandage.csv");
    csvBandage << "Segment,Color\n";
    vector< vector<vertex_descriptor> > chains;
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        // cout << "Working on component " << componentId << " of " << components.size() << endl;
        SearchGraph& component = components[componentId];
        transitiveReductionAny(component);
        component.removeBranches();
        // component.writeGraphviz("SearchGraph-" + to_string(componentId) + ".dot");

        findLinearVertexChains(component, chains);

        // For each chain, generate a new AssemblyGraphEdge,
        // without connecting it to the rest of the AssemblyGraph for now.
        for(const vector<vertex_descriptor>& chain: chains) {
            if(chain.size() < 2) {
                continue;
            }

            // Get some information about the first and last AssemblyGraph edge in the chain.
            const SearchGraph::vertex_descriptor sv0 = chain.front();
            const SearchGraph::vertex_descriptor sv1 = chain.back();
            const edge_descriptor e0 = component[sv0].e;
            const edge_descriptor e1 = component[sv1].e;

            const vertex_descriptor v0 = source(e0, assemblyGraph);
            const vertex_descriptor v1 = target(e1, assemblyGraph);

#if 0
            // Create new vertices for the new AssemblyGraphEdge, so it
            // will stay isolated from the rest of the AssemblyGraph for now.
            const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
            const AssemblyGraphEdge& edge1 = assemblyGraph[e1];
            const AnchorId anchorId0 = edge0.firstAnchorId();
            const AnchorId anchorId1 = edge1.lastAnchorId();
            const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
            const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
#endif

            // Create the new AssemblyGraphEdge.
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
            csvBandage << edgeNew.id << ",Red\n";



            // Concatenate the AssemblyGraphEdges of the chain, adding bridge steps where needed.
            for(uint64_t i=0; i<chain.size(); i++) {
                const SearchGraph::vertex_descriptor v = chain[i];
                const edge_descriptor e = component[v].e;
                const AssemblyGraphEdge& edge = assemblyGraph[e];

                // If necessary, add an AssemblyGraphEdgeStep to bridge.
                if(i != 0) {
                    const SearchGraph::vertex_descriptor vPrevious = chain[i - 1];
                    const edge_descriptor ePrevious = component[vPrevious].e;
                    const AssemblyGraphEdge& edgePrevious = assemblyGraph[ePrevious];
                    const AnchorId anchorIdPrevious = edgePrevious.lastAnchorId();
                    const AnchorId anchorIdNext = edge.firstAnchorId();
                    if(anchorIdPrevious != anchorIdNext) {
                        const AnchorPair bridgeAnchorPair = anchors.bridge(
                            edgePrevious.back().anchorPair,
                            edge.front().anchorPair,
                            options.aDrift, options.bDrift);
                        const uint64_t offset = bridgeAnchorPair.getAverageOffset(anchors);
                        edgeNew.push_back(AssemblyGraphEdgeStep(bridgeAnchorPair, offset));
                    }
                }

                // Now we can append this edge to the new edge.
                copy(edge.begin(), edge.end(), back_inserter(edgeNew));

            }

            // Now remove the edges of the chain.
            for(const SearchGraph::vertex_descriptor sv: chain) {
                const edge_descriptor e = component[sv].e;
                boost::remove_edge(e, assemblyGraph);
            }


            // Write this chain to the csv file.
            csv << edgeNew.id << ",";
            for(const vertex_descriptor v: chain) {
                const AssemblyGraph::edge_descriptor e = searchGraph[v].e;
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }
}


// Gather information on the oriented reads that appear
// in the AssemblyGraphSteps of an edge.
void AssemblyGraph::gatherOrientedReadInformationOnEdge(
    edge_descriptor e,
    vector<OrientedReadEdgeInfo>& orientedReadEdgeInfos) const
{
    const AssemblyGraph& assemblyGraph = *this;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    // A work vector to store information on each appearance of an OrientedRead
    // on a step of this edge.
    class WorkInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t positionInJourneyA;
        uint32_t positionInJourneyB;
        WorkInfo(
            OrientedReadId orientedReadId,
            uint32_t positionInJourneyA,
            uint32_t positionInJourneyB) :
            orientedReadId(orientedReadId),
            positionInJourneyA(positionInJourneyA),
            positionInJourneyB(positionInJourneyB)
        {}
        bool operator<(const WorkInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<WorkInfo> workInfos;



    // Loop over steps of this edge.
    for(const AssemblyGraphEdgeStep& step: edge) {
        const AnchorPair& anchorPair = step.anchorPair;
        const Anchor anchorA = anchors[anchorPair.anchorIdA];
        const Anchor anchorB = anchors[anchorPair.anchorIdB];

        // Triple loop over OrientedReadIds of anchorPair, anchorA, anchorB.
        auto itA = anchorA.begin();
        auto itB = anchorB.begin();
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            while(itA->orientedReadId < orientedReadId) {
                ++itA;
                SHASTA_ASSERT(itA != anchorA.end());
            }
            SHASTA_ASSERT(itA->orientedReadId == orientedReadId);
            while(itB->orientedReadId < orientedReadId) {
                ++itB;
                SHASTA_ASSERT(itB != anchorB.end());
            }
            SHASTA_ASSERT(itB->orientedReadId == orientedReadId);
            workInfos.emplace_back(orientedReadId, itA->positionInJourney, itB->positionInJourney);
        }
    }

    // Sort the WorkInfos by OrientedReadId.
    sort(workInfos.begin(), workInfos.end());



    // Each streak with the same OrientedReadId generates an OrientedReadEdgeInfo.
    orientedReadEdgeInfos.clear();
    for(auto streakBegin=workInfos.begin(); streakBegin!=workInfos.end(); /* Update later */) {
        const OrientedReadId orientedReadId = streakBegin->orientedReadId;

        auto streakEnd = streakBegin;
        while((streakEnd!= workInfos.end()) and (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }

        // Loop over this streak.
        OrientedReadEdgeInfo orientedReadEdgeInfo;
        orientedReadEdgeInfo.orientedReadId = orientedReadId;
        orientedReadEdgeInfo.stepCount = streakEnd - streakBegin;
        orientedReadEdgeInfo.minPositionInJourney = std::numeric_limits<uint32_t>::max();
        orientedReadEdgeInfo.maxPositionInJourney = std::numeric_limits<uint32_t>::max();
        for(auto it=streakBegin; it!=streakEnd; it++) {
            orientedReadEdgeInfo.minPositionInJourney = min(orientedReadEdgeInfo.minPositionInJourney, it->positionInJourneyA);
            orientedReadEdgeInfo.maxPositionInJourney = max(orientedReadEdgeInfo.minPositionInJourney, it->positionInJourneyB);
        }
        orientedReadEdgeInfos.push_back(orientedReadEdgeInfo);

        // Prepare to process the next streak.
        streakBegin = streakEnd;
    }
}



// The detangling process can generate empty edges (edges without steps).
// This removes them by collapsing the vertices they join.
void AssemblyGraph::removeEmptyEdges()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // We need to find groups of vertices that need to be collapsed together.
    // Usually it will be just two vertices to be collapsed together,
    // but it could be larger groups if there are adjacent empty edges.
    // So we need to find the connected components generated by the empty edges.

    // Map vertices to integer.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }

    // Initialize the disjoint sets data structure.
    const uint64_t n = vertexIndexMap.size();
    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over the empty edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const uint64_t i0 = vertexIndexMap[v0];
            const uint64_t i1 = vertexIndexMap[v1];
            disjointSets.union_set(i0, i1);
            edgesToBeRemoved.push_back(e);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > groups(n);
    for(uint64_t i=0; i<n; i++) {
        groups[disjointSets.find_set(i)].push_back(i);
    }

    // Collapse each group into a single vertex.
    for(const vector<uint64_t>& group: groups) {
        if(group.size() < 2) {
            continue;
        }
        if(debug) {
            cout << "Found a group of " << group.size() << " vertices to be collapsed:";
            for(const uint64_t vertexIndex: group) {
                const vertex_descriptor v = vertexTable[vertexIndex];
                cout << " " << assemblyGraph[v].id << "/" << anchorIdToString(assemblyGraph[v].anchorId);
            }
            cout << endl;
        }

        // Create the collapsed vertex.
        const AnchorId anchorId = assemblyGraph[vertexTable[group.front()]].anchorId;
        const vertex_descriptor vNew = add_vertex(
             AssemblyGraphVertex(anchorId, assemblyGraph.nextVertexId++), assemblyGraph);



        // For each vertex in this group, reroute all incoming/outgoing non-empty edges
        // to/from the new vertex.
        for(const uint64_t vertexIndex: group) {
            const vertex_descriptor v = vertexTable[vertexIndex];

            // Loop over non-empty incoming edges.
            BGL_FORALL_INEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = source(eOld, assemblyGraph);

                // Create the new edge, with target vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(u, vNew, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }


            // Loop over non-empty outgoing edges.
            BGL_FORALL_OUTEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = target(eOld, assemblyGraph);

                // Create the new edge, with source vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(vNew, u, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }
        }
    }

    deduplicate(edgesToBeRemoved);
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }

    // Remove any vertices that were left isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}



// Count the distinct target vertices among
// all the out-edges of a given vertex.
// So this counts the number of parallel edge sets
// outgoing from the given vertex.
uint64_t AssemblyGraph::countDistinctTargetVertices(vertex_descriptor v) const
{
    const AssemblyGraph& assemblyGraph = *this;

    boost::container::small_vector<vertex_descriptor, 6> targetVertices;
    BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
        targetVertices.push_back(target(e, assemblyGraph));
    }

    std::ranges::unique(targetVertices);
    return targetVertices.size();
}



// Count the distinct source vertices among
// all the in-edges of a given vertex.
// So this counts the number of parallel edge sets
// incoming into the given vertex.
uint64_t AssemblyGraph::countDistinctSourceVertices(vertex_descriptor v) const
{
    const AssemblyGraph& assemblyGraph = *this;

    boost::container::small_vector<vertex_descriptor, 6> sourceVertices;
    BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
        sourceVertices.push_back(source(e, assemblyGraph));
    }

    std::ranges::unique(sourceVertices);
    return sourceVertices.size();
}



// Detangle short tangles.
uint64_t AssemblyGraph::detangleShortTangles(Detangler& detangler)
{
    const AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t shortTangleMaxStepCount = 10;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }
    const uint64_t n = vertexIndexMap.size();

    // Initialize the disjoint set data structures.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }



    // Compute connected components of the AssemblyGraph using only short edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {

        // If not a short edge, skip it.
        if(assemblyGraph[e].size() > shortTangleMaxStepCount) {
            continue;
        }

        // Find the index corresponding to the source vertex.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const auto it0 = vertexIndexMap.find(v0);
        SHASTA_ASSERT(it0 != vertexIndexMap.end());
        const uint64_t i0 = it0->second;

        // Find the index corresponding to the target vertex.
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const auto it1 = vertexIndexMap.find(v1);
        SHASTA_ASSERT(it1 != vertexIndexMap.end());
        const uint64_t i1 = it1->second;

        // Update the disjoint sets for this edge.
        disjointSets.union_set(i0, i1);
    }



    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > components(n);
    for(uint64_t i=0; i<n; i++) {
        components[disjointSets.find_set(i)].push_back(i);
    }



    // Each connected components with size greater than 1 generates a short tangle.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    for(const vector<uint64_t>& component: components) {
        if(component.size() <= 1) {
            continue;
        }

        // Gather the tangle vertices.
        vector<vertex_descriptor> tangleVertices;
        for(const uint64_t i: component) {
            tangleVertices.push_back(vertexTable[i]);
        }

        // Add it to our detangling candidates.
        detanglingCandidates.push_back(tangleVertices);
    }

    // Do the detangling.
    const uint64_t successCount = detangleLowLevel(detanglingCandidates, detangler);

    cout << "Of " << detanglingCandidates.size() <<
        " short tangles, " << successCount << " were detangled successfully." << endl;
    compress();

    return successCount;
}
