#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Anchor.hpp"
#include "Assembler.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "Base.hpp"
#include "CycleAvoider.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "diploidBayesianPhase.hpp"
#include "ExternalAnchors.hpp"
#include "extractKmer128.hpp"
#include "findConvergingVertex.hpp"
#include "LikelihoodRatioDetangler.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultithreadedObject.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing.hpp"
#include "ReadSummary.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
#include "testSubsetGraph.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta2, shasta2Module)
{

	// Class Options.
    class_<Options>(shasta2Module, "Options")

        // Constructor from the name of a configuration file.
        .def(pybind11::init<const string&>(),
            arg("configurationFileName") = "shasta2.conf")
        .def_readwrite("threadCount", &Options::threadCount)
        .def_readwrite("k", &Options::k)
        .def_readwrite("markerDensity", &Options::markerDensity)
        .def_readwrite("maxMarkerErrorRate", &Options::maxMarkerErrorRate)
        .def_readwrite("externalAnchorsName", &Options::externalAnchorsName)
        .def_readwrite("minAnchorCoverage", &Options::minAnchorCoverage)
        .def_readwrite("maxAnchorCoverage", &Options::maxAnchorCoverage)
        .def_readwrite("maxAnchorRepeatLength", &Options::maxAnchorRepeatLength)
        .def_readwrite("minAnchorGraphEdgeCoverage", &Options::minAnchorGraphEdgeCoverage)

		// Options defined in OptionsDefine.hpp
		#define SHASTA2_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
			.def_readwrite(#name, &Options::name)
		#include "OptionsDefine.hpp"
		#undef SHASTA2_OPTION_DEFINE
        ;



    // Class Assembler.
    class_<Assembler>(shasta2Module, "Assembler")

        // Constructor.
        .def(pybind11::init<const string&>(),
            "Assembler constructor.",
            arg("largeDataFileNamePrefix") = "Data/")

        // Reads.
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")
        .def("createReadLengthDistribution",
            &Assembler::createReadLengthDistribution)
        .def("accessReadLengthDistribution",
            &Assembler::accessReadLengthDistribution)
        .def("accessReadSummaries",
            &Assembler::accessReadSummaries)

        // K-mer checker.
        .def("createKmerChecker",
            &Assembler::createKmerChecker,
            arg("k"),
            arg("markerDensity"),
            arg("threadCount") = 0)
        .def("accessKmerChecker",
            &Assembler::accessKmerChecker)

         // Markers.
        .def("accessMarkers",
            &Assembler::accessMarkers)
        .def("createMarkers",
            &Assembler::createMarkers,
            arg("threadCount") = 0)

         // Marker k-mers.
        .def("createMarkerKmers", &Assembler::createMarkerKmers)
        .def("accessMarkerKmers",
            &Assembler::accessMarkerKmers)
        .def("computeMarkerErrorRates",
            &Assembler::computeMarkerErrorRates)

        // Anchors.
       .def("createAnchors",
           &Assembler::createAnchors,
           arg("minAnchorCoverage"),
           arg("maxAnchorCoverage"),
           arg("maxHomopolymerLength"),
           arg("threadCount") = 0)
       .def("readExternalAnchors", &Assembler::readExternalAnchors)
       .def("accessAnchors",
           &Assembler::accessAnchors,
           arg("writeAccess") = false)

       // Journeys.
      .def("createJourneys",
          &Assembler::createJourneys,
          arg("threadCount") = 0)
      .def("accessJourneys",
          &Assembler::accessJourneys)
      .def("storeAnchorGaps",
          &Assembler::storeAnchorGaps)
      .def("readFollowing",
          &Assembler::readFollowing)

      // AnchorGraph.
      .def("createAnchorGraph",
          &Assembler::createAnchorGraph)
      .def("accessAnchorGraph",
          &Assembler::accessAnchorGraph)
      .def("saveAnchorGraph",
          &Assembler::saveAnchorGraph)
      .def("anchorGraphTransitiveReduction",
          &Assembler::anchorGraphTransitiveReduction)
      .def("createSimpleAnchorGraph",
          &Assembler::createSimpleAnchorGraph)

      // AssemblyGraph.
      .def("createAssemblyGraph",
          &Assembler::createAssemblyGraph)
      .def("getAssemblyGraph",
          &Assembler::getAssemblyGraph, return_value_policy::reference)
    ;



    class_<AssemblyGraph> assemblyGraphClass(shasta2Module, "AssemblyGraph");
    assemblyGraphClass
        .def(pybind11::init<
            const Assembler&,
            const Options&,
            const AssemblyGraph&,
            const vector< vector<AssemblyGraph::edge_descriptor> >&
            >())
        .def_readwrite("compressDebugLevel", &AssemblyGraph::compressDebugLevel)
        .def("prune", &AssemblyGraph::prune)
        .def("bubbleCleanupIteration", &AssemblyGraph::bubbleCleanupIteration)
        .def("detangleVertices", &AssemblyGraph::detangleVertices)
        .def("detangleVerticesIteration", &AssemblyGraph::detangleVerticesIteration)
        .def("detangleEdges", &AssemblyGraph::detangleEdges)
        .def("detangleEdgesIteration", &AssemblyGraph::detangleEdgesIteration)
        .def("detangleShortTangles", &AssemblyGraph::detangleShortTangles)
        .def("detangleHighLevel", &AssemblyGraph::detangleHighLevel)
        .def("detangleLowLevel", &AssemblyGraph::detangleLowLevel)
        .def("detangleShortTangles", &AssemblyGraph::detangleShortTangles)
        .def("compress", &AssemblyGraph::compress)
        .def("removeEmptyEdges", &AssemblyGraph::removeEmptyEdges)
        .def("assembleAll", &AssemblyGraph::assembleAll)
        .def("clearSequence", &AssemblyGraph::clearSequence)
        .def("phaseSuperbubbleChains", &AssemblyGraph::phaseSuperbubbleChains)
        .def("simplifySuperbubbles", &AssemblyGraph::simplifySuperbubbles)
        .def("colorStrongComponents", &AssemblyGraph::colorStrongComponents)
        .def("write", &AssemblyGraph::write)
        .def("writeFasta", &AssemblyGraph::writeFasta)
        .def("computeJourneys", &AssemblyGraph::computeJourneys)
        .def("findOrientedReadEdgeInformation", &AssemblyGraph::findOrientedReadEdgeInformation)
        .def("writeOrientedReadEdgeInformation", &AssemblyGraph::writeOrientedReadEdgeInformation)
        .def("findEdgePairs", &AssemblyGraph::findEdgePairs)
        .def("testSearch", &AssemblyGraph::testSearch)
        .def("testLocalSearch", &AssemblyGraph::testLocalSearch)
        .def("createSearchGraph", &AssemblyGraph::createSearchGraph)
        .def("findAndConnectAssemblyPaths", &AssemblyGraph::findAndConnectAssemblyPaths)
        ;

    // Expose AssemblyGraph vertex_descriptor and edge_descriptor.
    class_<AssemblyGraph::vertex_descriptor>
        assemblyGraphVertexDescriptorClass(assemblyGraphClass, "AssemblyGraphVertexDescriptor");
    class_<AssemblyGraph::edge_descriptor>
        assemblyGraphEdgeDescriptorClass(assemblyGraphClass, "AssemblyGraphEdgeDescriptor");

    class_<AssemblyGraphPostprocessor>(shasta2Module, "AssemblyGraphPostprocessor",
        pybind11::base<AssemblyGraph>())
        .def_readonly("vertexMap", &AssemblyGraphPostprocessor::vertexMap)
        .def_readonly("edgeMap", &AssemblyGraphPostprocessor::edgeMap)
        .def("getId", &AssemblyGraphPostprocessor::getId)
        ;



    // Detangler classes.
    class_<Detangler>(shasta2Module, "Detangler")
        .def_readwrite("debug", &Detangler::debug)
        ;
    class_<LikelihoodRatioDetangler>(shasta2Module, "LikelihoodRatioDetangler", pybind11::base<Detangler>())
        .def(init<double, double, double, uint64_t, bool, bool, bool, bool>())
        ;



    // Class ReadFollowing.
    class_<ReadFollowing::Graph>(shasta2Module, "ReadFollowing")
        .def(pybind11::init<const AssemblyGraph&>())
        .def("writePath", &ReadFollowing::Graph::writePath)
        .def("writePaths", &ReadFollowing::Graph::writePaths)
        ;

    // Class ExternalAnchors.
    class_<ExternalAnchors>(shasta2Module, "ExternalAnchors")
        .def(pybind11::init<const string&>())
        .def("beginNewAnchor", &ExternalAnchors::beginNewAnchor)
        .def("addOrientedRead", &ExternalAnchors::addOrientedRead)
        ;



    // Non-member functions exposed to Python.
    shasta2Module.def("openPerformanceLog",
        openPerformanceLog
        );
    shasta2Module.def("testMultithreadedObject",
        testMultithreadedObject
        );
    shasta2Module.def("testMemoryMappedVector",
        testMemoryMappedVector
        );
    shasta2Module.def("testBase",
        testBase
        );
    shasta2Module.def("testShortBaseSequence",
        testShortBaseSequence
        );
    shasta2Module.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shasta2Module.def("testExtractKmer128",
        testExtractKmer128
        );
    shasta2Module.def("testSplitRange",
        testSplitRange
        );
    shasta2Module.def("testDeduplicateAndCount",
        testDeduplicateAndCount
        );
    shasta2Module.def("testBitReversal",
        testBitReversal
        );
    shasta2Module.def("mappedCopy",
        mappedCopy
        );
    shasta2Module.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shasta2Module.def("testDiploidBayesianPhase",
        testDiploidBayesianPhase
        );
    shasta2Module.def("testSubsetGraph",
        testSubsetGraph
        );
    shasta2Module.def("testFindConvergingVertex",
        testFindConvergingVertex
        );
    shasta2Module.def("anchorIdToString",
        anchorIdToString
        );
    shasta2Module.def("anchorIdFromString",
        anchorIdFromString
        );
    shasta2Module.def("testCycleAvoider",
        testCycleAvoider
        );
}

#endif


