#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "Base.hpp"
#include "ChiSquareDetangler.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "diploidBayesianPhase.hpp"
#include "extractKmer128.hpp"
#include "findConvergingVertex.hpp"
#include "globalMsa.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultithreadedObject.hpp"
#include "PermutationDetangler.hpp"
#include "performanceLog.hpp"
#include "ShortBaseSequence.hpp"
#include "SimpleDetangler.hpp"
#include "splitRange.hpp"
#include "testSpoa.hpp"
#include "testSubsetGraph.hpp"
#include "TrivialDetangler.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta2, shasta2Module)
{
    // Class AssemblerOptions::AssemblyGraphOptions.
    class_<AssemblerOptions::AssemblyGraphOptions>(shasta2Module, "AssemblyGraphOptions")
        .def_readonly("transitiveReductionThreshold", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionThreshold)
        .def_readonly("transitiveReductionA", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionA)
        .def_readonly("transitiveReductionB", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionB)
        .def_readonly("detangleMinCommonCoverage", &AssemblerOptions::AssemblyGraphOptions::detangleMinCommonCoverage)
        .def_readonly("detangleLowCoverageThreshold", &AssemblerOptions::AssemblyGraphOptions::detangleLowCoverageThreshold)
        .def_readonly("detangleHighCoverageThreshold", &AssemblerOptions::AssemblyGraphOptions::detangleHighCoverageThreshold)
        .def_readonly("detangleEpsilon", &AssemblerOptions::AssemblyGraphOptions::detangleEpsilon)
        .def_readonly("detangleMaxLogP", &AssemblerOptions::AssemblyGraphOptions::detangleMaxLogP)
        .def_readonly("detangleMinLogPDelta", &AssemblerOptions::AssemblyGraphOptions::detangleMinLogPDelta)
        ;

    // Class AssemblerOptions.
    class_<AssemblerOptions>(shasta2Module, "AssemblerOptions")

        // Constructor from the name of a configuration file.
        .def(pybind11::init<const string&>(),
            arg("configurationFileName") = "shasta2.conf")
        .def_readonly("threadCount", &AssemblerOptions::threadCount)
        .def_readonly("k", &AssemblerOptions::k)
        .def_readonly("markerDensity", &AssemblerOptions::markerDensity)
        .def_readonly("minAnchorCoverage", &AssemblerOptions::minAnchorCoverage)
        .def_readonly("maxAnchorCoverage", &AssemblerOptions::maxAnchorCoverage)
        .def_readonly("minAnchorGraphEdgeCoverageNear", &AssemblerOptions::minAnchorGraphEdgeCoverageNear)
        .def_readonly("minAnchorGraphEdgeCoverageFar", &AssemblerOptions::minAnchorGraphEdgeCoverageFar)
        .def_readonly("assemblyGraphOptions", &AssemblerOptions::assemblyGraphOptions)
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
        .def("createMarkerKmers",
            &Assembler::createMarkerKmers,
            arg("threadCount") = 0)
        .def("accessMarkerKmers",
            &Assembler::accessMarkerKmers)

        // Anchors.
       .def("createAnchors",
           &Assembler::createAnchors,
           arg("minAnchorCoverage"),
           arg("maxAnchorCoverage"),
           arg("maxHomopolymerLength"),
           arg("threadCount") = 0)
       .def("accessAnchors",
           &Assembler::accessAnchors,
           arg("writeAccess") = false)

       // Journeys.
      .def("createJourneys",
          &Assembler::createJourneys,
          arg("threadCount") = 0)
      .def("accessJourneys",
          &Assembler::accessJourneys)

      // AnchorGraph.
      .def("createAnchorGraph",
          &Assembler::createAnchorGraph)
      .def("accessAnchorGraph",
          &Assembler::accessAnchorGraph)

      // AssemblyGraph.
      .def("createAssemblyGraph",
          &Assembler::createAssemblyGraph,
          arg("assemblerOptions"),
          arg("threadCount") = 0)
      .def("getAssemblyGraph",
          &Assembler::getAssemblyGraph, return_value_policy::reference)
    ;



    class_<AssemblyGraph>(shasta2Module, "AssemblyGraph")
        .def("detangleVertices", &AssemblyGraph::detangleVertices)
        .def("detangleEdges", &AssemblyGraph::detangleEdges)
        .def("compress", &AssemblyGraph::compress)
        .def("assembleAll", &AssemblyGraph::assembleAll)
        .def("analyzeSuperbubbles", &AssemblyGraph::analyzeSuperbubbles)
        .def("colorStrongComponents", &AssemblyGraph::colorStrongComponents)
        .def("write", &AssemblyGraphPostprocessor::write)
        ;

    class_<AssemblyGraphPostprocessor>(shasta2Module, "AssemblyGraphPostprocessor",
        pybind11::base<AssemblyGraph>());



    // Detangler classes.
    class_<Detangler>(shasta2Module, "Detangler")
        .def_readwrite("debug", &Detangler::debug)
        ;
    class_<TrivialDetangler>(shasta2Module, "TrivialDetangler", pybind11::base<Detangler>())
        .def(init<uint64_t>())
        ;
    class_<SimpleDetangler>(shasta2Module, "SimpleDetangler", pybind11::base<Detangler>())
        .def(init<uint64_t, uint64_t, uint64_t>())
        ;
    class_<PermutationDetangler>(shasta2Module, "PermutationDetangler", pybind11::base<Detangler>())
        .def(init<uint64_t, double, double, double>())
        ;
    class_<ChiSquareDetangler>(shasta2Module, "ChiSquareDetangler", pybind11::base<Detangler>())
        .def(init<uint64_t, double, double, double>())
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
    shasta2Module.def("testSpoa",
        testSpoa
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
    shasta2Module.def("globalMsaPython",
        globalMsaPython
        );
}

#endif
