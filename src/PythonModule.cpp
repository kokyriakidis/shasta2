#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
#include "diploidBayesianPhase.hpp"
#include "globalMsa.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultithreadedObject.hpp"
#include "performanceLog.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
#include "testSpoa.hpp"
#include "testSubsetGraph.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta2, shasta2Module)
{
    // Expose class AssemblerOptions to Python.
    class_<AssemblerOptions>(shasta2Module, "AssemblerOptions")

        // Constructor from the name of a configuration file.
        .def(pybind11::init<const string&>(),
            arg("configurationFileName") = "shasta2.conf")
        .def_readonly("k", &AssemblerOptions::k)
        .def_readonly("markerDensity", &AssemblerOptions::markerDensity)
        ;



    // Expose class Assembler to Python.
    class_<Assembler>(shasta2Module, "Assembler")

        // Constructor.
        .def(pybind11::init<const string&>(),
            "Assembler constructor.",
            arg("largeDataFileNamePrefix") = "Data/")

        // Reads
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")

        // Marker K-mers.
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
        .def("createMarkerKmers",
            &Assembler::createMarkerKmers,
            arg("threadCount") = 0)
        .def("accessMarkerKmers",
            &Assembler::accessMarkerKmers)
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
    shasta2Module.def("globalMsaPython",
        globalMsaPython
        );
}

#endif
