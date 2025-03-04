#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "Base.hpp"
#include "ConfigurationTable.hpp"
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



PYBIND11_MODULE(shasta, shastaModule)
{

    // Expose class Assembler to Python.
    class_<Assembler>(shastaModule, "Assembler")

        // Constructor.
        .def(pybind11::init<const string&, bool, size_t>(),
            "Assembler constructor.",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("createNew") = false,
            arg("largeDataPageSize") = 2*1024*1024)

        // Reads
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")

        // K-mers.
        .def("accessKmerChecker",
            &Assembler::accessKmerChecker)

         // Markers.
        .def("accessMarkers",
            &Assembler::accessMarkers)
        .def("createMarkers",
            &Assembler::createMarkers,
            "Find markers in reads.",
            arg("threadCount") = 0)
        .def("createMarkerKmers",
            &Assembler::createMarkerKmers,
            arg("threadCount") = 0)
        .def("accessMarkerKmers",
            &Assembler::accessMarkerKmers)
    ;



    // Expose portions of class AssemblerOptions to Python.
    class_<AssemblerOptions>(shastaModule, "AssemblerOptions")
        .def(pybind11::init<const string&>())
        .def_readonly("assemblyOptions", &AssemblerOptions::assemblyOptions)
        ;
    class_<AssemblyOptions>(shastaModule, "AssemblyOptions")
        .def_readonly("mode3Options", &AssemblyOptions::mode3Options)
        ;
    class_<Mode3AssemblyOptions>(shastaModule, "Mode3AssemblyOptions");



    // Non-member functions exposed to Python.
    shastaModule.def("openPerformanceLog",
        openPerformanceLog
        );
    shastaModule.def("testMultithreadedObject",
        testMultithreadedObject
        );
    shastaModule.def("testMemoryMappedVector",
        testMemoryMappedVector
        );
    shastaModule.def("testBase",
        testBase
        );
    shastaModule.def("testShortBaseSequence",
        testShortBaseSequence
        );
    shastaModule.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shastaModule.def("testSplitRange",
        testSplitRange
        );
    shastaModule.def("testSpoa",
        testSpoa
        );
    shastaModule.def("testDeduplicateAndCount",
        testDeduplicateAndCount
        );
    shastaModule.def("testBitReversal",
        testBitReversal
        );
    shastaModule.def("mappedCopy",
        mappedCopy
        );
    shastaModule.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shastaModule.def("testDiploidBayesianPhase",
        testDiploidBayesianPhase
        );
    shastaModule.def("testSubsetGraph",
        testSubsetGraph
        );
    shastaModule.def("globalMsaPython",
        globalMsaPython
        );
}

#endif
