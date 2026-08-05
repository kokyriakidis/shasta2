#!/usr/bin/python3

import shasta2


shasta2.openPerformanceLog("Python-performance.log")

# Get the options from shasta2.conf.
options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

# Load the specified assembly stage and assemble sequence.
assemblyGraph = assembler.getAssemblyGraph("A", options)
assemblyGraph.detangleVertices()
assemblyGraph.write("A-Detangled")





