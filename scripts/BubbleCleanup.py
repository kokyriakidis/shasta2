#!/usr/bin/python3

from shasta2 import *

openPerformanceLog("Python-performance.log")
options = Options()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessAnchors()
assembler.accessJourneys()


assemblyGraph = assembler.getAssemblyGraph("A", options)
assemblyGraph.bubblePairCleanupIterationMultithreaded([])



