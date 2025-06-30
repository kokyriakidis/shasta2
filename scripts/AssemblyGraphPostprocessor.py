#!/usr/bin/python3

from shasta2 import *

openPerformanceLog("Python-performance.log")
options = Options()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()

"""
assemblyGraph = assembler.getAssemblyGraph("C", options)
assemblyGraph.phaseSuperbubbleChains()
"""

assemblyGraph = assembler.getAssemblyGraph("D", options)
detangler = LikelihoodRatioDetangler(
    options.detangleMinCommonCoverage,
    options.detangleEpsilon,
    options.detangleMaxLogP,
    options.detangleMinLogPDelta)
assemblyGraph.detangle(10, 1000000000, detangler);
assemblyGraph.write("X")
