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


assemblyGraph = assembler.getAssemblyGraph("A", options)
detangler = LikelihoodRatioDetangler(
    options.detangleMinCommonCoverage,
    options.detangleEpsilon,
    options.detangleMaxLogP,
    options.detangleMinLogPDelta,
    options.detangleHighCoverageThreshold)
assemblyGraph.detangle(10, 1000000000, detangler);
assemblyGraph.write("X")


"""
assemblyGraph = assembler.getAssemblyGraph("E", options)
detangler = SimpleDetangler(1, 2)

edgeLengthThreshold = 1000

detangledVertexCount = assemblyGraph.detangleVerticesIteration(detangler);
print(detangledVertexCount, "successful vertex detagling operations.")
assemblyGraph.compress()
assemblyGraph.write("V1")

detangledVertexCount = assemblyGraph.detangleVerticesIteration(detangler);
print(detangledVertexCount, "successful vertex detangling operations.")
assemblyGraph.compress()
assemblyGraph.write("V2")


detangledEdgeCount = assemblyGraph.detangleEdgesIteration(edgeLengthThreshold, detangler);
print(detangledEdgeCount, "successful edge detagling operations.")
assemblyGraph.compress()
assemblyGraph.write("E1")

detangledEdgeCount = assemblyGraph.detangleEdgesIteration(edgeLengthThreshold, detangler);
print(detangledEdgeCount, "successful edge detagling operations.")
assemblyGraph.compress()
assemblyGraph.write("E2")

detangledEdgeCount = assemblyGraph.detangleEdgesIteration(edgeLengthThreshold, detangler);
print(detangledEdgeCount, "successful edge detagling operations.")
assemblyGraph.compress()
assemblyGraph.write("E3")

detangledEdgeCount = assemblyGraph.detangleEdgesIteration(edgeLengthThreshold, detangler);
print(detangledEdgeCount, "successful edge detagling operations.")
assemblyGraph.compress()
assemblyGraph.write("E4")
"""

"""
assemblyGraph = assembler.getAssemblyGraph("D", options)
assemblyGraph.search()
"""

"""
assemblyGraph = assembler.getAssemblyGraph("Z", options)
detangler = LikelihoodRatioDetangler(
    options.detangleMinCommonCoverage,
    options.detangleEpsilon,
    options.detangleMaxLogP,
    options.detangleMinLogPDelta)
detangler.debug = True
detangledEdgeCount = assemblyGraph.detangleEdgesIteration(1000000000, detangler);
"""

"""
assemblyGraph = assembler.getAssemblyGraph("B", options)
assemblyGraph.computeJourneys()
"""
