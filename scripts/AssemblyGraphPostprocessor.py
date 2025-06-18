#!/usr/bin/python3

from shasta2 import *

openPerformanceLog("Python-performance.log")
assemblerOptions = AssemblerOptions()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()

# Get the assemblyGraph3Postprocessor at a chosen assembly stage.
assemblyGraph = assembler.getAssemblyGraph("C", assemblerOptions)
assemblyGraph.phaseSuperbubbleChains(10)
# assemblyGraph.colorStrongComponents()


"""
detangler = ChiSquareDetangler(
    assemblerOptions.detangleMinCommonCoverage,
    assemblerOptions.detangleEpsilon,
    assemblerOptions.detangleMaxLogP,
    assemblerOptions.detangleMinLogPDelta)
detangler.debug = True  
assemblyGraph.detangleVertices(detangler)
assemblyGraph.compress()
# assemblyGraph.write("X")
"""
