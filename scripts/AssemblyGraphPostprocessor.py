#!/usr/bin/python3

from shasta2 import *

openPerformanceLog("Python-performance.log")
assemblerOptions = AssemblerOptions()
assemblyGraphOptions = assemblerOptions.assemblyGraphOptions

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()

# Get the assemblyGraph3Postprocessor at a chosen assembly stage.
assemblyGraph = assembler.getAssemblyGraph("Z", assemblerOptions)
# assemblyGraph.analyzeSuperbubbles(10)
assemblyGraph.colorStrongComponents()


"""
detangler = ChiSquareDetangler(
    assemblyGraphOptions.detangleMinCommonCoverage,
    assemblyGraphOptions.detangleEpsilon,
    assemblyGraphOptions.detangleMaxLogP,
    assemblyGraphOptions.detangleMinLogPDelta)
detangler.debug = True  
assemblyGraph.detangleVertices(detangler)
assemblyGraph.compress()
# assemblyGraph.write("X")
"""
