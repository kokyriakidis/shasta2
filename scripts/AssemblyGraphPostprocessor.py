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
assemblyGraph3 = assembler.getAssemblyGraph3("Z", assemblerOptions)
# assemblyGraph3.analyzeSuperbubbles(10)
assemblyGraph3.colorStrongComponents()


"""
detangler = ChiSquareDetangler(
    assemblyGraphOptions.detangleMinCommonCoverage,
    assemblyGraphOptions.detangleEpsilon,
    assemblyGraphOptions.detangleMaxLogP,
    assemblyGraphOptions.detangleMinLogPDelta)
detangler.debug = True  
assemblyGraph3.detangleVertices(detangler)
assemblyGraph3.compress()
# assemblyGraph3.write("X")
"""
