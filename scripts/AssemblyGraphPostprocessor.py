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


# Get the assemblyGraph3Postprocessor at a chosen assembly.
assemblyGraph3 = assembler.getAssemblyGraph3("D", assemblerOptions)
assemblyGraph3.assembleAll(0)
assemblyGraph3.write("Z1");

"""
detangler = TrivialDetangler(assemblyGraphOptions.minCommonCoverage)

for iteration in range(3):
	assemblyGraph3.detangleVertices(detangler)
	assemblyGraph3.compress()
	
	assemblyGraph3.detangleEdges(detangler)
	assemblyGraph3.compress()

assemblyGraph3.write("C")
"""

