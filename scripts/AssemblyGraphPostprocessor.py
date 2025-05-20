#!/usr/bin/python3

from shasta2 import *

assemblerOptions = AssemblerOptions()
assemblyGraphOptions = assemblerOptions.assemblyGraphOptions

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


# Get the assemblyGraph3Postprocessor at a chosen assembly.
assemblyGraph3Postprocessor = assembler.getAssemblyGraph3("B", assemblerOptions);

detangler = TrivialDetangler(assemblyGraphOptions.minCommonCoverage)	
assemblyGraph3Postprocessor.detangleEdges(detangler)
# assemblyGraph3Postprocessor.write("C")


