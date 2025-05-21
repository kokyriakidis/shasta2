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
assemblyGraph3 = assembler.getAssemblyGraph3("B", assemblerOptions);

detangler = TrivialDetangler(assemblyGraphOptions.minCommonCoverage)	
assemblyGraph3.detangleEdges(detangler)
assemblyGraph3.compress()
assemblyGraph3.write("C")


