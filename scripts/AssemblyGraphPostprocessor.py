#!/usr/bin/python3

from shasta2 import *

assemblerOptions = AssemblerOptions()
assemblyGraphOptions = assemblerOptions.assemblyGraphOptions

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


# Get the AssemblyGraph3 at a chosen assembly.
assemblyGraph3 = assembler.getAssemblyGraph3("B", assemblerOptions);

# detangler = TrivialDetangler(assemblyGraphOptions.minCommonCoverage)	
# assemblyGraph2.detangleVertices(detangler)
assemblyGraph3.write("C")


