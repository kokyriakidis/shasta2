#!/usr/bin/python3

from shasta2 import *

assemblerOptions = AssemblerOptions()
assemblyGraphOptions = assemblerOptions.assemblyGraphOptions

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


# Get the AssemblyGraph2 at a chosen assembly stage and do something with it.
assemblyGraph2 = assembler.getAssemblyGraph2("B", assemblerOptions);
detangler = TrivialDetangler(assemblyGraphOptions.minCommonCoverage)
	
assemblyGraph2.detangleVertices(detangler)
assemblyGraph2.write("C")


