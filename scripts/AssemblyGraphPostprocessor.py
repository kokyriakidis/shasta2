#!/usr/bin/python3

from shasta2 import *

assemblerOptions = AssemblerOptions()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()

# Get the AssemblyGraph at a chosen assembly stage and do something with it.
assemblyGraph = assembler.getAssemblyGraph("B");
detangler = TrivialDetangler()
assemblyGraph.detangleVertices(detangler)