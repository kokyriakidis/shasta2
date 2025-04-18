#!/usr/bin/python3

from shasta2 import *

options = AssemblerOptions()
assemblyGraphOptions = options.assemblyGraphOptions

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


# Get the AssemblyGraph at a chosen assembly stage and do something with it.
assemblyGraph = assembler.getAssemblyGraph("Z");
detangler = PermutationDetangler(assemblyGraphOptions.minCommonCoverage)
assemblyGraph.detangleEdges(detangler)


"""
assemblyGraph = assembler.getAssemblyGraph("A");
assemblyGraph.transitiveReduction(
    assemblyGraphOptions.transitiveReductionThreshold,
    assemblyGraphOptions.transitiveReductionA,
    assemblyGraphOptions.transitiveReductionB);
"""