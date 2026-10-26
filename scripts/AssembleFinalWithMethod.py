#!/usr/bin/python3

# Reassembles the sequence of an assembly's Final AssemblyGraph with a chosen
# LocalAssembly7 Method, so whole-assembly accuracy can be compared between
# consensus methods on an identical graph (only the local assemblies differ).
# Run in the assembly directory. The Method and msaRepair come from the
# environment, as LocalAssembly7 reads them:
#   SHASTA2_LA7_METHOD   Adaptive (the default), Abpoa, TheseusOnly, TheseusAll or Plegma
#   SHASTA2_LA7_REPAIR   1 or 0: msaRepair on or off
#   SHASTA2_PLEGMA       the plegma executable, for Method Plegma
# Writes Assembly-Final-<name>.fasta.
#
# usage: AssembleFinalWithMethod.py <name> [homopolymer model name]

import shasta2
import sys

name = sys.argv[1]
homopolymerModelName = sys.argv[2] if len(sys.argv) > 2 else ""

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.createHomopolymerModel(homopolymerModelName)
assemblyGraph = assembler.getAssemblyGraph("Final", options)
assemblyGraph.clearSequence()
assemblyGraph.assembleAllStrandSymmetric()
assemblyGraph.writeFasta("Final-" + name)
