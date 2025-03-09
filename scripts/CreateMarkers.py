#!/usr/bin/python3

import shasta2

assemblerOptions = shasta2.AssemblerOptions()
assembler = shasta2.Assembler()
assembler.createKmerChecker(assemblerOptions.k, assemblerOptions.markerDensity)
assembler.createMarkers()

