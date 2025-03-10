#!/usr/bin/python3

import shasta2

assemblerOptions = shasta2.AssemblerOptions()
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessMarkerKmers()
assembler.createAnchors(assemblerOptions.minAnchorCoverage, assemblerOptions.maxAnchorCoverage)

