#!/usr/bin/python3

import shasta2

shasta2.openPerformanceLog("Python-performance.log")
options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.accessAnchorGraph()
assembler.createAssemblyGraph(options)
