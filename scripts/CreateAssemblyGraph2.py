#!/usr/bin/python3

import shasta2

shasta2.openPerformanceLog("Python-performance.log")
assemblerOptions = shasta2.AssemblerOptions()
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.createAssemblyGraph2(assemblerOptions, assemblerOptions.threadCount)
