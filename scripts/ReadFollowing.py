#!/usr/bin/python3

from shasta2 import *


# Get the arguments.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('anchorId0', type=str)
parser.add_argument('direction', type=int)
parser.add_argument('minCommon', type=int)
arguments = parser.parse_args()

openPerformanceLog("Python-performance.log")
options = Options()
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


anchorId0 = "100+"
anchorId1 = assembler.readFollowing(anchorIdFromString(arguments.anchorId0), arguments.direction, arguments.minCommon)
print(anchorIdToString(anchorId1));
