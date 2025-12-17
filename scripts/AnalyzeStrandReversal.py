#!/usr/bin/python3

# Get the argument.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('readId', type=int)
arguments = parser.parse_args()

import shasta2
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.analyzeStrandReversal(arguments.readId, True)

