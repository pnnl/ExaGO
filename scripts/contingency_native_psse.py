#!/usr/bin/env python
# -*- mode: python -*-
# -------------------------------------------------------------
# file: contingency_format.py
#
# A hack script to convert ExaGO NATIVE contingencies to PSS/E
# -------------------------------------------------------------

import sys
import os
from optparse import OptionParser
import csv
from dataclasses import dataclass
import enum

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])

# -------------------------------------------------------------
# Outage
# -------------------------------------------------------------


class Outage(enum.Enum):
    GEN = 1
    BR = 2
    TR = 3
    LOAD = 4

# -------------------------------------------------------------
# Contingency
# -------------------------------------------------------------


@dataclass
class Contingency:
    num: int
    type: Outage
    bus: int
    fbus: int
    tbus: int
    id: str
    status: int
    prob: float

    def WritePSSE(self, f):
        f.write("CONTINGENCY ")
        if self.type == Outage.GEN:
            self.WritePSSE_GEN(f)
        elif self.type == Outage.BR:
            self.WritePSSE_BR(f)
        else:
            raise ValueError("Unsupported outage type: %s" % self.type)
        f.write("END\n")

    def WritePSSE_GEN(self, f):
        unit = 1
        name = ("G_%d_%d_%d" % (self.num, unit, self.bus))
        f.write("%s\nREMOVE UNIT %d FROM BUS %d\n" %
                (name, unit, self.bus))

    def WritePSSE_BR(self, f):
        circuit = 1
        name = ("L_%d-%d_%d_%d" %
                (self.fbus, self.tbus, circuit, self.num))
        f.write("%s\nOPEN BRANCH FROM BUS %d TO BUS %d CIRCUIT %d\n" %
                (name, self.fbus, self.tbus, circuit))


# -------------------------------------------------------------
# handle command line
# -------------------------------------------------------------
usage = "Usage: %prog [options] [file.cont]"
parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose",
                  dest="verbose", action="store_true", default=False,
                  help="show what's going on")

parser.add_option("-o", "--output", type="string",
                  dest="output", action="store")

(options, args) = parser.parse_args()

doverbose = options.verbose

if options.output:
    output = open(options.output, mode='w')
else:
    output = sys.stdout


if (len(args) > 0):
    inputfile = args[0]
    fd = open(inputfile, newline='')
    if fd is None:
        sys.stderr.write("%s: %s: error: cannot open\n" % (program, inputfile))
        exit(2)
else:
    fd = sys.stdin

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

reader = csv.reader(fd)
l = 0
for row in reader:
    l = l + 1
    c = Contingency(
        int(row[0]),
        Outage(int(row[1])),
        int(row[2]),
        int(row[3]),
        int(row[4]),
        row[5],
        int(row[6]),
        float(row[7]))
    c.WritePSSE(output)
