#!/usr/bin/env python
# -------------------------------------------------------------
# file: contingency_psse_pti.py
#
# A hack script to convert PSS/E contingencies to ExaGO NATIVE 
# -------------------------------------------------------------
# -------------------------------------------------------------

import sys, os
from optparse import OptionParser
from pyparsing import *
import enum

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])

class Outage(enum.Enum):
    GEN = 1
    BR = 2
    TR = 3
    LOAD = 4

Name = Word(printables)

intid = Word(nums)

From = Keyword("FROM", caseless=True)
To = Keyword("TO", caseless=True)
Open = Keyword("OPEN", caseless=True)

Bus = Keyword("BUS", caseless=True)
Branch = Keyword("BRANCH", caseless=True)
Circuit = Keyword("CIRCUIT", caseless=True)
Remove = Keyword("REMOVE", caseless=True)
Unit = Keyword("UNIT", caseless=True)

LineOutage = (
    Open.suppress() + Branch.suppress() + 
    From.suppress() + Bus.suppress() + intid.setName("from_bus") +
    To.suppress() + Bus.suppress() + intid.setName("to_bus") +
    Circuit.suppress() + intid.setName("circuit")
).setParseAction(lambda t: (Outage.BR, t[0:]) ) 

GenOutage = (
    Remove.suppress() + Unit.suppress() +
    intid.setName("gen") +
    From.suppress() + Bus.suppress() + intid.setName("bus")
).setParseAction(lambda t: (Outage.GEN, t[0:]) ) 

Contingency = (
    Keyword("CONTINGENCY", caseless=True).suppress() +
    Name +
    OneOrMore( LineOutage | GenOutage ) +
    Keyword("END", caseless=True).suppress()
).setParseAction(lambda t: { t[0] : t[1:] })

ContingencyList = OneOrMore(Contingency)

# -------------------------------------------------------------
# handle command line
# -------------------------------------------------------------
usage = "Usage: %prog [options] file.con"
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

conts = ContingencyList.parseFile(fd)

id = 0
for (c) in conts:
    for n in c:
        id = id + 1
        for o in c[n]:
            output.write("%d," % (id))
            (typ, dat) = o
            if typ == Outage.BR:
                output.write("%d,0,%d,%d," %
                             (int(typ.value), int(dat[0]), int(dat[1])))
            elif typ == Outage.GEN:
                output.write("%d,%d,0,0," %
                             (int(typ.value), int(dat[1])))
            output.write("'1 ',1,0,0.01\n")
                         
