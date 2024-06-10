#!/usr/bin/env python3
import os
import sys
import pyBigWig
import argparse
from collections import defaultdict
from ningchao.nSys import trick, fix, system
desc = '''Convert the chromosome names of a bigWig file.'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("input", metavar="input.bigWig", help="Input bigWig file")
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

# read in the name map
bw = pyBigWig.open(args.input)
# Make a new header
if 'chr' in list(bw.chroms().keys())[0]:
    print ('chr already in {}'.format(args.input), file = sys.stderr)
    exit()
cwd = os.path.abspath( args.input )
output = fix.fix(cwd).insert('chrInsert')
hdr = [ ('chr'+chrom, length) for chrom, length in bw.chroms().items() if 'chr' not in chrom ]

bwOutput = pyBigWig.open( output, "w")
bwOutput.addHeader(hdr)
for chrom, length in bw.chroms().items():
    ints = bw.intervals(chrom, 0, length)
    if len(ints):
        bwOutput.addEntries(['chr' + chrom] * len(ints),
                            [x[0] for x in ints],
                            ends=[x[1] for x in ints],
                            values=[x[2] for x in ints])
bw.close()
bwOutput.close()
system.run('mv {} {}'.format( output, args.input), shell = True, noReturn = True )





