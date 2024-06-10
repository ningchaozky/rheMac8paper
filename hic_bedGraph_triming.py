#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick
from ningchao.nBio import chromosome

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='bedGraph file trimming', formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('bedGraph', nargs='?',help='bedGraph file for trimming')
parser.add_argument('-o', nargs='?', help ='reference')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

bedGraph = args.bedGraph
ofh = sys.stdout
if args.o:
    ofh = open(args.o,'w')


chroms = chromosome.chr('rh8').chr
def lst2string(lst):
    return '\t'.join([ str(i) for i in lst ]) + '\n'
def parse(line):
    line_arr = line.strip().split('\t')
    line_arr[1] = int(line_arr[1])
    line_arr[2] = int(line_arr[2])
    return line_arr
bfh = open(bedGraph)
for line in bfh:
    if line.startswith('track'):
        ofh.write(line)
        continue
    chrom,start,end,val = parse(line)
    chrom_len = chroms[chrom]
    if start <= chrom_len <= end:
        end = chrom_len
    ofh.write(lst2string([chrom,start,end,val]))
    


















