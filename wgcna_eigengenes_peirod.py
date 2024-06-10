#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick,system

example = '''NA.yellow.txt -marker K4 K27 -peirod E50'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('-eigs','-e', nargs='*', help = 'eigengenes group' )
parser.add_argument('-marker','-m', nargs='*', help = 'markers' )
parser.add_argument('-peirod','-p', nargs='*', help = 'peirods' )
parser.add_argument('-reverse','-r', action='store_false', help = 'if not in print' )
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()
marker = ' '.join(args.marker)
peirod = ' '.join(args.peirod)
xls = '/dataB/ftp/pub/rheMac3/prefrontal/binary/meta/chip_mark_tss.xls'
cmd = 'binary_analysis_methods.py -xls %s -marker %s -peirod %s' % (xls, marker, peirod)
stdout,stderr,stopcode = system.run( cmd )
genes = []
for line in stdout:
    genes.append(line.strip())

def reverse(gene, genes, rev):
    if rev and gene in genes:
        print(gene)
    if not rev and gene not in genes:
        print(gene)

eigs = args.eigs
for eig in eigs:
    efh = open(eig)
    next(efh)
    rev = args.reverse
    for line in efh:
        line_arr = line.strip().split('\t')
        reverse( line_arr[1], genes, rev)
    efh.close()























