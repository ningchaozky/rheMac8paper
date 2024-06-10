#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick
from collections import defaultdict
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('bedpes', nargs='+', help = 'hicup_dedup.sort.bedpe bedpe file')
parser.add_argument('p', nargs='?', help = 'peirod E50 E90')
parser.add_argument('t', nargs='?', help = 'rep r1 r2 combine merge')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


outfiles = defaultdict( str )

for bedpe in args.bedpes :
    with open( bedpe ) as f :
        for line in f :
            chrom = line.rstrip().split('\t')[0]
            if chrom not in outfiles :
                outfiles[chrom] = open( '{}.bedpe.raw'.format( chrom), 'w')
            print ( line.rstrip(), file = outfiles[chrom] )
for chrom in outfiles:
    outfiles[chrom].close()
    outputfilename = 'cis.{}.PFC.{}.{}.hicup.dedup.sort.bedpe'.format(args.p, args.t, chrom)
    cmd = 'sort -k1,1 -k2,2n {} > {}'.format( outfiles[chrom].name, outputfilename)
    output = system.run( cmd, shell = True )
    for line in output :
        print ( line )




























