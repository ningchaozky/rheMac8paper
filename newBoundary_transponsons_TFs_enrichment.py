#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
from collections import defaultdict
from ningchao.nSys import trick, system, fix
from ningchao.nBio import transposon
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'bed', nargs='?', help = 'boundary for find transposons and do enrichment')
parser.add_argument( '-t', nargs='?', help = 'transposon file', default = '')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def get_rsmk_bed():
    intron = '/home/soft/data/genome/rheMac8/intron_uniq.bed'
    rsmk_beds = transposon.transponson( 'rheMac8' ).no_intron_transposon( intron, level = [1] )
    return rsmk_beds
def get_new_boundary_rsmk( nb, beds ):
    dit = defaultdict( str )
    for bed in beds:
        cmd = 'bedtools intersect -a {} -b {}'.format( nb, beds[bed] )
        nb_noIntron = open( fix.fix(nb).insert( bed, 'noIntron'), 'w' )
        for line in system.run( cmd, shell = True):
            print ( line.strip(), file = nb_noIntron )
        nb_noIntron.close()
        dit[bed] = nb_noIntron.name
    return dit
def do_enrichment( beds ):
    cmd = 'bed_motif_enrichment.py {}'
    for bed in beds:
        if 'SVA' not in bed:
            continue
        for line in system.run( cmd.format( beds[bed] ), shell = True ):
            os.system( line )


if __name__ == '__main__':
    rsmk_beds = get_rsmk_bed()
    nb_noIntron = get_new_boundary_rsmk( args.bed, rsmk_beds)
    do_enrichment( nb_noIntron )
























