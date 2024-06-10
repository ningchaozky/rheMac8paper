#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import argparse
from tempfile import NamedTemporaryFile
from ningchao.nSys import trick, system
from ningchao.nBio import rheMac
example = '''hic_interaction_gene_pick_region.py results pick peak'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'interact_region', nargs='?', help ='output from hic_interaction_gene_pick_region.py' )
parser.add_argument( 'peak', nargs='?', help = 'peak bed file' )
parser.add_argument( '-p', nargs='?', help = 'peirod use for region extract' )
parser.add_argument( '-otype', choices = ['raw','peak'], nargs='?', help = 'peirod use for region extract', default = 'peak' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    peak, peirod, iregion, otype = args.peak, args.p, args.interact_region, args.otype
    if not peirod:
        peirods = rheMac.rheMac().peirods_combine_reps()
        peirods.extend(rheMac.rheMac().peirods())
        for key in peirods:
            if key in peak:
                peirod = key
                break
    return peak, iregion, peirod, otype

def iregion_infor( iregion, peirod ):
    ofh = NamedTemporaryFile( prefix='hic_interaction', suffix='.bed', dir='/tmp', delete = False )
    fh = open( iregion )
    header = fh.next().strip().split('\t')
    keys = [ i for i in header if peirod in i ]
    for line in fh:
        tmp = dict( list(zip( header, line.strip().split('\t') )) )
        gene = tmp['gene']
        for each in trick.dit(tmp).get(keys):
            if not each.strip():
                continue
            each_dit = eval(each)
            for egene in each_dit :
                bed_string = each_dit[egene]
                for each in bed_string.split(','):
                    ofh.write( each + '\t' + gene +'.'+egene +'\n')
    ofh.close()
    return ofh.name
def pbed( peak ):
    fh,ofh = open(peak), NamedTemporaryFile( prefix='hic_interaction_peak', suffix='.bed', dir='/tmp', delete = False )
    for line in fh:
        if 'track' in line:
            continue
        line_arr = line.strip().split('\t')
        ofh.write( '\t'.join(line_arr[0:3]) +'\n' )
    ofh.close()
    return ofh.name

def write( line, otype ):
    if not line.strip():
        return ''
    line_arr = line.strip().split('\t')
    if otype == 'peak':
        index = [4, 5, 6 , 3]
        line_arr = trick.lst(line_arr).get(index)
    return '\t'.join( line_arr )


if __name__ == '__main__':
    peak, iregion, peirod, otype = parse( args )
    iregion_bed = iregion_infor( iregion, peirod)
    peak_bed = pbed( peak )
    cmd = 'bedtools intersect -b %s -a %s -wo' % (peak_bed, iregion_bed)
    sys.stderr.write( cmd + '\n' )
    stdout,stderr,returncode = system.run( cmd )
    for line in stdout:
        if line.strip():
            print(write( line, otype ))























