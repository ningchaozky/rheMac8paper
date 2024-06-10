#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick
from ningchao.nBio import rheMac

example = '''hic file name rename'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('mat', nargs='?', help ='matrix')
parser.add_argument('-header', nargs='?', help ='header and row need to deal with', default = 0, type = int)
parser.add_argument('-row', nargs='+', help ='header and row need to deal with', default =  [ ], type = int)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

peirods = rheMac.rheMac().peirods()

dit = dict(zip(['20year','P0','4year','LV'], ['20Y','0M','45Y','liver']))
peirods.extend( dit.keys() )
fh = open( args.mat )


def deal( line_arr ):
    if isinstance( line_arr, str ):
        line_arr = [ line_arr ]
    for i,v in enumerate( line_arr ):
            for peirod in peirods:
                if peirod in v:
                    if peirod in dit:
                        peirod = dit[peirod]
                    if 'R1' in v:
                        line_arr[i] = '{}.rep1'.format( peirod)
                    elif 'R2' in v:
                        line_arr[i] = '{}.rep2'.format( peirod)
            line_arr[i] = rheMac.trick( line_arr[i] ).short()
    if len( line_arr ) == 1 :
        return line_arr[0]
    else :
        return line_arr


for j, line in enumerate( fh ):
    line_arr = line.rstrip().split('\t')
    if j == args.header :
        line_arr = deal( line_arr )
    if args.row:
        arr = trick.lst( line_arr ).get( args.row )
        out = deal( arr )
        for i in args.row:
            line_arr[ i ] = out[i]
    trick.write( line_arr )




















