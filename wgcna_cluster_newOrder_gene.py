#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, system

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'order', nargs = '?', help = 'reoder list' )
#parser.add_argument( '-p', nargs = '?', help = 'prefix for the wgcna output file', default = 'select' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



def parse( args ):
    #fls = list ( system.dir( os.getcwd() ).fls( '%s.*txt' % args.p, fpattern = 'Cytoscape' ) )
    fls = list ( system.dir( os.getcwd() ).fls( 'txt$', fpattern = 'Cytoscape' ) )
    fls = dict( list(zip( [ i.split('.')[-2] for i in fls ], fls )))
    return fls, args.order

def order( fls, order_fl ):
    colors = [ i.strip() for i in open(order_fl).readlines() ]
    for color in colors:
        color = color.replace('ME','')
        if color in fls:
            lines = [ i.strip().split('\t')[-1] for i in open(fls[color]).readlines() if 'modProbes' not in i ]
            print('\n'.join( lines ))



def main( args ) :
    fls, order_fl = parse( args )
    order( fls, order_fl )
if __name__ == '__main__':
    main( args )


























