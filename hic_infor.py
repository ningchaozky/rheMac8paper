#!/usr/bin/env python3
import os
import sys
import fanc
import argparse
from ningchao.nSys import trick
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'hic', nargs = '?', help = 'hic file')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()




hic = fanc.load( args.hic )
print ( hic.chromosomes() )
print ( ['VC', 'VC_SQRT', 'KR'] )
print ( [ 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000] )


























