#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
from collections import defaultdict
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='?', help ='tab file to deal' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



with open( args.tab ) as f:
    header, infor = next( f ).strip(), defaultdict( list )
    print ( header )
    for line in f:
        line_arr = line.strip().split('\t')
        symbol = line_arr[0].split('.')[0]
        infor[symbol].append( line_arr[1:] )
    for symbol in infor :
        val = np.array( infor[symbol] ).astype(float).sum( axis = 0 )
        print ( symbol, *val, sep = '\t' )



































