#!/usr/bin/env python3
import os
import sys
import argparse
import pyBigWig
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'bw', nargs = '?', help = 'bw file')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



bw = pyBigWig.open(args.bw)
print ( bw.chroms() )
bw.close()




















