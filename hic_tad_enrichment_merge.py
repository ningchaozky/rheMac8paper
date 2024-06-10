#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
from ningchao.nSys import trick, system
from collections import defaultdict
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('tabs', nargs='+', help = 'tabs for merge')
parser.add_argument('-o', nargs='?', help = 'output file', default = sys.stdout)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


infor = defaultdict( lambda : defaultdict( float ) )
periods = system.dir.periods( prefrontal_raw = True )
for tab in args.tabs:
    with open( tab ) as f :
        header = f.readline().rstrip().split('\t')
        p1,p2 = [ system.dir.str_map_period(i) for i in header[1:3] ]
        for line in f :
            line_arr = line.rstrip().split('\t')
            gene = line_arr[0].split('.')[0]
            tdit = dict(zip([p1,p2], line_arr[1:3]))
            for p in periods:
                if p in tdit :
                    infor[p][gene] = tdit[p]
                else :
                    infor[p][gene] = 0
df = pd.DataFrame(infor)
df.to_csv( args.o, sep = '\t', index_label = 'symbol')




























