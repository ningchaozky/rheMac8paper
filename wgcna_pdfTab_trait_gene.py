#!/usr/bin/env python3

import os
import sys
import argparse
from ningchao.nSys import trick,num,system
example = ''' mini.m1.addNum.quantile_normmoduleEigengenes.tab '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='?', help = 'output tab form Eigengenes pdf' )
parser.add_argument( '-p', nargs='?', help = 'peirods you want to extract' )
parser.add_argument( '-c', nargs='?', help = 'cut off for anaysis', default = 0.2, type = float )
parser.add_argument( '-prefix', nargs='?', help = 'prefix for the files' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


tab,tfh,prefix,cut = args.tab, open(args.tab), args.prefix, args.c
header = tfh.next().rstrip().split('\t')
header[0], peirod = 'module', args.p
peirods = header[1:]
wd,tname = os.path.split(os.path.abspath(tab))
if not prefix :
    prefix = tname.replace('moduleEigengenes.tab','')
fls = list(system.dir( wd ).fls(pattern = prefix, fpattern = 'Cytoscape|png|tab|pdf|^%s$' % prefix))
keys = [ i.split('.')[-2] for i in fls ]
fls = dict(list(zip(keys,fls)))
if peirod not in peirods:
    exit('Should in %s' % str(peirods))
for line in tfh:
    line_arr = line.strip().split('\t')
    line_arr[1:] = [ num.num(i).round() for i in line_arr[1:] ]
    dit = dict(list(zip(header, line_arr)))
    if peirod in dit:
        if dit[peirod] >= cut :
            color = dit['module'].replace('ME','')
            cfh = open(fls[color])
            for line in cfh:
                line_arr = line.rstrip().split('\t')
                if len(line_arr) >= 2 :
                    print('\t'.join( [ line_arr[1], color ] ))
            cfh.close()




























