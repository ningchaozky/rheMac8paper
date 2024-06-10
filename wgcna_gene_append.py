#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,system

example = '''matrx'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'matrix', nargs='?', help ='matrix want to append' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ) :
    return args.matrix

def fl2color():
    fls, dit = system.dir('.').fls(level = 0, pattern = r'%s.*txt' % matrix, fpattern = 'Cytoscape|go' ),{}
    for key in fls:
        fl = key
        key = key.split('.')[-2]
        fh = open(fl)
        next(fh)
        for line in fh:
            line_arr = line.strip().split('\t')
            trick.dinit( dit, line_arr[-1], [] ) 
            dit[ line_arr[-1] ].append( key )
    return dit
if __name__ == '__main__':
    matrix = parse( args )
    dit = fl2color()
    fh = open( matrix )
    print(fh.next().strip() + '\t' + 'cluster')
    for line in fh:
        line_arr = line.strip().split('\t')
        gene = line_arr[0]
        if gene in dit:
            line_arr.append( '.'.join( dit[gene] ) )
        else :
            line_arr.append( 'None' )
        print('\t'.join( line_arr ))


























