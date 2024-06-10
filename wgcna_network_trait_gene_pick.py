#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick

example = '''matrix'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'matrix', nargs= '?', help ='matrix for trait gene pick' )
parser.add_argument( '-peirod','-p', nargs= '?', help = 'peirod you want to pick', required = True )
parser.add_argument( '-maxCut', nargs= '?', help = 'max cut off for the hard pick', default = 1, type = float )
parser.add_argument( '-nfc', nargs= '?', help = 'neighbor fold change cut', default = 1.2, type = float )
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

matrix = args.matrix
mfh = open(matrix)
peirod = args.peirod
header = mfh.next().strip()
print(header)
index = [ i for i,v in enumerate(header.split('\t')) if peirod in v ][0]

def com(bval, val, aval):
    nfc = args.nfc
    if not aval :
        aval = 1
    if not bval:
        bval = 1
    if aval == 1 and val / bval > nfc:
        return True 
    elif bval == 1 and val / aval > nfc:
        return True
    elif aval != 1 and bval != 1 :
        if  val / aval > nfc or val / bval > nfc :
            return True 
    else :
        return False 

def check( line_arr, index, maxCut) :
    try :
        bval = float ( line_arr[ index - 1 ] )
    except :
        bval = 1
    try :
        aval = float ( line_arr[ index + 1 ] )
    except :
        aval = 1
    val = float ( line_arr[index] )
    arr = [ float(i) for i in line_arr[1:] ]
    amax = max(arr)
    if val == amax and amax >= maxCut and com( bval, val, aval ):
    #if amax >= 0.0001 and com( bval, val, aval ):
        return True
    else :
        return False 

for line in mfh:
    line_arr = line.strip().split('\t')
    if check( line_arr, index, args.maxCut) :
        print(line, end=' ')



























