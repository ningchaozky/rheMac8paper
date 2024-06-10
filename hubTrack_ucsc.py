#!/usr/bin/env python
import sys
import copy
import argparse
import ningchao.usage as uTookit
import ningchao.nBio.hub as hubTookit

parser = argparse.ArgumentParser(prog = sys.argv[0],description='generate trackHub.txt from patter file')
parser.add_argument('-k', nargs='+', help ='key words for filter file split by "."')
parser.add_argument('-g', nargs='?', default = 'rheMac3', help ='genome for bw')
parser.add_argument('-p', nargs='?', help ='file prefix should different with other track')
parser.add_argument('-w', nargs='?', default = '.',help ='which directory you want to find bw files')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

pattern = args.k
genome = args.g
prefix = '.'.join(pattern)
if 'p' in args:
	prefix = args.p
marker = copy.copy(pattern)
work_dir= args.w
print('\n'.join(hubTookit.hub(genome,prefix,work_dir).db(pattern,marker)))
