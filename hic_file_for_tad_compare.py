#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, system, fix
from ningchao.nBio import chromosome

desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('hic', nargs='?', help = '.hic file')
parser.add_argument('-res', nargs='?', help = 'resolution', default = 10000, type = int)
parser.add_argument('-bed', nargs='?', help = 'bed file for split')
parser.add_argument('-d', nargs='?', help = 'outputdir')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


print ( 'Deal with file: ', args.hic, file = sys.stderr )
cmd = '''straw NONE {} {} {} BP {} > {}/{}'''
csize = chromosome.chr('rh8').size()
period = system.dir.str_map_period( args.hic )
outputdir = system.dir(args.d if args.d else period ).check()
for chrom in csize :
    print ( cmd.format( args.hic, chrom.replace('chr',''), chrom.replace('chr',''), args.res, outputdir, fix.fix(period).append( chrom, args.res, 'txt')))

if args.bed :
    cmd = 'tad_split_chroms.py {} {}'.format(args.bed, outputdir)
    print ( cmd )



























