#!/usr/bin/env python3
import os
import re
import sys
import argparse
import subprocess as sp
import multiprocessing as mp
from ningchao.nSys import system
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('cmd', nargs='?', help ='excel you want to mini', default = sys.stdin, type = argparse.FileType('r'))
parser.add_argument( '-c', nargs = '?', help = 'cpu number', type = int, default = 3)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def srun( cmd, check = '||'):
    f = system.run( cmd, shell = True, noReturn = True )
    for line in f:
        print ( line )
def multi_run( cmds ):
    with mp.Pool( args.c ) as p:
        p.map( srun, cmds )
    p.close()
    p.join()
cmds = re.split(r'\n{3,}', ''.join(args.cmd.readlines() ), re.M )
multi_run( cmds )





















