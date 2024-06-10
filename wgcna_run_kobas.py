#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,system
example = '''dir'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help ='dir for run kobas')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

fls = system.dir(args.dir).fls(pattern = 'txt$', fpattern = 'Cy', level = 0 )
for fl in fls:
    #fl = fl.strip('./')
    #print('run_kobas_ningch.py -i %s 2 -o %s.go -t gene -nt ncbigene ' % (fl,fl))
    print ('run_kobas_ningch_v2.py {},index_col,2'.format( fl ) )
    print('run_kobas.py -i {}.kobas.all.id -o {}.kobas.all.id.go -s hsa -t id:ncbigene'.format(fl,fl))



























