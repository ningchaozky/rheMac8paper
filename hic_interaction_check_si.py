#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s ' % os.path.basename(sys.argv[0]), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help ='dir for the pdf|tab files', default = '.')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

fls = os.listdir(args.dir)
def check(fl):
    fh = open(fl)
    next(fh)
    for line in fh:
        line_arr = line.strip().split('\t')
        if 'NA' in line:
            continue
        val = [ float(i) for i in line_arr[3:] if float(i) >= 3 ]
        if val :
            return fl
    return ''
    
for each in fls:
    prefix = ''
    if each.endswith('tab'):
        prefix = each.replace('.tab','')
    if os.path.exists(prefix + '.pdf'):
        tab = check(each)
        if tab:
            print('mv %s ./signal/' % tab)
            print('mv %s ./signal/' % tab.replace('.tab','.pdf'))

    



























