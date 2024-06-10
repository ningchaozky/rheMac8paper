#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s xxoo.tab' % os.path.basename(sys.argv[0]), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'dir', nargs='?', help ='tab files in this dir' )
parser.add_argument( '-b', nargs='?', help = 'bed file before get the signal', default = '/home/ningch/data/genome/rheMac8/snp/snp.bed' )
parser.add_argument( '-sc', nargs='?', help = 'signal cut off', default = 3 , type = int )
parser.add_argument( '-db', action='store_true', help = 'debug the signal less than cutoff' )
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()
tabs =[ i for i in os.listdir(args.dir) if i.endswith('tab') ]
bpos = {}
bfh = open(args.b)
signal_cutoff = args.sc
for line in bfh:
    line_arr = line.strip().split('\t')
    trick.set1dict(bpos,line_arr[3],line_arr)


for tab in tabs:
    snp = tab.replace('.tab','')
    tfh = open(tab)
    tinfor = bpos[snp]
    header = tfh.next().strip().split('\t')[3:]
    lst = []
    for line in tfh:
        line_arr = line.strip().split('\t')
        chrom = line_arr.pop(0)
        start = line_arr.pop(0)
        end = line_arr.pop(0)
        for i,each in enumerate(line_arr):
            if each == 'NA':
                line_arr[i] = 0
            else :
                line_arr[i] = float(line_arr[i])
        hic_signal = dict(list(zip(header,line_arr)))
        for each in hic_signal:
            esignal = hic_signal[each]
            if esignal >= signal_cutoff :
                lst.append(int(start))
                lst.append(int(end))
    if len(lst) >= 1:
        rstart = min(lst)
        rend = max(lst)
        tinfor.insert(0,rend)
        tinfor.insert(0,rstart)
        tinfor.insert(0,chrom)
        print('\t'.join([ str(i) for i in tinfor ] ))
    elif args.db:
        sys.stdout.write(snp +': have no signal big than %d \n' % signal_cutoff )























