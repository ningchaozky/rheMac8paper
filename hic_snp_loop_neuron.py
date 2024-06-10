#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s ' % os.path.basename(sys.argv[0]), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'loop', nargs='?', help ='loop file' )
parser.add_argument( '-n', nargs='?', help ='neuron gene bed', default = '/home/ningch/data/genome/rheMac8/exon/ref_ensemble_xeonRefGene_neuron.bed')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()
loop = args.loop
nfile = args.n
nfh = open(nfile)
infor = {}
for line in nfh:
    line_arr = line.strip().split('\t')
    chrom, start, end, name, ph, chain = line_arr
    if chain == '+' :
        gbin = int(start) / 40000 + 1
    else :
        gbin = int(end) / 40000 + 1
    trick.set2dict(infor,chrom,gbin,[])
    infor[chrom][gbin].append(line)

nfh.close()

lfh = open(loop)
for line in lfh:
    line_arr = line.strip().split('\t')
    chrom = line_arr.pop(0)
    start = int(line_arr.pop(0))
    end = int(line_arr.pop(0))
    sbin = start / 40000
    ebin = end / 40000
    bins = list(range(sbin,ebin))
    bins.append(ebin)
    bins.insert(0,sbin - 1)
    lst = list(range(start,end))
    lst.append(end)
    lst.insert(0, start - 1 )
    for each in bins:
        if each in infor[chrom]:
            for gene in infor[chrom][each]:
                chrom, start, end, name, ph, chain = gene.strip().split('\t')
                if chain == '+' :
                    if int(start) in lst :
                        print(line.strip() + '\t' + gene, end=' ')
                else :
                    if int(end) in lst :
                        print(line.strip() + '\t' + gene, end=' ')




















