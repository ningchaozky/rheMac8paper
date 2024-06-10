#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick
from ningchao.nBio import chromosome,geneKit,bedKit

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='this script already concert the chain\nbed_peaks_mapping.py -bed /home/ningch/data/genome/rheMac8/exon/ref.ens.xeon.coding.gene -peak K4.bed.bw.tab', formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('-bed', nargs='?', help ='bed file you want to extract peaks', default = '/home/ningch/data/genome/rheMac8/exon/ref.ens.xeon.coding.gene')
parser.add_argument('-peak', nargs = '?', help ='peak file for the mapping', required = True )
parser.add_argument('-s', '-span', help = 'span for the mapping, default 500000', nargs = '?', type = int , default = 500000 )
parser.add_argument('-onlyMapping', '-om', help = 'only output mapping peak', action='store_false')
parser.add_argument('-ignore_promoter', action='store_true', help = 'ignore promoter ')
parser.add_argument('-o', nargs = '?', help =' output file' )
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

bed = open(args.bed)
peak = open(args.peak)
out = sys.stdout
span = args.s
pspan = 2500
debug = False
if args.o :
    out = open(args.o,'w')
onlyMapping = args.onlyMapping

ignore_promoter = args.ignore_promoter

def get_bin( pos, span):
    pos = float(pos)
    span = float(span)
    float_return = pos / span
    int_return = int(pos) / int(span)
    if float_return == int_return :
        return int_return
    else :
        return int_return + 1


def gene_in_span(region, gstart, span):
    lst = list(range(region[0], region[1]))
    lst.append( region[1] )
    for i in lst:
        if abs(i - gstart) == span :
            return True
    return False



infor = {}
bed_infor = {}
chroms = chromosome.chr('rh8').chr
gene_infor = {}
for line in bed:
    line_arr = line.strip().split('\t')
    chrom = line_arr[0]
    if chrom not in chroms:
        continue
    if len(line_arr) >= 5 :
        chain = line_arr[5]
        gene = geneKit.name(line_arr[3])
    else :
        print('should be the bed6')
    start = int(line_arr[1])
    end = int(line_arr[2])
    if chain == '+':
        pos_bin = get_bin(start, span)
    else :   
        pos_bin = get_bin(end, span)

    trick.set2dict( infor, chrom, pos_bin, [])
    key = ','.join([chrom,gene,chain,str(pos_bin)])
    infor[chrom][pos_bin].append(key)
    trick.set3dict(gene_infor,chrom,key,'start',start)
    trick.set3dict(gene_infor,chrom,key,'end',end)
    trick.set3dict(gene_infor,chrom,key,'chain',chain)


def enh(pick_gene, espan):
    out_line = []
    if pick_gene:
        sort_dit =  sorted(list(pick_gene.items()), lambda x, y: cmp(x[1], y[1]), reverse= False )
        distance =  sort_dit[0][1]
        for each in sort_dit:
            distance =  each[1]
            if distance < espan :
                out_line.append(each[0])
                out_line.append(str(distance))
                return out_line
    return ['no match','']


def promoter_peak(pstart, pend, gstart, key):
    bed_set = set(range(start,end))
    gene_set = set(range(gstart - pspan, gstart + pspan))
    lst = bed_set.intersection(gene_set)
    if len(lst) >= 1 :
        elst = [str(i) for i in (len(lst), key)]
        sys.stderr.write('\t'.join(elst) + '\n')
        return True 
    else :
        return False 




for line in peak:
    if 'track' in line :
        continue
    line_arr = line.strip().split('\t')
    chrom = line_arr[0]
    if chrom == 'chr' :
        header = line_arr
        header.insert(3,'gene')
        header.insert(4,'mid2midLen')
        if not ignore_promoter :
            pass
#            header.append('PromoterPeak')
        print('\t'.join(header))
    if chrom not in chroms:
        continue
    start = int(line_arr[1])
    end = int(line_arr[2])
    out_line = []
    last = line_arr[0:3]
    pick_gene = {}
    pos_bin_start = get_bin(start,span)
    pos_bin_end = get_bin(end,span)
    lst = list(range(pos_bin_start,pos_bin_end))
    lst.append(pos_bin_end)
    for i in range(1,3):
        lst.append(pos_bin_start - i)
        lst.append(pos_bin_end + i)
    lst = set(lst)
    #if '\t'.join(line_arr[0:3]) == 'chr11\t54735607\t54737821':
    peak_index = False
    if 1 :
        for i in lst :
            if i in infor[chrom]:
                keys = infor[chrom][i]
                for key in keys:
                    chain = gene_infor[chrom][key]['chain']
                    gene_start = gene_infor[chrom][key]['start']
                    gene_end = gene_infor[chrom][key]['end']
                    nearest_infor = {}
                    if chain == '+':
                        abs_span = abs((start + end)/2 - gene_start)
                        if ignore_promoter :
                            if not promoter_peak(start, end, gene_start, key):
                                trick.set1dict(pick_gene,key,abs_span)
                            else :
                                peak_index = True 
                                break 
                        else :
                            trick.set1dict(pick_gene,key,abs_span)
                    elif chain == '-' :
                        abs_span = abs((start + end)/2 - gene_end)
                        if ignore_promoter :
                            if not promoter_peak(start, end, gene_end, key):
                                trick.set1dict(pick_gene,key,abs_span)
                            else :
                                peak_index = True 
                                break 
                        else :
                            trick.set1dict(pick_gene,key,abs_span)
                    if debug:
                        print(keys)
                        print(start,end,gene_start,gene_end)
                        print(pick_gene)
            if peak_index :
                break 
    else :
        continue
    if peak_index:
        sys.stderr.write(line)
        continue
    #pick_gene = set(pick_gene)
    out_line.extend(last)
    out_line.extend(enh(pick_gene,span))
    out_line.extend(line_arr[3:])    
    out.write('\t'.join(out_line) + '\n')
        
