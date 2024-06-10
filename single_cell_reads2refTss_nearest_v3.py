#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import pysam
import queue
from ningchao.nSys import trick,fix
from multiprocessing.dummy import Pool as ThreadPool

example = '''bam refbed'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('bam', nargs='?', help ='mapping bam file')
parser.add_argument('bed', nargs='?', help ='bed file')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

bed = args.bed
bam = args.bam
span = 2000
refDit = {}
cfix = fix.fix(os.path.basename(bam))
line_num = 0
def ref2dit(bed):
    dit = dict()
    bfh = open(bed)
    next(bfh)
    for line in bfh:
        line_arr = line.strip().split('\t')
        chrom = line_arr[2]
        chain = line_arr[-4]
        transcript_tss_bin = str(int(line_arr[-3])/span)
        key = '.'.join([chrom, chain, transcript_tss_bin])
        dit[key] = line.strip()
    return dit

def check(chrom, chain, start, end, dit):
    sbin = int(start)/span
    tmp_trans = {}
    if chain == '+':
        sbin = int(start)/span
    else :
        sbin = int(end)/span
    for each in range(sbin-1,sbin+2):
        key = '.'.join([chrom, chain, str(each)])
        if key in dit:
            line_arr = dit[key].split('\t')
            transId = line_arr[1]
            trans_tss = int(line_arr[-3])
            if chain == '+':
                nspan = abs( int(start) - trans_tss )
            else :
                nspan = abs( int(end) - trans_tss )
            if nspan < span :
                tmp_trans.update({transId:nspan})
    if tmp_trans:
        return sorted(list(tmp_trans.items()), lambda x, y: cmp(x[1], y[1]), reverse= False )[0]
    else :
        return False

tss = ref2dit(bed)

samfile = pysam.AlignmentFile(bam, "rb")
outPutNotMatch = pysam.AlignmentFile( cfix.append('noMatch.bam'), "wb", template=samfile )
transDit = {}
def line( reads):
    dit = transDit
    for read in reads:
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_start
        reverse = read.flag & 16
        if reverse :
            chain = '-1'
        else :
            chain = '1'
        transcript = check(chrom, chain, start, end, tss)
        if transcript:
            tid,tspan = transcript
            trick.set1dict( dit, tid, 0)
            dit[tid] += 1
        else :
            outPutNotMatch.write( read )
    return dit
def log_result(result):
    results.append(result)

reads = []
transDit = {}
q = queue.Queue()
results = []
pool = ThreadPool(5)
for read in samfile:
    line_num += 1
    if line_num % 10000 == 0:
        #record
        sys.stderr.write('line num is: %d\n' % line_num)
        
        result = pool.apply_async(line, args = (reads, ))
        results.append(result)
        reads = []
        reads.append(read)
    else :
        reads.append(read)

    if line_num % 50000 == 0 :
        pool.close()
        pool.join()
        
        print([ i.get() for i in results ])
        pool = ThreadPool(5)
        results = []
        result = pool.apply_async(line, args = (reads, ))
        results.append(result)


outPutNotMatch.close()
ofh = open(cfix.append('infor'),'w')
for key in transDit:
    ofh.write('\t'.join([key, str(transDit[key])+'\n']))
ofh.close()






