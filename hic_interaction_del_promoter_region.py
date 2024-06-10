#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,system
from tempfile import NamedTemporaryFile
example = ''' gene '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'gene', nargs='?', help ='gene file from bw_matrix_insert_gene.py. chr,start,end,gene,[values...]' )
parser.add_argument( '-p', nargs='?', help ='promoter file which was used to exclude the enhancer peaks', default = '/home/ningch/data/genome/rheMac8/exon/gene.sort.tss1K')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    gene = args.gene
    promoter_file = args.p
    return gene,promoter_file

def intersect( gene, promoter_file ):
    genefh, infor = NamedTemporaryFile(prefix='gene_intersect', suffix='.bed', dir='/tmp', delete = False ), {}
    for line in open(gene):
        if 'start' in line:
            continue
        genefh.write(line)
    genefh.close()
    cmd = 'bedtools intersect -a %s -b %s -wo' % (promoter_file, genefh.name)
    sys.stderr.write(cmd + '\n')
    stdout,stderr,returncode = system.run('bedtools intersect -a %s -b %s -wo' % (promoter_file, genefh.name))
    for line in stdout:
        key = '\t'.join(line.strip().split('\t')[6:9])
        trick.dit(infor).set(key, 0)
        infor[key] += 1
    return infor


if __name__ == '__main__':
    gene,promoter_file = parse( args )
    del_peaks = intersect( gene, promoter_file )
    gfh = open( gene )   
    for line in gfh:
        line_arr = line.strip().split('\t')
        key = '\t'.join(line_arr[0:3])
        if key in del_peaks:
            continue
        print(line, end=' ')






