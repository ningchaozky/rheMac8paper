#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import pandas as pd
from ningchao.nSys import trick, system
example = '''ordered cortab and pval tab get the ran ranked'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'cor', nargs = '?', help ='sorted correrlation matrix' )
parser.add_argument( 'pval', nargs = '?', help ='pval order the same with cor' )
parser.add_argument( 'tpm', nargs = '?', help = 'expression vales include all the genes')
parser.add_argument( '-p', nargs = '?', help = 'prefix for the cluster', default = 'select' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    return args.pval, args.cor, args.tpm, args.p

def head( df ):
    fh, out = open( df ), {}
    header = fh.next().rstrip('\n').split('\t')
    return header, fh 

def cluster( prefix ):
    cwd, out = os.getcwd(), {}
    fls = system.dir( cwd ).fls(pattern = '%s.*txt' % prefix, fpattern = 'Cytoscape', level = 0 )
    for fl in fls:
        color = 'ME' + os.path.basename( fl ).split('.')[ -2 ]
        genes = [ i.strip().split('\t')[1] for i in open(fl).readlines() if 'modProbes' not in i ]
        trick.dinit( out, color, genes )
    return out

def cor_dit( cor ):
    header, fh = head( cor )
    out = {}
    for line in fh:
        line_arr = line.strip().split('\t')
        color = line_arr.pop(0)
        line_arr = [ float(i) for i in line_arr ]
        if trick.lst(line_arr).uniq() == [ 0 ]:
            continue
        out.update( {color: dict(list(zip( header, [ float(i) for i in line_arr ])))} )
    return out

def meaning( cor, pval ):
    corDit = cor_dit( cor )
    header, fh = head( pval )
    meaning = { }
    for line in fh:
        line_arr = line.strip().split('\t')
        color = line_arr.pop(0)
        row_dit = dict( list(zip( header, [ float(i) for i in line_arr ] )) )
        row_dit_sort = sorted( list(row_dit.items()), lambda x, y: cmp(x[1], y[1]), reverse=False)
        for i,j in row_dit_sort:
            if j < 0.05 and corDit[color][i] > 0 :
                trick.dinit( meaning, color, [])
                meaning[color].append( i )
    return meaning

def rna_meaning( tpm ):
    dit = cor_dit( tpm )
    for gene in dit:
        print(gene)


def main( ):
    pval, cor, tpm, prefix = parse( args )
    cluster_with_gene = cluster( prefix )
    sense = meaning( cor, pval )
    rna = cor_dit( tpm )
    for color in sense:
        genes = cluster_with_gene[ color ]
        select_peirod = sense[ color ]
        for gene in genes:
            if gene not in rna:
                continue
            arr = sorted(list(rna[ gene ].items()), lambda x, y: cmp(x[1], y[1]), reverse = True )
            cut_arr = arr[ 0 : len(select_peirod) ]
            #print gene, select_peirod, color, cut_arr
            RNA_peirod = [ i.split('.')[0] for i,j in cut_arr ] 
            select_peirod = [i.split('.')[0] for i in select_peirod ]
            if trick.lst(RNA_peirod).sublist(select_peirod):
                print(color + '\t' + gene)
if __name__ == '__main__':
    main()
































