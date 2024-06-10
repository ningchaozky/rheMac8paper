#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('prefix', nargs='?', help ='prefix')
parser.add_argument('-g', nargs='*', help ='genes', required = True )
parser.add_argument('-stw', action='store_false', help ='startswith' )
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()
files = os.listdir('.')
cwd = os.getcwd()
prefix = args.prefix
files = [ i for i in files if i.startswith(prefix) and i.endswith('.txt') and 'Cytoscape' not in i ]
infile = cwd+'/' + args.g[0]

def check( genes ):
    for gene in args.g:
        fls = []
        for fl in files:
            fh = open(fl)
            for line in fh :
                if 'modProbes' in line:
                    continue
                lgene = line.strip().upper().split('\t')[1]
                gene = gene.upper()
                index = args.stw and gene in lgene or lgene.startswith(gene)
                #print index,  gene in lgene, lgene.startswith(gene)
                if index :
                    if fl not in fls:
                        fls.append('\t'.join([gene, fl, line]))
        print('\n'.join( fls ))

if os.path.exists(cwd+'/' + args.g[0]) :
    gfh = open(infile)
    genes = []
    for each in gfh:
        genes.append(each.strip().split('\t')[0])
    gfh.close()
    check(genes)

else :
    check(args.g)
























