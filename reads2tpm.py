#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
from ningchao.nSys import trick, system
from bioinfokit.analys import norm, get_data
example = '''this need gene_length file, which can get from gtftools, excel_concat.py and excel_mini.py pipline\n'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'reads', nargs = '?', help = 'reads counts matrix with mean gene_length')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



df = pd.read_csv( args.reads, sep = '\t', header = 0, index_col = 0)
nm = norm()
nm.tpm(df=df, gl = 'mean')
tpm = nm.tpm_norm.round(4)
tpm.to_csv( '{}.tpm'.format( args.reads), sep = '\t' )
























