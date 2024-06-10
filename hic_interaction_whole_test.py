#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
import math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from ningchao.nSys import trick,fix,dataframe,status,line_num
from ningchao.nSys import parse as Parse
from ningchao.nBio import chromosome, order, rheMac34
from tempfile import NamedTemporaryFile
from statsmodels.sandbox.stats.multicomp import multipletests

example = '''one normal matrix get 2M interaction results rotio'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( '-peirod', choices=['E50','E80','E90','E120','0M','4M','45Y','20Y'], help = 'peirods', nargs = '*' )
parser.add_argument( '-res', nargs = '?', type = int, help = 'span for the', default = 2000000 )

if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()


def loadDivisionRemote( arrs, span = 2000000 ):
    step, dit = span / 40000, {}
    for peirod in arrs:
        matrixs = arrs[peirod]['chrom']
        wd = arrs[peirod]['work_dir']
        for chrom in matrixs:
            signal_file = matrixs[chrom]
            sys.stderr.write( signal_file + '\n' )
            array = np.array(pd.read_csv( signal_file, header = None, sep = '\t'))
            shape = array.shape
            distal, local = 0, 0
            for i, line in enumerate( array ):
                if i + step >= shape[0]:
                    continue
                local_line_sum = np.take( line, np.arange(i, i + step ) ).sum()
                distal_line_sum = line.sum() - local_line_sum
                distal += float( distal_line_sum )
                local += float( local_line_sum )
            trick.dinit( dit, peirod, chrom, 'distal_sum', distal )
            trick.dinit( dit, peirod, chrom, 'local_sum', local )
    return dit

def plot( infor ):
    df, chroms = {}, []
    for peirod in infor:
        for chrom in infor[peirod]:
            distal_sum = infor[peirod][chrom]['distal_sum']
            local_sum = infor[peirod][chrom]['local_sum']
            # peirod, chrom, distal_sum, local_sum, distal_sum/local_sum 
            trick.dinit( df, peirod, chrom, [] )
            chroms.append( chrom )
            df[peirod][chrom].append( distal_sum/local_sum )
    df = pd.DataFrame( df )
    print(df)

def main():
    peirods = args.peirod
    arrs = rheMac34.hic( peirods ).dir2dit( )
    infor = loadDivisionRemote( arrs, args.res )
    plot( infor )
if __name__ == '__main__':
    main()
















