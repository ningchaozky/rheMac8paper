#!/usr/bin/env python
import os
import sys
import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')
import pandas as pd
from pandas import read_csv
import itertools
import matplotlib.pyplot as plt
import ningchao.nSys.fix as fixKit
from scipy.stats import wilcoxon, levene, ttest_ind
from collections import defaultdict
import seaborn as sns;sns.set(color_codes=True)
import ningchao.nSys.trick as trKit

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='print useage for script', formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('matrix_labal_pair', nargs = '+', help ='matrix and label pair' )
parser.add_argument('-y', nargs='?', type = float, help ='ymax for the boxplot', default = 25)

if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()

lst = iter(args.matrix_labal_pair )
dfs = []
for i in lst:
    data = read_csv( i, sep = '\t', header = 0, index_col = 0)
    data['label'] = next( lst  )
    dfs.append ( data.melt( id_vars=['label'] ) )
df = pd.concat( dfs, ignore_index = True )
#ax=sns.boxplot(x='variable',y='value', hue = 'label', data = df, showfliers = False )
plt.ylim((0, args.y))
df['value'] = df['value'].apply( lambda x: x if x < 80 else 80 )
df = df[df['value'] > 2]
ax = sns.violinplot( x='variable', y='value', hue = 'label', data = df)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
ax.set_ylim([0, 150]) 
plt.savefig(('prefix' + '.pdf'), format='pdf',  bbox_inches='tight', pad_inches=+1 )

