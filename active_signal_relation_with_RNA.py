#!/usr/bin/env python3
import os
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
plt.tight_layout()
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from collections import defaultdict
#plt.savefig( figname, dpi=250, transparent=True, facecolor=fig.get_facecolor(), edgecolor='none')
#plt.tick_params( axis = 'both', left = True, labelleft = True, which = 'both', bottom = True, top = False, labelbottom = True, direction = 'in' )
#sns.despine( ax = ax[2]  )#remove top and right axis
#subplot define, also polar plot
#with sns.axes_style("whitegrid", {'xtick.top': True, 'ytick.left': True, 'axes.grid': False, 'ytick.color': 'black', 'ytick.direction': 'in' }):
    #fig, ax = plt.subplots( 6, figsize=(10,40), subplot_kw=dict(projection=None))
	#fig,ax = plt.subplots( 3, figsize=(10,30), sharex=True, sharey=True, subplot_kw=dict(projection='polar')); ax[0],ax[1]
plt.style.use('ggplot')
plt.subplots_adjust(wspace=0, hspace=0)
import seaborn as sns;sns.set(color_codes=True)
sns.set_style("ticks")
#sns.barplot( palette="Set3" )
#fig = plt.figure( figsize=( 18, 24) )
import argparse
from ningchao.nSys import trick, system
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='?', help = 'active signal with period append')
parser.add_argument( 'rna', nargs='?', help = 'RNA data')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()





def parse_signal(tab):
    periods = system.dir.periods( prefrontal_raw = True )
    infor = defaultdict( lambda : defaultdict( list ) )
    with open( tab ) as f :
        header = next(f).rstrip().split('\t')
        header_match = [ system.dir.str_map_period(i) for i in header ]
        for line in f :
            line_arr = line.rstrip().split('\t')
            period = line_arr.pop(-1)
            symbol = line_arr[0].split('.')[0]
            for p,v in zip(header_match, line_arr ):
                if p == period :
                    infor[period]['exp'].append( float(v) )
    return infor

def rna_pick_no_exp( infor, rna):
    rna_exp = defaultdict( list )
    with open( rna ) as f :
        header = next(f).rstrip().split('\t')
        header_match = [ system.dir.str_map_period(i) for i in header ]
        for line in f :
            line_arr = line.rstrip().split('\t')
            symbol = line_arr[0].split('.')[0]
            rna_exp[ symbol ].append( line )
    for symbol in rna_exp :
        rna_exp[symbol] = trick.lst( rna_exp[symbol] ).pick( pt = 'sum' )
    for period in infor :
        for rna_symbol in rna_exp :
            if rna_symbol not in infor[period] :
                infor[ period ][ 'noexp' ].append( rna_exp[rna_symbol][header_match.index(period)])
    return infor
def main():
    infor = parse_signal(args.tab)
    infor = rna_pick_no_exp(infor, args.rna)
    lst = []
    for period in infor :
        exp = infor[period]['exp']
        noexp = infor[period]['noexp']
        for e in exp :
            lst.append([ period, 'active', np.log10(e+1) ])
        for noe in noexp :
            lst.append([period, 'noActive', np.log10(noe+1)])
    df = pd.DataFrame( lst, columns = ['period','typ','log10(tpm)'] )
    ax = sns.violinplot( x = 'period', y='log10(tpm)', hue = 'typ', data = df )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.ylim((0, 6))
    plt.savefig( 'active_signal_relation_with_RNA.pdf', dpi=250, transparent=True, edgecolor='none')

if __name__ == '__main__':
    main()



















