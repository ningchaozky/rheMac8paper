#!/usr/bin/env python3
import os
import sys
import argparse
import random
import math
import pickle
import numpy as np
import pandas as pd
from scipy.stats import ranksums, ttest_rel, ttest_ind, wilcoxon, mannwhitneyu
from ningchao.nSys import trick, system, fix
from collections import defaultdict
desc = ''' merge output of hic_interaction_si_plot_v2.py results'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('directory', nargs='?', help = 'output directory of hic_interaction_si_plot_v2.py')
parser.add_argument('-dp', nargs='?', help = 'directory name for find uniq gene name, 10K|40K', default = '10K')
parser.add_argument('-c', nargs='?', help = 'cutoff for signal', default = 10, type = float)
parser.add_argument('-f', choices=['pickle','tab'], help = 'force generate results')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def create_int_defaultdict():
    return defaultdict( int )
def create_list_defaultdict():
    return defaultdict( list )


def get_infor( **kwargs ) :
    if os.path.exists( pick_file ) and os.path.getsize( pick_file ) > 1000 and 'pickle' not in args.f:
        print ('Already own {}... loading it...'.format( pick_file ), file = sys.stderr)
        return pickle.load( open(pick_file, 'rb'))
    else :
        signal_infor, num_infor = defaultdict( create_list_defaultdict ), defaultdict( create_int_defaultdict )
        fls = system.dir( args.directory ).fls( r'{}.*tab$'.format(args.dp), depth = [2,8] )
        i = 0
        for fl in fls:
            i += 1
            if not i % 1000:
                print ( '!!!gene: ', i, fl, file = sys.stderr)
            fl_arr = fl.split('/')
            peirod, bfn = fl_arr[1], fl_arr[-1]
            bfn_arr = bfn.split('.')
            key = '.'.join( bfn_arr[3:6] )
            with open( fl ) as f :
                next(f)
                for line in f:
                    line_arr = line.strip().split('\t')
                    signal = float( line_arr[-1] )
                    signal_infor[peirod][ key ].append( signal )
                    if signal > args.c:
                    #if 1 :
                        num_infor[peirod][ key ] += 1
        data = { 'num': num_infor, 'signal': signal_infor }
        with open( pick_file, 'wb'  ) as f :
            pickle.dump( data, f, protocol = pickle.HIGHEST_PROTOCOL )
        return data
if __name__ == '__main__':
    pick_file = fix.fix( args.dp ).append( 'c{}.pickle'.format(args.c) )
    infor = get_infor( **vars( args ) )
    signal_infor, num_infor = infor.get('signal'), infor.get('num')
    order = system.dir.sort( list(signal_infor.keys()), up = '', down = '')
    combs = trick.lst.combs( order, k = 2, sort = 'peirod', up = '', down = '')
    if 'tab' in args.f:
        for peirod in signal_infor :
            for key in signal_infor[peirod] :
                val = [ i for i in signal_infor[peirod][key] if i > args.c  ]
                #signal_infor[peirod][key] = np.sum( [ i for i in signal_infor[peirod][key] if i > args.c ])
                if val :
                    signal_infor[peirod][key] = max( val )
                else :
                    signal_infor[peirod][key] = 0
    ntab, stab = '{}.c{}.num.tab'.format(args.dp, args.c), '{}.c{}.signal.tab'.format(args.dp, args.c)
    if not os.path.exists( stab ) or 'tab' in args.f:
        print ('generate tab...')
        df_signal, df_num = pd.DataFrame( signal_infor ).fillna(0), pd.DataFrame( num_infor ).fillna(0)
        df_signal, df_num = df_signal.reindex( order, axis=1), df_num.reindex( order, axis=1 )
        random = ''.join( random.sample( [ chr(i) for i in range(65,91)  ], k = 4  ) )
        df_signal.to_csv( stab, sep = '\t', index_label = 'symbol')
        df_num.to_csv( ntab, sep = '\t', index_label = 'symbol')
        exit('Signal merge one, do not continue...')
    #symbols, kwargs = [], vars( args )
    if 0 :
        for key in symbols:
            if 'EOMES' in key :
                #for p1,p2 :
                p1, p2 = signal_infor['E50'][key], signal_infor['E120'][key]
                print (key, ranksums( p1, p2 ))
                #print (key, ttest_rel ( p1, p2 ) )
                print (key, ttest_ind( p1, p2 ) )
                #signal_infor[ peirod  ][ key ] = np.median( signal_infor[ peirod ][key] )
    signal_infor = infor.get('signal')
    infor = defaultdict( lambda : defaultdict( float ) )
    for peirod1, peirod2 in combs:
        for key in signal_infor[peirod1]:
            if 'GLI3' in key :
                p1, p2 = signal_infor[ peirod1 ][ key ], signal_infor[ peirod2 ][ key ]
                #p1, p2 = [ i for i in p1 if i > args.c ], [ i for i in p2 if i > args.c ]
                #pvalues = [ ranksums( p1, p2  ).pvalue, ranksums( p1, p2  ).pvalue, ttest_rel ( p1, p2  ).pvalue, ttest_ind( p1, p2).pvalue #]
                pvalues = [ ranksums( p1, p2  ).pvalue, ttest_ind( p1, p2).pvalue, mannwhitneyu( p1, p2   ).pvalue ]
                pvalues = [ ranksums( p1, p2  ).pvalue, ttest_ind( p1, p2).pvalue, mannwhitneyu( p1, p2   ).pvalue ]
                infor[peirod1][peirod2] = 1 - min( pvalues )
                infor[peirod2][peirod1] = infor[peirod1][peirod2]
        infor[peirod2][peirod2] = 0
        infor[peirod1][peirod1] = 0
    df = pd.DataFrame(infor)
    df = df.reindex( order, axis=0  )
    df.to_csv( 'eomes.tab', sep = '\t', index_label = 'symbol' )
    os.system('heatmap_seaborn.py eomes.tab -cmap Greys')
