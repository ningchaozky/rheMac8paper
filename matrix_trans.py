#!/usr/bin/env python3
import os
import re
import csv
import ast
import sys
import math
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from numpy import mean,std
from scipy.stats import variation 
from scipy.spatial import distance
from ningchao.nSys import env,trick,lineKit,fileKit,fix,quantile_norm,excel,xls
#fig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#subplot define, also polar plot
fig,ax = plt.subplots( 1, figsize=(10,30), sharex=True, sharey=True, subplot_kw=dict(projection=None))
plt.style.use('ggplot')
plt.subplots_adjust(wspace=0, hspace=0)
import seaborn as sns;sns.set(color_codes=True)
sns.barplot( palette="Set3" )


parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='print useage for script')
parser.add_argument( 'matrix', nargs='*', help = 'matrix|matrixs you want to calculate the value. K27.mini,index_col,1,header,1 ' )
parser.add_argument( '-method','-m', choices = [ 'pow2', 'exp', 'entropy', '+1log2', '+0log10', 'zcol', 'zrow', 'n1','normalCol','quantile_normalization','+','-','*','/','cv','js' , 'ts','sum1','rev', 'assign', 'qcut','binary','discrete'], help ='js: jensenshannon, ts: transposion, sum1: sum all is 1, rev mean reverse all the value', required = True )
parser.add_argument( '-rev', action = 'store_true', help = 'revese the normal. Like K27')
parser.add_argument( '-dtype','-d', choices = ['discrete','continuous'], help = 'revese the normal. Like K27', default = 'continuous')
parser.add_argument( '-cut', nargs='?', help = 'cutoff for express value. if there is no more than this value all will be set to zero Like K27.max value file' )
parser.add_argument( '-cutoff', nargs='?', help = 'cutoff for tab max value', default = -1000000, type = float )
parser.add_argument( '-value','-v', nargs='+', help = 'value can use any want to get. eg: delete more than 1 entropy')
parser.add_argument( '-fcol','-f', nargs='+', help = 'fcol not deal', type = int)
parser.add_argument( '-weight','-w', nargs='?', help = 'weight for the matrix', default = 1.0, type = float )
parser.add_argument( '-output','-o', nargs='?', help ='output file for the entropy value' )
#parser.add_argument( '-cat', nargs = '*', help ='val_fiter(entropy fiter),single_fiter(fignal val samll than this will ignore)', default = [])
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def parse ( **kwargs):
    #matrix, method, rev, cut, weight, cutoff
    ofh = kwargs.get('output') and open( kwargs.get('output'),'w') or sys.stdout
    dfs = xls.excel( kwargs.get('matrix')).dfs( row_cutoff = kwargs.get('cutoff'))
    if kwargs.get('cut') :
        max_df = xls.excel( kwargs.get('cut') ).dfs()
        if len(max_df) == 1 : max_df = max_df[0]
    else :
        max_df = pd.DataFrame()
    kwargs.update({'dfs': dfs, 'method': args.method, 'ofh': ofh})
    return kwargs
def correct( express ): # matrix for the sample absolute express
    dfs = xls.excel ( express )
    print(dfs)



if 0 :
    if method == 'normalCol':
        for j,each in enumerate(line_arr):
            line_arr[j] = float(line_arr[j]) * (max_sum_value/ float(zinfor[str(j)]['sum']))


def deal( df, **kwargs ):
    method, ofh = [ kwargs.get(i) for i in [ 'method', 'ofh'] ]
    if method == 'quantile_normalization':
        qdf = quantile_norm.column(df).round(4)
        qdf.to_csv( ofh, sep = '\t', header = True, index = True, index_label = qdf.index.names )
    elif method == 'zrow':
        zdf = pd.DataFrame( stats.zscore(df, axis=1, ddof=1).round(4), index = df.index.values)
        zdf.columns = df.columns
        if args.rev :
            zdf = (0 - zdf).round(4)
        zdf.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif method == 'zcol':
        vmax = float( args.value[0] )
        for peirod in df.columns:
            df[peirod] = np.where( df[peirod] > vmax, vmax, df[peirod] )
        index_values = list(df.index.values)
        if isinstance( index_values[0], tuple ):
            index = pd.MultiIndex.from_tuples( df.index.values )
        elif isinstance( index_values[0], str):
            index = index_values
        index_label = df.index.names
        df = pd.DataFrame( stats.zscore(df, axis = 0, ddof=1).round(4), index = index, columns = df.columns )
        #zdf.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names,  escapechar=" ", quoting = csv.QUOTE_NONE)
        if len( args.value ) >= 2:
            normal_val = float( args.value[1] )
            df = df + normal_val
        df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = index_label )
    elif method == 'entropy':
        #print stats.entropy(df)
        print ('#Delete all 0 line')
        df = df[~df.eq( 0 ).all(1)]
        df = df.T
        pA = df / df.sum( axis = 0 )
        print ('#Now is entropy method, row cutoff is {}, use flot method'.format(kwargs.get('cutoff')))
        if kwargs.get('dtype') == 'continuous':
            Shannon2 = -np.sum( pA*np.log2(1/df), axis = 0 )
        else :
            Shannon2 = -np.sum( pA*np.log2(1/pA), axis = 0 )
        Shannon2 = Shannon2.sort_values( 0, ascending = False )
        if args.value :
            entropy_cut = float( kwargs['value'][0] )
            Shannon2 = Shannon2[Shannon2 > entropy_cut]
            print ( '#{}'.format(Shannon2.shape), file = sys.stderr )
        Shannon2.plot.line()
        plt.savefig( 'testfafaf.pdf', dpi=250, transparent=True, facecolor=fig.get_facecolor(), edgecolor='none')
        Shannon2.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'log' in method:
        p = re.compile(r'\+|log')
        av,lv = [ int(i) for i in p.split(method) if i ]
        ldf = df + av
        ldf = eval( 'np.log%d(ldf)' % lv )
        ldf.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'cv' in method:
        var = variation(df, axis=0)
        cv =  pd.DataFrame(np.std(df, axis = 1) / np.mean(df, axis = 1))
        cv.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'normalTo1' in method or 'n1' in method:
        adf =  np.array( df )
        adfmin = adf.min( axis = 1 )
        scale = adf.max(axis = 1) - adf.min(axis = 1)
        scale[ scale == 0 ] = 1
        relative = adf - np.outer(adf.min( axis = 1 ), np.ones(adf.shape[1]))
        val = pd.DataFrame( relative / np.outer(scale, np.ones(adf.shape[1])), index = df.index.values).round(4)
        val.columns = df.columns
        #cutoff for selfvalue c1|c5
        if not args.cut:
            max_df = df
        cut = max_df.max( axis = 1 ) >= 0
        #cut = adf.max( axis = 1 ) > cutoff
        val['cut'] = cut
        val.loc[ ( val.cut == False ), val.columns ] = 0.0
        #val.loc[ args.rev & ( val.cut == True ), val.columns ] = (1 - val ) .round(4)
        if args.rev: val.loc[ ( val.cut == True ), val.columns ] = (1 - val ) .round(4)
        val.pop('cut')
        val = val * args.weight
        val.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names)
    elif 'js' in method:
        pdit = peirod_matrix( df.columns )
        gene = df.ix[ 1573 ]
        print(gene)
        for peirod in pdit:
            print(peirod, distance.jensenshannon( gene, pdit[peirod] ))
    elif 'ts' in method:
        df.T.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'sum1' in method:
        df = df.T
        normal = df.div( df.sum( axis = 1), axis = 0 )
        normal = pd.DataFrame( normal, index = df.index.values)
        normal.T.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'exp' in method :
        np.exp(df).to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif 'pow2' in method :
        df.pow(2).to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif method in ['+','-','*','/']:
        if not kwargs.get('value'):
            exit ('#Must set value for {}'.format(method))
        else :
            df = eval( 'df{}{}'.format( method, kwargs.get('value')[0]))
            df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif method in ['assign']:
        vals = iter( args.value  )
        for col in vals:
            df[col] = float( next( vals ) )
        df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names  )
    elif method in ['qcut']:
        times = int(args.value[0])
        for peirod in df.columns:
            qcut_colmun = '{}.qcut'.format(peirod)
            df[ qcut_colmun ] = df[peirod].rank(method='first')
            df[qcut_colmun] = pd.qcut(df[qcut_colmun].values, times).codes
            df[qcut_colmun] = np.where( df[peirod] == 0, 1, df[qcut_colmun] )
        df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif method in ['binary']:
        if not args.value:
            print ('##Must set args.value for binary cutoff', file = sys.stderr)
            exit()
        condition_list = [ df >= float( args.value[0]  ), df < float( args.value[0]  )]
        choices_list = [ 1, 0 ]
        df = pd.DataFrame( np.select(condition_list, choices_list), index = df.index, columns = df.columns)
        df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names )
    elif method in ['discrete']:
        df = df.astype('float64')
        labels  = [ 'low','median','hight']
        labels  = [ 0, 1, 2]
        df = df.apply( pd.cut, axis = 1, bins = [ 0, 0.2, 0.8, 1.5 ], right = False, labels = labels)
        df.to_csv( ofh, sep = '\t', header = True, index = True, index_label = df.index.names  )
    else :
        print ( 'Wrong input method {}'.format(method) )



def jitter(a_series, noise_reduction=1000000):
    return (np.random.random(len(a_series))*a_series.std()/noise_reduction)-(a_series.std()/(2*noise_reduction))



def peirod_matrix( columns ):
    pdit = {}
    for i,peirod in enumerate( columns ):
        arr = [0] * len( columns )
        arr[i] = 1
        trick.dinit(pdit, peirod, arr)
    return pdit

def deals( **kwargs ):
    dfs, method = kwargs.get('dfs'), kwargs.get('method')
    if method in ['+','-','*','/']:
        xls.df_operate( dfs ).opt( kwargs.get('ofh'), typ = method )



if __name__ == '__main__' :
    conf = parse ( **vars( args ) )
    if len(conf['dfs']) == 1 :
        df = conf['dfs'][0]
        #deal( df, method, ofh, rev, max_df = max_df, cutoff = cut )
        if args.fcol:
            df.drop(df.columns[args.fcol],axis=1,inplace=True)
        deal( df, **conf )
    else :
        deals( **conf )





