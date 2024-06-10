#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import re
import sys
import argparse
import pickle
import random
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.tight_layout()
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
from collections import defaultdict
from ningchao.nSys import trick,system
from ningchao.nBio import entrez,geneKit
from ningchao.nBio import rna,geneKit,neuron,entrez
import numpy as np
import math
from collections import defaultdict
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import seaborn as sns;sns.set(color_codes=True)
sns.set_style("ticks")
sns.barplot( palette="Set3" )
import matplotlib.pyplot as plt
from matplotlib import dates
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from ningchao.nSys import trick, fix, parse, excel
from ningchao.nBio import neuron,rheMac,geneKit
desc='''!!!you can edit the label_deal fun to show the elegant label\n{} include the alias information'''.format(os.path.basename(sys.argv[0]))
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep, description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='*', help = 'tab and col for the gene. if len big than 3 will use like gene.ESNG0000')
parser.add_argument( '-s', nargs='?', help = 'pattern for select', default = '.*')
parser.add_argument( '-f', nargs='?', help = 'pattern for select', default = 'fdafasdfasdfdasf122212' )
parser.add_argument( '-n', action = 'store_true', help = 'only the name ' )
parser.add_argument( '-ih', action = 'store_false', help = 'ignore header ' )
parser.add_argument( '-order', nargs = '?', help = 'order file' )
parser.add_argument( '-not_include_alias','-nia', action = 'store_false', help = 'not parse the input gene in alias' )
parser.add_argument( '-cut', nargs = '?', help = 'pvalue cut off', default = 1.3, type = float )
parser.add_argument( '-diff_num', nargs = '?', help = 'diff gene num you set', type = int )
parser.add_argument( '-pt', choices = ['h','v','heatmap'], help = 'h|v', required = True)
parser.add_argument( '-cc', action = 'store_true', help = 'h|v')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()


def label_deal( dit ):
    out = defaultdict( list )
    for k,v in dit.items() :
        karr = iter(k.split('fdafasdfasdfafaf'))
        if 1 :
            if 1 :
                key = next(karr)
            if 0 :
                for each in karr :
                    if re.search('C\d+', each ) :
                        key = '|'.join([each, next(karr)])
        else :
            key = k
        out[key] = dit[k]
    return out

def sort_plot( label, values):
    lst = zip( label, values )
    #lst_sorted = sorted( lst, key = lambda x: int( x[0].split('|')[0].replace('C','')), reverse = False)
    lst_sorted = sorted( lst, key = lambda x: x[1], reverse = False)
    return [ i[0] for i in lst_sorted ], [ i[1] for i in lst_sorted ]



def symbol_deal( symbol, **kwargs):
    p2 = re.compile(r'\.|_|\|')
    p,fp = kwargs.get('replace_pattern'), kwargs.get('fp')
    #if not kwargs.get('alias').get( p.sub( '', symbol) )
    alias = kwargs.get('alias').get( symbol )
    if alias :
        symbol = p.sub( '', kwargs.get('alias').get( symbol ))
        symbol = p2.split(symbol)[0]
    symbol = symbol.strip().upper()
    return symbol
def parse( args ):
    kwargs = vars( args )
    kwargs.update({'replace_pattern': re.compile(r'\.\d+$|\-AS\d+$|\_\d+$|\-\d+$', re.I), 'select_pattern': re.compile(r'{}'.format(args.s), re.IGNORECASE)})
    kwargs.update({'alias': entrez.entrez(), 'fp': re.compile(r'{}'.format(args.f, re.I)) })
    infor, tmpinfor = defaultdict( list ), defaultdict( set )
    fls = system.dir('/home/soft/data/genome/rheMac8/neuron/').fls('init.pickle', depth = 10)
    for fl in fls :
        print ('Load: {}'.format( fl )  )
        with open( fl, 'rb') as handle:
            b = pickle.load(handle)
            b = { k:v for k,v in b.items() if 'GO:' not in k }
            #if 'NenadSestan_SpatioTemporal' in fl:
            #    for each in b.keys() :
            #        print ( fl, each )
            #        print ( b[each] )
            tmpinfor.update(b)
    if args.s :
        for sp in args.s.split('|'):
            ps, fp = re.compile(sp, re.I), kwargs.get('fp')
            for k,v in tmpinfor.items():
                if ps.search( k ) and not fp.search( k ):
                    print ( 'select', k, len(v), file = sys.stderr )
                    infor[ k ] = [ symbol_deal(i, **kwargs) for i in tmpinfor[k] ]
    else :
        for k,v in tmpinfor.items():
            print ( k, len(v) )
    return infor, kwargs, kwargs.get('alias')

def input_symbols( alias ):
    lst = []
    symbols = defaultdict( list )
    tabs = excel.xls( args.tab ).parse( simple = False)
    for fl, flkwargs in tabs:
        print ( flkwargs )
        flkwargs.pop('header')
        flkwargs.pop('usecols')
        if args.ih : flkwargs.update({'header': None })
        df = pd.read_csv(fl, **flkwargs)
        fl_symbols_deal = [ i.split('.')[0] for i in list( df.index ) ]
        for symbol in fl_symbols_deal:
            symbol = symbol.upper()
            fl_symbol_deal = alias.get( symbol )
            if fl_symbol_deal :
                symbols[fl].append( fl_symbol_deal )
            else :
                symbols[fl].append( symbol )
    return symbols


def soure_only_num( source_cellType_dit ):
    source_cellType_dit_only_num = {}
    for key in list(source_cellType_dit.keys()):
        gene_num = len(source_cellType_dit[key])
        #key = key.replace('xcx.Given.Age.','')
        trick.dinit( source_cellType_dit_only_num, key, gene_num)
    return source_cellType_dit_only_num

def tab_cellType( tab ):
    genes = []
    col = len(tab) > 1 and int(tab[1]) - 1 or 0
    for line in open(tab[0]):
        if 'modProbes' in line:
            continue
        line_arr = line.strip().split('\t')
        name_sps = line_arr[col].split('.')
        if 0 :
            rmIndex = trick.lst( name_sps ).index( ['copy'], regular = True )
            gene = '.'.join([ v for i,v in enumerate(name_sps) if i not in rmIndex ])
        else :
            gene = name_sps[0]
        genes.append( gene )
    return genes


def prepare( genes, soure_dit, **kwargs):
    soure_only_num_dit, dit = soure_only_num( soure_dit ), defaultdict( float )
    for typ in soure_dit:
        whole_num = len( soure_dit[ typ ] )
        diff_num = args.diff_num if args.diff_num else len( genes )
        intersect_genes = list( set([ symbol_deal(i,**kwargs) for i in soure_dit[ typ ]]) & set([ symbol_deal(i,**kwargs) for i in genes]) )
        gene_num = len( intersect_genes )
        if not gene_num :
            continue
        #pvalue = 1 - stats.hypergeom( 25000,  whole_num, diff_num ).cdf( gene_num )
        pvalue = stats.hypergeom( M = 25000, n = whole_num, N = diff_num ).sf( gene_num - 1)
        print ( typ, pvalue, 'whole_num:', whole_num, 'diff_num:', diff_num, 'intersect_genes_num: ', intersect_genes)
        #default match
        pvalue = -math.log10( pvalue ) if pvalue else 10
        #if typ == 'PrimateSpcificGenes':
        #    pvalue = 2
        if intersect_genes and pvalue > args.cut :
            print('postive:', typ, 'length for diff:', diff_num, pvalue, intersect_genes )
            dit[typ] = pvalue
    return dit

def goid2desc( lst ):
    fh, dit, out = open('/home/ningch/data/GO/htmls.gos.append'), {}, []
    for line in fh:
        line_arr = line.strip().split('\t')
        goid, desc = line_arr[0], line_arr[-1]
        trick.dinit( dit, goid, desc )
    for each in lst :
        if 'GO:' in each:
            if each in dit:
                out.append(dit[each])
                continue
        out.append( each )
    return out


def plot( plot_dit, prefix, bar = 'v'):
    '''h|v'''
    #tab = args.tab
    x, y= [], []
    for tab in plot_dit:
        for typ, val in sorted( plot_dit[tab].items(), key = lambda x:x[1], reverse=False):
            x.append( typ )
            y.append( val )
        #y.append( plot_dit[ typ ] )
    #pdit = dict(list(zip(x,y)))
    #x = rheMac.trick().sort(list(pdit.keys()))
    #label = x
    #labels = [i.replace('ternal','') for i in label ]
    #y = trick.dit(pdit).get(x)
    #label, y = sort_plot( label, y )
    x_int, label = [ i + 1 for i,v in enumerate( x ) ], x

    if bar == 'h' :
        fig = plt.figure( )
        ax = fig.add_subplot(111)
        ax.plot( x, y, 'bo')
        ax.vlines( x, 0, y, lw=2)
        ax.set_xlabel( tab[0] )
        ax.set_ylabel('-log10-Pvalue')
        plt.xticks( np.arange(len(x)) + 1, [i.replace('ternal','') for i in label], rotation = 90 )
        plt.axhline(y=1.3, color='r', linestyle='--')
        plt.show( )
        fig.savefig( pdf, bbox_inches='tight',pad_inches=+1, width = 12, height = 12 )
    elif bar == 'v':
        pval,typ = y,x
        typ.reverse()
        pval.reverse()
        x = np.arange(len(x))
        #fig = plt.figure(figsize = tuple(args.s))
        #formatter = FuncFormatter(pc)
        fig, ax = plt.subplots(1,2,figsize = tuple( (24, 120 )))
        fig.subplots_adjust(wspace=0, top=1, right=1, left=0, bottom=0.1)
        plt.barh(x, pval, alpha=1)
        ax[0].axis('off')
        #ax[0].xaxis.set_visible(False)
        #ax[0].xaxis.set_major_locator(plt.NullLocator())
        #labels = [i.replace('ternal','') for i in label ]
        plt.yticks(np.arange(len(x)), label)
        plt.axvline(x=1.3)
        fig.savefig( pdf, bbox_inches='tight',pad_inches=+1 )


def heatmap_plot( df ):
    #if args.order :
    #    with open( args.order ) as f :
    #        orders = [ i.replace('ME','') for i in f.readlines() ]
    #        print ( orders )
    df.to_csv( fix.fix(pdf).change(name + '.tab'), sep = '\t', index = True, index_label = 'module')
    #sns.clustermap( df, cmap='OrRd', method="ward", linewidths = 0.5, col_cluster = args.cc)
    df_corr = df.T.corr()
    DF_dism = 1 - df_corr
    linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
    #method weighted ward average single complete centroid
    sns.clustermap( df, method = "ward", linewidths = 0.5, col_cluster = args.cc, figsize=( 24, 20), cmap='OrRd' )
    #sns.clustermap( df, row_linkage=linkage, linewidths = 0.5, col_cluster = args.cc, figsize=( 24, 20), cmap='OrRd' )

if __name__ == '__main__':
    infor, kwargs, alias = parse( args )
    infor = label_deal( infor )
    source_dit = infor
    soure_only_num_dit = soure_only_num( infor )
    genes = input_symbols( alias )
    prefix = fix.fix( args.tab[0] )
    name = ''.join( random.sample( [ chr(i) for i in range(65,91) ], k = 4 ) )
    pdf = prefix.append( 'enrichment', name, 'pdf').replace(',','_')
    #if cat and len(cat) >= 2 :
    #    source_dit = trick.dit(source_dit).keyCat( cat )
    plot_dit = defaultdict( list )
    for fl in genes :
        plot_dit[fl] = prepare( genes[fl], infor, **kwargs)
    if len(genes) == 1 :
        #print ( plot_dit )
        #plot( list(plot_dit.values())[0], prefix, bar = args.pt)
        plot( plot_dit, prefix, bar = args.pt)
        if args.pt == 'heatmap':
            print ('Warning: heatmap only for multi not for single ...')
            exit()
    else :
        df = pd.DataFrame( plot_dit )
        df = df.fillna(0)
        print ( df )
        heatmap_plot( df )
        #sns.heatmap(df, annot=True, cmap='RdYlGn', linewidths=0.5)
        #sns.clustermap( df, cmap='OrRd', method="ward", linewidths = 0.5)
        #df.style.background_gradient(cmap='Blues')
    plt.savefig( pdf, dpi = 250, transparent = True, edgecolor = 'none' )
    print ( 'get {} postive'.format( df.shape, file = sys.stderr ) )
    for line in system.run('link_generate.py {}'.format(pdf), shell = True ):
        print ( line.strip() )

