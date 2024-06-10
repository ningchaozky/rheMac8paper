#!/usr/bin/env python3
import os
import sys
import time
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import argparse
import pandas as pd
from ningchao.nSys import trick, OONMF, OONMFhelpers, fix
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
plt.tight_layout()
#plt.tick_params( axis = 'both', left = True, labelleft = True, which = 'both', bottom = True, top = False, labelbottom = True, direction = 'in' )
#sns.despine( ax = ax[2]  )#remove top and right axis
#subplot define, also polar plot
#with sns.axes_style("whitegrid", {'xtick.top': True, 'ytick.left': True, 'axes.grid': False, 'ytick.color': 'black', 'ytick.direction': 'in' }):
    #fig, ax = plt.subplots( 6, figsize=(10,40), subplot_kw=dict(projection=None))
            #fig,ax = plt.subplots( 3, figsize=(10,30), sharex=True, sharey=True, subplot_kw=dict(projection='polar')); ax[0],ax[1]
            #plt.style.use('ggplot')
            #plt.subplots_adjust(wspace=0, hspace=0)
import seaborn as sns;sns.set(color_codes=True)
sns.set_style("ticks")
sns.barplot( palette="Set3" )
            #fig = plt.figure( figsize=( 18, 24) )
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('bins', nargs='?', help ='bins for cluster')
parser.add_argument('-num', nargs='?', help ='nums you want to cluster', required = True, type = int )
parser.add_argument('-p', nargs='?', help ='prefix for output')
parser.add_argument('-newone', nargs='?', help ='new one, or keep same', default = 'keep')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


if 1 :
    prefix = args.p if args.p else args.bins
    df = pd.read_table( args.bins , header = 1, index_col = 1)
    A = df.T
    #A = A.fillna(2, inplace = True)
    #A = A.fillinfinite(2, inplace = True)
    t = np.any(np.isnan(A))
    print ('NAN: {}'.format(t), file = sys.stderr)
    t = np.all(np.isfinite(A))
    print ('Infinite: {}'.format(t), file = sys.stderr)
    a = OONMF.NMFobject(theNcomps= args.num)
    #a = OONMF.NMFobject()
    #a.performNMF(data=A, randomseed=seed, theinit='nndsvd')
    if args.newone :
        a.performNMF(data=A, theinit='random')
    else :
        seed = 200
        a.performNMF(data=A, randomseed=seed, theinit='random')
    Basis, Mixture = fix.fix( prefix ).append( 'dim'+str(args.num), 'Basis.npy'), fix.fix( prefix ).append( 'dim'+str(args.num), 'Mixture.npy')
    a.writeNMF(Basis_foutname= Basis, Mixture_foutname = Mixture)
    Basis_tsv, Mixture_tsv = fix.fix( Basis ).change( 'tsv' ), fix.fix( Mixture ).change( 'tsv' )
    print ( Basis_tsv )
    a.writeNMF_CSV(Basis_foutname= Basis_tsv, Mixture_foutname = Mixture_tsv, sep = '\t', index_label = 'pos')

if 0 :
    df = pd.read_table('NNDSVD_Basis.tsv', header = 1 )
    #df = pd.read_table('NNDSVD_Mixture.tsv', header = 1 )
    time_start = time.time()
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(df)
    print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start), file = sys.stderr)
    df['tsne-x'] = tsne_results[:,0]
    df['tsne-y'] = tsne_results[:,1]
    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-x", y="tsne-y",
        #hue="y",
        palette=sns.color_palette("hls", 10),
        data=df,
        legend="full",
        alpha=0.3
    )

    plt.savefig( 'test.pdf', dpi=250, transparent = True)
    exit()


























