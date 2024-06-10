#!/usr/bin/env python3
import os
import sys
import argparse
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.tight_layout()
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex, colorConverter
from collections import defaultdict
import ningchao.nSys.fix as fixKit
from pandas import read_csv
import numpy as np
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
import seaborn as sns;sns.set(color_codes=True)
import ningchao.nSys.trick as trKit
import matplotlib.patches as mpatches
from ningchao.nSys import excel,fix,parse,convert
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster

sys.setrecursionlimit(1000000)
sns.set(context="paper", font="monospace")
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='print useage for script', formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('i', nargs='?', help ='matrix for plot')
parser.add_argument('-p','-outputPrefix', nargs='?', help ='output prefix')
parser.add_argument('-c','-center', nargs='?', help ='center for the color. float',type = float )
parser.add_argument('-vmax', nargs='?', help ='vmax. float', type = float )
parser.add_argument('-vmin', nargs='?', help ='vmin. float', type = float )
parser.add_argument('-k','-key', nargs='*', help ='key for heatmap, can think is the ignor cols', default = [1])
parser.add_argument('-f','-fiter', nargs='*', help ='fiter cols for the heatmap', type = int )
parser.add_argument('-sw', action='store_true', help ='swith for look the header order infor')
parser.add_argument('-log2', action='store_true', help ='+1log2 transfrom')
parser.add_argument('-2n', action='store_true', help ='2**val')
parser.add_argument('-std', choices = [ 0, 1, None ], type = int, help ='standard_scale. 0 for rows, 1 for columns. default none', default = None )
parser.add_argument('-z_score','-z',choices = [ '0', '1', 'None' ], help ='0 for rows, 1 for columns.', default = '0')
parser.add_argument('-m', choices = ['single','average','complete','ward','single','mcquitty','median','centroid'], help ='output prefix' )
parser.add_argument('-cluster_num', nargs='?', help ='cluster_num for input matrix, default is length of columns. input num is add to the whole', type = int )
parser.add_argument('-cmap', nargs='?', help ='cmap define Reds|OrRd is good' )
#parser.add_argument('-cmap', choices = ['Accent','Accent_r','Blues','Blues_r','BrBG','BrBG_r','BuGn','BuGn_r','BuPu','BuPu_r','CMRmap','CMRmap_r','Dark2','Dark2_r','GnBu','GnBu_r','Greens','Greens_r','Greys','Greys_r','OrRd','OrRd_r','Oranges','Oranges_r','PRGn','PRGn_r','Paired','Paired_r','Pastel1','Pastel1_r','Pastel2','Pastel2_r','PiYG','PiYG_r','PuBu','PuBuGn','PuBuGn_r','PuBu_r','PuOr','PuOr_r','PuRd','PuRd_r','Purples','Purples_r','RdBu','RdBu_r','RdGy','RdGy_r','RdPu','RdPu_r','RdYlBu','RdYlBu_r','RdYlGn','RdYlGn_r','Reds','Reds_r','Set1','Set1_r','Set2','Set2_r','Set3','Set3_r','Spectral','Spectral_r','Vega10','Vega10_r','Vega20','Vega20_r','Vega20b','Vega20b_r','Vega20c','Vega20c_r','Wistia','Wistia_r','YlGn','YlGnBu','YlGnBu_r','YlGn_r','YlOrBr','YlOrBr_r','YlOrRd','YlOrRd_r','afmhot','afmhot_r','autumn','autumn_r','binary','binary_r','bone','bone_r','brg','brg_r','bwr','bwr_r','cool','cool_r','coolwarm','coolwarm_r','copper','copper_r','cubehelix','cubehelix_r','flag','flag_r','gist_earth','gist_earth_r','gist_gray','gist_gray_r','gist_heat','gist_heat_r','gist_ncar','gist_ncar_r','gist_rainbow','gist_rainbow_r','gist_stern','gist_stern_r','gist_yarg','gist_yarg_r','gnuplot','gnuplot2','gnuplot2_r','gnuplot_r','gray','gray_r','hot','hot_r','hsv','hsv_r','icefire','icefire_r','inferno','inferno_r','jet','jet_r','magma','magma_r','mako','mako_r','nipy_spectral','nipy_spectral_r','ocean','ocean_r','pink','pink_r','plasma','plasma_r','prism','prism_r','rainbow','rainbow_r','rocket','rocket_r','seismic','seismic_r','spectral','spectral_r','spring','spring_r','summer','summer_r','tab10','tab10_r','tab20','tab20_r','tab20b','tab20b_r','tab20c','tab20c_r','terrain','terrain_r','viridis','viridis_r','vlag','vlag_r','winter','winter_r'], help ='default. summer_r' )
parser.add_argument('-cl','-cluster', choices = ['none','row','col','both'],help = 'cluster', default = 'none')
parser.add_argument('-dist', choices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'], help = 'dist method for col')
parser.add_argument('-yColorFile', '-y', help = 'color lst file')
parser.add_argument('-ylabel', action = 'store_true', help = 'ylabel or not', default = False)
parser.add_argument('-split', help = 'split the color', action='store_true')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

sns.set_palette('Set1', 10, 0.65)
palette = sns.color_palette()
set_link_color_palette(list(map(rgb2hex, palette)))
sns.set_style('white')
vmax = args.vmax
vmin = args.vmin
fiter = False
if args.f:
    fiter = [i -1 for i in args.f ]
yColorFile = args.yColorFile

#def and classes
class Clusters(dict):
    def _repr_html_(self):
        html = '<table style="border: 0;">'
        for c in self:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
            '<td style="background-color: {0}; ' \
                       'border: 0;">' \
            '<code style="background-color: {0};">'.format(hx)
            html += c + '</code></td>'
            html += '<td style="border: 0"><code>'
            html += repr(self[c]) + '</code>'
            html += '</td></tr>'

        html += '</table>'

        return html

def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = Clusters()
    for c, l in list(cluster_idxs.items()):
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes

#####



fl = args.i
sw = args.sw 
if sw:
    print(excel.xls(fl).lstinfor())
    exit()
method = args.m
center = args.c
metric = out = ''
if args.dist:
    metric = args.dist
col_cluster = False
row_cluster = False
if args.cl == 'row':
    row_cluster = True
elif args.cl == 'col':
    col_cluster = True
elif args.cl == 'both':
    row_cluster = col_cluster = True
elif args.cl == 'none':
    col_cluster = False
    row_cluster = False
else :
    exit('wrong input cluster\n')


if args.p:
    prefix = args.p
else :
    prefix = fl
if prefix.endswith('pdf'):
    out = prefix
else :
    out = prefix +'.pdf'

def add(key,val,num = 0):
    out = ',' + '%s=' % key + '\'' + str(val) + '\''
    if num:
        return out.replace('\'','')
    else :
        return out

index = args.k
index = [ int(i) - 1 for i in index ]
print('file %s' % fl)
print('index cols is: ', end=' ')
print(index)


### before read in add change the index name to change the color
## name and type in
ckcols = [1,2,3]
cvcols = [5,6]
if yColorFile:
    fl_fh = open(fl)
    tmp_fl = fl + '.tmp'
    tmp_fh = open( tmp_fl, 'w')
    #postion = parse.parse(yColorFile).tab2dict(kcols = [1,2,3], vcols = [5,6], other = 6)
    postion = parse.parse(yColorFile).tab2dict(kcols = ckcols, vcols = cvcols)
    for line in fl_fh:
        line_arr = line.strip().split('\t')
        key_lst = [ line_arr[i - 1] for i in ckcols ]
        key = '\t'.join(key_lst)
        if key in postion:
            line_arr[0] = postion[key]
            line_arr[1] = '0'
            line_arr[2] = '1'
        tmp_fh.write('\t'.join(line_arr) + '\n' )
    tmp_fh.close()
    x = read_csv(tmp_fl,sep = '\t', header = 0, index_col = index)
    if fiter :
        length = len(index) + len(x.columns)
        use_cols = list(range(length))
        use_cols = [ i for i in use_cols if i not in fiter ]
        x = read_csv(tmp_fl,sep = '\t', header = 0, index_col = index, usecols = use_cols)
else :
    yTick = False
    x = read_csv(fl,sep = '\t', header = 0, index_col = index)
    if fiter :
        length = len(index) + len(x.columns)
        use_cols = list(range(length))
        use_cols = [ i for i in use_cols if i not in fiter ]
        x = read_csv(fl,sep = '\t', header = 0, index_col = index, usecols = use_cols)
print ( x )
x = x.dropna()
if args.log2:
    x = x + 1
    x = np.log2(x)
    print ( x )
header = x.columns

cmap = args.cmap and args.cmap or None
split = args.split
#x, row_labels, col_labels = make_data(fl)
#sns.heatmap(flights, annot=True, fmt="d")
#cmap = sns.diverging_palette(220, 20, s=85, l=25, n = 7, as_cmap=True)
#pal = sns.diverging_palette(240, 60, s= 99, l = 80, n = 100, center="dark", as_cmap = True )
pal = sns.diverging_palette( 194, 60, s = 60, l = 85, n = 100, center="dark", as_cmap = True )
if yColorFile:
    yTick = True
    if 1 :
        yTickLst = []
        yTick = x.index.values
        #for color
        dit = parse.parse(yColorFile).tab2dict(kcols = ckcols, vcols = [5])
        for i,each in enumerate(yTick):
            key = '\t'.join([str(i) for i in each ])
            if key in dit:
                #name = [key]
                name = []
                name.extend(dit[key])
                yTickLst.append(','.join(name))
            else :
                yTickLst.append('\t'.join([str(j) for j in each] ))
        yTick = yTickLst
else :
    yTick = False 

#f, ax = plt.subplots(figsize=(20, 8))
#mcmap = ListedColormap(sns.color_palette(pal).as_hex()) center=x.loc["January", 1955])
z_score = None if args.z_score == 'None' else int(args.z_score)
h = '''sns.clustermap( x, col_cluster=col_cluster,row_cluster=row_cluster,standard_scale=args.std,z_score = z_score, robust=True, linewidths = 0, linecolor = None '''
#h += add('cmap', pal, num = 1)
if method :
    h += add('method',method)
if metric :
    h += add('metric',metric)
if center :
    h += add('center',center, num = 1)
if vmax:
    h += add('vmax', vmax, num = 1)
if vmin:
    h += add('vmin', vmin, num = 1)
if not cmap:
    h += ',cmap = pal'
else :
    h += ''',cmap = '{}' '''.format(cmap)

if yColorFile:
    h += ',yticklabels' + '=' + 'yTick'
else :
    h += add('yticklabels', args.ylabel , num = 1 )
h += ')'
print(h)
h = eval(h)
#h.ax_heatmap.set_yticklabels( , rotation = 0)

print('eval done')

# define the color for the specify row
def change_key(dit):
    out = {}
    for key in dit:
        value = dit[key]
        key = key.replace('\t',',')
        key = '\t'.join([key, key.split(',')[0]])
        out.update({key : value})
    return out
m = 0
if yColorFile and 1 :
    #dit = parse.parse(args.yColorFile).tab2dict(kcols = [5,6], vcols = [4], vsep = ',')
    dit = parse.parse(yColorFile).tab2dict(kcols = cvcols, vcols = [4])
    dit = change_key(dit)
    print(dit) 
    #dit2 = parse.parse(args.yColorFile).tab2dict(kcols = [4,6], vcols = [4], vsep = ',', ksep=',')
    for tick_label in h.ax_heatmap.axes.get_yticklabels():
        tick_text = tick_label.get_text()
        key = tick_text.replace('-','\t')
        key = key.replace('\t0\t1','')
        key = convert.convert(key)
        if key in dit:
            print('set color for %s' % key)
            tick_label.set_color([ float(i) for i in dit[key].split(',')])
        else :
            tick_label.set_alpha(0)
        #species_name = species.loc[int(tick_text)]




#h.ax_cbar.set_position((0.8, .2, .03, .4))
prefix = fix.fix(out)
norfh = open(prefix.append('normal.xls').get(), 'w')
odfh = open(prefix.append('cluster.order.xls').get(), 'w')
odfh2 = open(prefix.append('cluster2.xls').get(), 'w')
print ('geneId', *header, sep = '\t', file = norfh)
print ('geneId', *header, sep = '\t', file = odfh)
oldorder = list(x.index)
oldheader = list(x.axes[1])
#legend_TN = [mpatches.Patch(color=c, label=l) for c,l in df[['tissue type','label']].drop_duplicates().values]
#l2=g.ax_heatmap.legend(loc='center left',bbox_to_anchor=(1.01,0.85),handles=legend_TN,frameon=True)
#l2.set_title(title='tissue type',prop={'size':10})
if row_cluster:
#   print dir(h.dendrogram_row.linkage)
#   den = dendrogram(h.dendrogram_row.linkage , labels=x.index, above_threshold_color='#AAAAAA', leaf_rotation=0, orientation="left", color_threshold=240)
    #import pylab
    #fig = pylab.figure(figsize=(20,8))
    #ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    plt.xticks(rotation=45)
    h.savefig( fix.fix(out).insert('raw'),format='pdf')
    den = dendrogram( h.dendrogram_row.linkage, labels=x.index, above_threshold_color='#AAAAAA', leaf_rotation=0, orientation="left",no_labels = True )
#   print h.dendrogram_row._calculate_linkage_fastcluster()

#new method for cluster extract 
    L = h.dendrogram_row.linkage
    if args.cluster_num:
        cluster_index = fcluster(L, len(x.dtypes) + args.cluster_num, criterion='maxclust')
    else :
        cluster_index = fcluster(L,len(x.dtypes),criterion='maxclust')

    x['cluster'] = cluster_index
    #x = x.sort(['cluster'] ,ascending=[True])
    x = x.sort_values( by = 'cluster', ascending=[True])
    x.to_csv(prefix.append('cluster3.xls'),sep='\t')


#
    h.ax_row_dendrogram = den
    clusters = get_cluster_classes( h.ax_row_dendrogram )
    for clu in clusters:
        print(clu +' length is: %d' %  len(clusters[clu]))
        keys = clusters[clu]
        for i,key in enumerate( keys ):
            if isinstance( key, list) :
                key = [str(k) for k in list(key)]
                keys[i] = '\t'.join(key)
        odfh2.write(clu + '\t' + ','.join([str(i) for i in keys])+'\n')
    for i,each in enumerate(h.dendrogram_row.reordered_ind):
        line = []
        key = oldorder[each]
        if isinstance( key, list ):
            key = [str(k) for k in key]
            key = '\t'.join(key)
        line.append(key)
        line.extend([ str(j) for j in h.dendrogram_row.array[i]])
#       print dir(h.dendrogram_row)
#       den = dendrogram(h.dendrogram_row.linkage, labels=x.index, above_threshold_color='#AAAAAA')
        odline = [ str(key),'\t'.join([ str(k) for k in x.loc[oldorder[each]].values ])]
        odfh.write('\t'.join(odline) + '\n')
        norfh.write('\t'.join([ str(i) for i in line]) + '\n')

if col_cluster:
    print([ oldheader[m] for m in h.dendrogram_col.reordered_ind ])

#for _, spine in h.spines.items():
    #spine.set_visible(False)
#h.ax_heatmap.grid( visible = False, which='major')
#h.ax_heatmap.grid( visible = False, which='minor')
#h.ax_heatmap.set_xticks( range(len(header)), rotation=45)
h.ax_heatmap.set_xticklabels( header, rotation=40 )
h.savefig(out,format='pdf')
os.system('link_generate.py {}'.format(out))
