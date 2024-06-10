#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 11:29:04 2019

@author: yyzou

this script to get the subtraction of two hic matrix

input : .cool file with a certain resolution

---------

"""


import matplotlib as mpl
mpl.use('Agg')
import random
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
plt.tight_layout()
import matplotlib.pyplot as plt
import numpy as np
import cooler
import click
import os
from os.path import exists
from math import log
from scipy.stats import zscore
import seaborn as sns;sns.set(color_codes=True)
from ningchao.nSys import fix

def read_in(file_path):
    return cooler.Cooler(file_path)

"""
  get matrix data then calculate input1 - input2

  normalized by log10() with positive and negative signs unchanged


"""

import math
import numpy as np
def get_Subtraction(c1,c2,area_s,area_e, chrom):
    area_s, area_e = int(area_s), int(area_e)
    print ( chrom, c1,c2,area_s,area_e)
    mat1 = c1.matrix(balance=False, sparse=True).fetch('{}:{}-{}'.format( chrom, area_s , area_e ))
    mat2 = c2.matrix(balance=False, sparse=True).fetch('{}:{}-{}'.format( chrom, area_s, area_e))
    #mat1 = (c1.matrix(balance=False, sparse=True)[area_s:area_e,area_s:area_e]).toarray()
    #mat2 = (c2.matrix(balance=False, sparse=True)[area_s:area_e,area_s:area_e]).toarray()
    mat1 = mat1.toarray()
    mat2 = mat2.toarray()
    x, y = mat1.shape
    s1, e1 = x - 25, y -18
    #s2, e2 = 80, 92
    s2, e2 = 80, y
    #np.random.randint(1, 5, size=(4, 4))
    #mat1[ 80:92, 80:92 ] = [[ random.randint( 10, 15 ) for i in range(20)] for j in range(20)]
    select_data = mat1[ x-e1: y-s1, s2 : e2]
    raw_data = mat1[ x-e1 - 6: y-s1 - 6, s2: e2]
    #mat1[ x-e1 - 6: y-s1 - 6, s2: e2] = np.random.randint( np.amin(select_data), np.amax(select_data), size=(e1-s1, e2-s2))
    mat1[ x-e1 - 6: y-s1 - 6, s2: e2] = select_data
    mat1[ x-e1: y-s1, s2 : e2] = raw_data
    #
    mat2[ x-e1: y-s1, s2 : e2] = mat2[ x-e1-6: y-s1-6, s2 : e2]
    #mat1[ x-e1: y-s1, s2 : e2] = np.random.randint( np.amin(select_data), np.amax(select_data)/6, size=(e1-s1, e2-s2))
    #mat1[ x-e:y-s, s:e] = np.random.randint(1, 30, size=(e-s, e-s))
    test = mat1[ x-e1:y-s1, s2:e2]
    arr1 = (mat1/mat1.sum())*(mat1.sum()+mat2.sum())
    arr2 = (mat2/mat2.sum())*(mat1.sum()+mat2.sum())
    #arr_sub = (np.log(arr1 + 1) + 1)/(np.log(arr2 + 1) + 1)
    arr_sub = arr1 - arr2
    arr_list = np.array([[x*log(x+1 if x>=0 else -x+1,10)/(abs(x)+1) for x in row ] for row in arr_sub])
    arr = arr_list.reshape(arr_list.shape[0],arr_list.shape[1])
    lst, labels = [ arr, np.log(arr1 + 1), np.log(arr2 + 1), test], [ '{}-{}'.format(c1.filename,c2.filename), c1.filename, c2.filename, 'test']
    return list( zip( lst, labels ) )


def get_fig( datas, area_s,area_e,outdir,outfig):
    figs_num = len( datas )
    fig = plt.figure(figsize=(10, 10 * figs_num))
    for i, data in enumerate( datas ):
        data, label = data
        ax = fig.add_subplot( figs_num, 1, i + 1 )
        #im = ax.matshow( zscore( data, axis =None), cmap='seismic',interpolation='none') ## vmin vmax determin the range of colorbar
        #im = ax.matshow( data, cmap='seismic',interpolation='none') ## vmin vmax determin the range of colorbar
        #im = ax.matshow( data, cmap='seismic',interpolation='none') ## vmin vmax determin the range of colorbar
        #sns.heatmap( data, cmap = 'bwr', robust = True, linewidths = 0, linecolor = None, ax = ax)
        center = 2.2 if i else 0.6
        sns.heatmap( data, cmap = 'seismic', robust = True, linewidths = 0, linecolor = None, ax = ax, center = center)
        #im = ax.matshow(data,cmap='YlOrRd',vmin=-1, vmax=1,interpolation='none')
        plt.ylabel( label )
        #fig.colorbar(im)
    pdf = fix.fix(outdir+outfig).append('pdf')
    fig.savefig( pdf )

def process(input1,input2,area_s,area_e,outdir,outfig, chrom):
    c1 = read_in(input1)
    c2 = read_in(input2)
    datas = get_Subtraction(c1,c2,area_s,area_e, chrom)
    get_fig( datas, area_s,area_e,outdir,outfig)


@click.command(name="hic_Subtraction")
@click.argument("input1")
@click.argument("input2")
@click.option("chrom","-c",
              #default=13,
              help="chrom for plot")
@click.option("area_s","-s",
              default=0,
              help="comparasion starts area")
@click.option("area_e","-e",
              default=None,
              help="comparasion starts area")
@click.option("--outdir", "-O",
    default="./",
    help="path to output files.")
@click.option("--outfig", "-fig",
    default="test.png",
    help="name of output figure.")

def main_(input1,input2, chrom, area_s,area_e, outdir, outfig):
    if not chrom :
        print ('!!!chrom must given')
    if not exists(outdir):
        os.mkdir(outdir)
    process(input1,input2, area_s,area_e,outdir,outfig, chrom)

if __name__ == "__main__":
    main_()
