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
from ningchao.nBio import chromosome,order
from tempfile import NamedTemporaryFile
from statsmodels.sandbox.stats.multicomp import multipletests

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'bed', nargs='?', help = 'bed file for plot')
parser.add_argument( '-rs', nargs='?', help = 'resolution', default = 40000 )
parser.add_argument( '-peirod', choices=['E50','E80','E90','E120','0M','4M','45Y','20Y'], help = 'peirods', nargs = '*' )
parser.add_argument( '-quiet','-q', action='store_true', help = 'quiet for the stderr' )
parser.add_argument( '-debug','-d', action='store_true', help = 'debug' )
parser.add_argument( '-adj', choices = ['bonferroni','fdr_bh','fdr_tsbh'], help = 'adj method', default = 'fdr_bh' )
parser.add_argument( '-dchr', nargs = '+', help = 'del which chr', default = ['chrX','chrY','chrMT','chrM'] )
parser.add_argument( '-includeReps', action = 'store_true', help = 'include the reps default not' )

if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()


def wdGetPeirod( wd ):
    p1,p2 = '',''
    dit = {'20Y':['20year','20Y'], '4M':['4month','4M'], 'E50': ['E50'], 'liver': ['liver'], 'E90':['E90']}
    dit.update({'E80':['E80'],'E120':['E120'],'0M':['0M'],'45Y':['4year']})
    reps = { 'rep1':['R1','rep1'],'rep2':['R2','rep2'] }
    for name in dit:
        allname = dit[name]
        for each in allname:
            if each in wd:
                p1 = name
    for name in reps:
        allname = reps[name]
        for each in allname:
            if each in wd:
                p2 = name
    if not p1:
        sys.stderr.write('check the name peirod %s rep %s\n' % (p1,p2))
        exit ( 1 )
    return '.'.join([p1,p2])


def dirs_lst(debug = False, peirods = [], inReps = False ):
    dirs = []
    #sun yao yu
    dirs.append('/dataE/rawdata/rheMac/hic/E80R1/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/dataE/rawdata/rheMac/hic/E80R2/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/dataE/rawdata/rheMac/hic/E80_mergerep/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/dataE/rawdata/rheMac/hic/E120R1/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/dataE/rawdata/rheMac/hic/E120R2/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/home/ningch/data/chip-seq/rheMac/rawdata/E120/prefrontal/hic/merge/Norm_matrix/ICE/40k')
    dirs.append('/dataE/rawdata/rheMac/hic/4monthR1/prefrontal_corte/Norm_matrix/ICE/')
    dirs.append('/dataE/rawdata/rheMac/hic/4monthR2/prefrontal_corte/Norm_matrix/ICE/' )
    dirs.append('/dataE/rawdata/rheMac/hic/4month_mergerep/prefrontal_corte/Norm_matrix/ICE/' )
    #wu shuai shuai results
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/20year/R1/SY-20YR1PLN1/7_R/3_Depth')
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/20year/R2/SY-RM21YR2PL/7_R/3_Depth/')
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/20year/Combine/PL/7_R/3_Depth')
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/4year/R1/SRM4YPLR1N3/7_R/3_Depth/')
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/4year/R2/SY-RM4YR2PL/7_R/3_Depth/')
    dirs.append('/dataE/rawdata/rheMac/hic/pfc/4year/Combine/PL/7_R/3_Depth')
    dirs.append('/dataE/rawdata/rheMac/hic/liver/R1/SY-LVR1N1/7_R/3_Depth')
    dirs.append('/dataE/rawdata/rheMac/hic/liver/R2/SY-LVR2/7_R/3_Depth/')
    dirs.append('/dataE/rawdata/rheMac/hic/liver/Combine/PL/7_R/3_Depth')
    #mine
    dirs.append('/dataB/ftp/pub/rheMac3/E50/prefrontal/hic/rep1/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/E50/prefrontal/hic/rep2/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/E50/prefrontal/hic/merge/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/E90/prefrontal/hic/rep1/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/E90/prefrontal/hic/rep2/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/E90/prefrontal/hic/merge/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/0M/prefrontal/hic/rep1/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/0M/prefrontal/hic/rep2/Norm_matrix/ICE/40k')
    dirs.append('/dataB/ftp/pub/rheMac3/0M/prefrontal/hic/merge/Norm_matrix/ICE/40k')
    p = re.compile(r'rep1|rep2|R1|R2')
    dirs = inReps and dirs or [ i for i in dirs if not p.search(i) ]
    if peirods:
        trans = { '20Y':'20year','45Y':'4year','4M':'4month' }
        peirods_trans = trick.lst( peirods ).ditReplace(trans)
        plot_peirods = '|'.join( peirods_trans )
        p1 = re.compile(r'%s' % plot_peirods)
        dirs = [ i for i in dirs if p1.search(i) ]
        out = []
        if 0 :
            for i,v in enumerate(dirs):
                if trick.string(v).fit(peirods, typ = 'or'):
                    out.append(v)
            return out
    if debug :
        dirs = dirs[0:2]
    return dirs

def dir2dit(dirs, load = []):
    dit = {}
    for each in dirs:
        each_arr = each.split('/')
        peirod = wdGetPeirod( each )
        files = os.listdir(each)
        fls = [ i for i in files if i.endswith('.ICE.scaleTo1.1w.txt') ]
        for fl in fls:
            chrom = fl.split('.')
            chrom = [ i for i in chrom if i.startswith('chr') ][0]
            abs_path = os.path.join(each,fl)
            trick.dinit(dit, peirod, 'work_dir', each)
            if chrom in load :
                sys.stderr.write('loading %s to matrix for extract data \n' % abs_path )
                arr = np.array(pd.read_csv(abs_path, header = None, sep = '\t'))
                trick.dinit(dit, peirod, 'chrom', chrom, arr)
            else :
                trick.dinit(dit, peirod, 'chrom', chrom, abs_path)
            
    return dit

def parse(fl, peirod ):
    dit = {}
    title = [ 'bin', 'chrom', 'shape', 'scale' ]
    fh = open(fl)
    for line in fh:
        wbin,chrom,shape,scale = line.strip().split('\t')
        # reverse the value
        wbin = abs(int(wbin) - 50 - 1 )
        trick.set3dict(dit, peirod, chrom, wbin, [ float(shape), float(scale) ])
    fh.close()
    return dit

def make_signal(dirs, quiet ):
    chroms,dit = list(chromosome.chr('rh8').schr),{}
    chroms = [ i for i in chroms if i != 'chrX' ]
    signal_file_dit = dir2dit(dirs)
    for peirod in signal_file_dit:
        for chrom in chroms:
            signal_file = signal_file_dit[peirod]['chrom'][chrom]
            if not quiet :
                sys.stderr.write('loading signal file %s...\n' % signal_file)
            work_dir = signal_file_dit[peirod]['work_dir']
            stat = Parse.ini(os.path.join(work_dir, 'readme.status.txt'),typ = 'tables').to_dict()
            ofl = '/'.join([work_dir, fix.fix( peirod ).append('%s.hic.signal' % chrom)])
            flag = '%s' % os.path.basename(ofl)
            trick.dinit( dit, peirod, chrom, ofl)
            if flag not in stat or not os.path.exists(ofl) or os.stat(ofl).st_size == 0 or line_num.line_num(ofl) != 50 :
                ofh = open(ofl,'w')
                sys.stderr.write('write the signal in %s\n' % ofl)
                array = np.array(pd.read_csv(signal_file, header = None, sep = '\t'))
                for i in range(50):
                    diagonal_i = np.take(array, np.arange(len(array)*(i+1), array.size, len(array) + 1))
                    selected = (( diagonal_i != 0 ) & (diagonal_i < np.percentile(diagonal_i, 95) ))
                    diagonal_i = diagonal_i[ selected ]
                    floc, shape, f0, scale = stats.exponweib.fit(diagonal_i, floc=0, f0=1)
                    ofh.write('\t'.join([str(i+1), chrom, str(shape), str(scale)]) + '\n')
                ofh.close()
                status.cmd(work_dir, flag).mark()
                if not quiet:
                    sys.stderr.write('write the signal in %s done...\n' % ofl)
            else :
                if not quiet :
                    sys.stderr.write('!Already done and ignore the %s\n' % ofl)
    return dit

def hic_dit(chroms, quiet):
    dit = {}
    signal_file_dit = make_signal( dirs, quiet)
    for peirod in signal_file_dit:
        for chrom in chroms:
            signal_file = signal_file_dit[peirod][ chrom ]
            if not quiet:
                sys.stderr.write( 'loading signal the file %s\n' % signal_file )
            #Too slow 
            #array = np.loadtxt(signal_file) 
            array = np.array(pd.read_csv(signal_file, header = None, sep = '\t'))
            for line in array:
                chrom = line[1]
                sbin = line[0]
                value = line[2:]
                trick.dinit(dit, peirod, chrom, sbin, value)
    return dit
#stats.weibull_min.logsf(95.8634189116823, 2.17771191111,loc=0, scale = 38.1917534204)


def bline_pick_pvalue( bline, signal, rs, array, genome):
    bfh =  open(bed)
    peirods, postive_bins = [],[]
    if bline.strip():
        values, pvalues = {}, {}
        coord = b2span(line, rs, genome)
        chrom, start, end, name, place_hold, strand = bline.strip().split('\t')
        for i,( rv, lv) in enumerate(coord):
            for peirod in array:
                if peirod not in peirods:
                    peirods.append(peirod)
                value = array[peirod]['chrom'][chrom][np.ix_([ rv ], [ lv ])][0][0]
                try :
                    value = array[peirod]['chrom'][chrom][np.ix_([ rv ], [ lv ])][0][0]
                except IndexError as err:
                    print(array[peirod]['chrom'][chrom],rv,lv)
                    print(peirod, chrom, 'not exists')
                    print(err)
                    exit( 1 )
                except :
                    print("Unexpected error:", sys.exc_info()[0])
                    raise
                trick.dinit( values, i, peirod, value )
        for pos in values:
            value = values[pos]
            bg_pos = pos + 1
            if bg_pos < 51 :
                bg_pos = 51 - bg_pos
            elif bg_pos > 51 :
                bg_pos = bg_pos - 51
            else :
                bg_pos = 1
            for peirod in peirods:
                shape, scale = signal[peirod][chrom][bg_pos]
                trick.dinit(pvalues, peirod, [])
                x = values[pos][peirod]
                if pos == 50 :
                    pvalue = 0
                else :
                    pvalue =  - stats.weibull_min.logsf( x, shape, loc=0, scale = scale)
                    print ( pvalue )
                #print(shape,scale,pvalue,bg_pos)
                #if pvalue > 30 :
                #    pvalue = 30.00
                #pvalues[peirod].append((pvalue, x, shape, scale))
                pvalues[peirod].append(round(pvalue, 4))
        if adj:
            for peirod in peirods:
                pval = [ math.pow(10, -x) for x in pvalues[peirod] ]
                '''
`bonferroni` : one-step correction
`sidak` : one-step correction
`holm-sidak` : step down method using Sidak adjustments
`holm` : step-down method using Bonferroni adjustments
`simes-hochberg` : step-up method  (independent)
`hommel` : closed method based on Simes tests (non-negative)
`fdr_bh` : Benjamini/Hochberg  (non-negative)
`fdr_by` : Benjamini/Yekutieli (negative)
`fdr_tsbh` : two stage fdr correction (non-negative)
`fdr_tsbky` : two stage fdr correction (non-negative)
'''
                qval = list(multipletests( pval, method = adj )[1])
                try :
                    pval = [ round(-math.log(x), 4 ) for x in qval  ]
                except :
                    sys.stderr.write('%s can not be use\n' % bline)
                pvalues[peirod] = pval
    return coord,pvalues, postive_bins

def b2span(bline, rs, genome):
    chrom, start, end, name, place_hold, strand = bline.strip().split('\t')
    max_bin = genome[chrom] // rs
    if strand == '+':
        sbin = int(start) // rs
    else :
        sbin = int( end ) // rs
    if sbin > max_bin:
        sys.stderr.write( 'promoter big than max_bin for line %s, max bin is %s\n' % (bline,sbin) )
        return []
    lst = list(range(-50, 0, 1))
    lst.extend(list(range(0, 51, 1)))
    pbin = [ i + sbin for i in lst ]
    span = [ i for i in pbin if i >= 0  and i < max_bin ]
    coord = [ [ sbin , i  ] for i in span ]
    return coord

def bedGraph( chrom, coord, values, rs, filename ):
    fh = open( filename, 'w')
    dit = {'chr':[],'start':[],'end':[]}
    for i,value in enumerate(coord):
        s,e = value
        #print i, s, e
        dit['chr'].append(chrom)
        dit['start'].append( str(e  * rs))
        dit['end'].append(str((e+1) * rs))
    dit.update(values)
    df = pd.DataFrame(dit)
    df = df.reindex(columns = order.orderLst(df.columns.values.tolist()))
    df.to_csv( filename, index = False, sep = '\t')
def need_chroms(beda):
    chroms = []
    fh = open(bed)
    for line in fh:
        line_arr = line.strip().split('\t')
        if dchr and line_arr[0] in dchr:
            continue
        if line_arr[0] not in chroms:
            chroms.append(line_arr[0])
    fh.close()
    return chroms


def already():
    fls = os.listdir('.')
    tabs = [ i for i in fls if i.endswith('tab') ]
    infor = {} 
    for tab in tabs:
        gene = tab.replace('.tab','')
        trick.dinit(infor, gene, [])
        lst = open(tab).readlines()
        if len(lst) == 102 :
            tregion = (lst[ 51 ].strip().split('\t')[0:3])
            infor[ gene ].append(tregion)
            sys.stderr.write('!Already for %s,%s\n' % (gene, ','.join(tregion)))
    return infor

def calready( chrom, start, end, name, ainfor ):
    if name in ainfor:
        for each in ainfor[name]:
            achrom, astart, aend = each
            astart,aend,start,end = [ int(i) for i in (astart,aend,start,end) ]
            if chrom == achrom and astart < start < aend  and astart < end < aend:
                return True 
    return False 






if __name__ == '__main__':
    rs, bed, quiet, inReps = args.rs, args.bed, args.quiet, args.includeReps
    adj,debug,dchr,speirods = args.adj,args.debug, args.dchr, args.peirod
    ainfor,already_lines, needLoadChroms = already(), [], []
    fh, dirs = open(bed), dirs_lst( debug, speirods, inReps)
    tfh = NamedTemporaryFile('w+t', delete = False, dir = '.')
    for line in fh:
        chrom, start, end, name, place_hold, strand = line.strip().split('\t')
        if calready( chrom, start, end, name, ainfor ):
            sys.stderr.write('!Already for %s' % line)
            already_lines.append(line)
        else :
            tfh.write(line)
            needLoadChroms.append( chrom )
    tfh.close()
    assert os.path.getsize(tfh.name)
    #make signal file for speed
    dit = make_signal(dirs, quiet)
    chroms = need_chroms(bed)
    signal = hic_dit( chroms, quiet)
    raw = dir2dit(dirs, load = needLoadChroms )
    genome = chromosome.chr('rh8').chr
    fh = open(tfh.name)
    for i,line in enumerate(fh):
        chrom, start, end, name, place_hold, strand = line.strip().split('\t')
        if chrom in dchr:
            continue
        coord, pvalues, postive_bins = bline_pick_pvalue( line, signal, rs, raw, genome)
        if coord == []:
            continue
        #write the result in the bedGraph
        tab = '%s.tab' % ('.'.join([name, chrom, start, end, strand]))
        pdf = tab.replace('.tab','.pdf')
        bedGraph( chrom, coord, pvalues, rs, tab )
        df = pd.DataFrame.from_dict(pvalues)
        peirods = list( df )
        #df = dataframe.df(df).insert(49, [ 0 for i in peirods ])
        ax = df.plot(style='.-')
        #ax.set_ylim(0, 10)
        ax.axvline(x= 50, linestyle = '-.')
        ax.axhline(y=3, linestyle = '-.')
        plt.savefig( pdf, format='pdf')
        plt.close()


















