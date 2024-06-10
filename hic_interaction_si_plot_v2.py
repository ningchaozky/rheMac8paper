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
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from ningchao.nSys import trick,fix,dataframe,status,line_num, system
from ningchao.nSys import parse as Parse
from ningchao.nBio import chromosome,order
from tempfile import NamedTemporaryFile
from statsmodels.sandbox.stats.multicomp import multipletests
from collections import defaultdict

example = '''if is to many you should put only one peirod'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'bed', nargs='?', help = 'bed file for plot')
parser.add_argument( '-peirod','-p', choices=['E50','E90','E120','0M','4M','45Y','20Y'], help = 'peirods', nargs = '*', required = True )
parser.add_argument( '-rs', nargs = '?', help = '10000|40000', default = 10000, type = int  )
parser.add_argument( '-span', nargs = '?', help = '6000000', default = 6000000, type = int  )
parser.add_argument( '-quiet','-q', action='store_true', help = 'quiet for the stderr' )
parser.add_argument( '-debug','-d', action='store_true', help = 'debug' )
parser.add_argument( '-adj', choices = ['bonferroni','fdr_bh','fdr_tsbh'], help = 'adj method', default = 'fdr_bh' )
parser.add_argument( '-dchr', nargs = '+', help = 'del which chr', default = ['chrX','chrY','chrMT','chrM'] )
parser.add_argument( '-includeReps', action = 'store_true', help = 'include the reps default not' )

if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()



def parse(fl, peirod, **kwargs):
    dit = {}
    title = [ 'bin', 'chrom', 'shape', 'scale' ]
    fh = open(fl)
    for line in fh:
        wbin,chrom,shape,scale = line.strip().split('\t')
        # reverse the value
        wbin = abs(int(wbin) - kwargs.get('span_bins') - 1 )
        trick.set3dict(dit, peirod, chrom, wbin, [ float(shape), float(scale) ])
    fh.close()
    return dit

def make_signal( **kwargs ):
    infor = kwargs.get('ice_fls')
    dit = defaultdict( lambda : defaultdict( str ) )
    for peirod in infor :
        for chrom in infor[peirod] :
            signal_file = infor[peirod][ chrom ]
            if not kwargs.get('quiet') :
                sys.stderr.write('loading distribution parameters file: {}, chrom: {}...\n\n'.format( signal_file, chrom ))
            work_dir = os.path.dirname( signal_file )
            stat = Parse.ini(os.path.join(work_dir, 'readme.status.txt') ).to_dict()
            if kwargs.get('span_bins') == 300 :
                ofl = fix.fix( signal_file ).append('{}.hic.signal'.format(chrom))
            else :
                ofl = fix.fix( signal_file ).append('{}.{}.hic.signal'.format(chrom, kwargs.get('span_bins')))
            flag = '%s' % os.path.basename(ofl)
            dit[peirod][chrom] = ofl
            #if flag not in stat or not os.path.exists(ofl) or os.stat(ofl).st_size == 0 or line_num.line_num(ofl) != 300 :
            if flag not in stat :
                ofh = open(ofl,'w')
                sys.stderr.write('write the signal in %s\n' % ofl)
                array = np.array( pd.read_csv(signal_file, header = None, sep = '\t') )
                for i in range( kwargs.get('span_bins') ):
                    diagonal_i = np.take(array, np.arange(len(array)*(i+1), array.size, len(array) + 1))
                    selected = (( diagonal_i != 0 ) & (diagonal_i < np.percentile(diagonal_i, 95) ))
                    diagonal_i = diagonal_i[ selected ]
                    floc, shape, f0, scale = stats.exponweib.fit(diagonal_i, floc=0, f0=1)
                    oarr = [str(i+1), chrom, shape, scale]
                    print ( *oarr, sep = '\t', file = ofh )
                ofh.close()
                status.cmd( flag = flag, wd = work_dir).mark()
                if not kwargs.get('quiet'):
                    sys.stderr.write('write the signal in %s done...\n' % ofl)
            else :
                if not kwargs.get('quiet') :
                    sys.stderr.write('!Already done and ignore the %s\n' % ofl)
    return dit

def hic_dit( **kwargs ):
    quiet, make_signal_dit = [ kwargs.get(i) for i in ('quiet','make_signal_dit')]
    dit = defaultdict( lambda : defaultdict( lambda : defaultdict( float )) )
    #signal_file_dit = make_signal( dirs, quiet)
    for peirod in make_signal_dit:
        for chrom in kwargs.get('need_chroms'):
            signal_file = make_signal_dit[peirod][ chrom ]
            if not quiet:
                print( 'loading distribution parameters the file: ', signal_file, sep = ' ', file = sys.stderr)
            #Too slow 
            #array = np.loadtxt(signal_file) 
            try :
                array = np.array(pd.read_csv(signal_file, header = None, sep = '\t'))
            except :
                print ( 'Check', signal_file, chrom, peirod, 'maybe not exists...',file = sys.stderr )
                exit()
            for line in array:
                sbin,chrom = line[0:2]
                dit[peirod][chrom][sbin] = line[2:]
                #trick.set3dict(dit, peirod, chrom, sbin, value)
    return dit
#stats.weibull_min.logsf(95.8634189116823, 2.17771191111,loc=0, scale = 38.1917534204)


def bline_pick_pvalue( bline, **kwargs ):
    signal, rs, array, genome = [ kwargs.get(i) for i in ('signal','rs','array','genome') ]
    #bfh =  open( kwargs.get('bed') )
    peirods, postive_bins = [],[]
    if bline.strip():
        values, pvalues = defaultdict( lambda : defaultdict( float ) ), defaultdict( list  )
        coord = b2span( line, **kwargs)
        chrom, start, end, name, place_hold, strand = bline.strip().split('\t')
        for i,( rv, lv) in enumerate(coord):
            for peirod in array:
                if peirod not in peirods:
                    peirods.append(peirod)
                value = array[peirod][chrom][np.ix_([ rv ], [ lv ])][0][0]
                try :
                    value = array[peirod][chrom][np.ix_([ rv ], [ lv ])][0][0]
                except IndexError as err:
                    print(array[peirod][chrom],rv,lv)
                    print(peirod, chrom, 'not exists')
                    print(err)
                    exit( 1 )
                except :
                    print("Unexpected error:", sys.exc_info()[0])
                    raise
                values[i][peirod] = value
                #trick.set2dict( values, i, peirod, value )
        for pos in values:
            value = values[pos]
            bg_pos = pos + 1
            if bg_pos < kwargs.get('span_bins') + 1:
                bg_pos = kwargs.get('span_bins') + 1 - bg_pos
            elif bg_pos > kwargs.get('span_bins') + 1 :
                bg_pos = bg_pos - ( kwargs.get('span_bins') + 1 )
            else :
                bg_pos = 1
            for peirod in peirods:
                shape, scale = signal[peirod][chrom][bg_pos]
                #trick.set1dict(pvalues, peirod, [])
                x = values[pos][peirod]
                if pos == kwargs.get('span_bins'):
                    pvalue = 0
                else :
                    pvalue =  - stats.weibull_min.logsf( x, shape, loc=0, scale = scale)
                #print(shape,scale,pvalue,bg_pos)
                #if pvalue > 30 :
                #    pvalue = 30.00
                #pvalues[peirod].append((pvalue, x, shape, scale))
                pvalues[peirod].append(round(pvalue, 4))
        if kwargs.get('adj'):
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
                qval = list(multipletests( pval, method = kwargs.get('adj') )[1])
                try :
                    pval = [ round(-math.log(x), 4 ) for x in qval  ]
                except :
                    print('{} can not be use'.format(bline), file = sys.stderr)
                pvalues[peirod] = pval
    return coord, pvalues, postive_bins

def b2span(bline, **kwargs):
    rs, genome = [ kwargs.get(i) for i in ('rs','genome') ]
    chrom, start, end, name, place_hold, strand = bline.strip().split('\t')
    max_bin = genome[chrom] // rs
    if strand == '+':
        sbin = int(start) // rs
    else :
        sbin = int( end ) // rs
    if sbin > max_bin:
        sys.stderr.write( 'promoter big than max_bin for line %s, max bin is %s\n' % (bline,sbin) )
        return []
    lst = list(range( - kwargs.get('span_bins'), 0, 1))
    lst.extend(list(range(0, kwargs.get('span_bins') + 1, 1)))
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


def already( **kwargs ):
    infor = defaultdict( lambda : defaultdict( lambda : defaultdict( list )))
    for tab in system.dir( kwargs.get('output_root_dir')  ).fls('.*{}.*spanBin{}.*tab$'.format( kwargs.get('rs_inKb'), kwargs.get('span_bins'))):
        #if os.path.getsize( tab ) < 10 :
        #    continue
        tab_arr = os.path.basename( tab ).split('.')
        chain = [ i for i in tab_arr if i in ['+','-'] ][0]
        copy_index, chrom_bin, chrom, chrom_index, gene_bed_bin = 0, 0, 0, '', 0
        for i,j in enumerate ( tab_arr ):
            if 'copy' in j :
                copy_index = i
            if 'spanBin' in j:
                chrom_bin = int( j.replace('spanBin','') )
            if 'chr' in j :
                chrom,chrom_index = j,i
            if 'gbin' in j :
                gbin = int(j.replace('gbin',''))
        start,end = [ tab_arr[i] for i in ( chrom_index + 1, chrom_index +2) ]
        gene = '.'.join( os.path.basename( tab ).split('.')[ 0: copy_index + 1] )
        infor[ chrom_bin ][chrom][gbin].append( tab )
        print ( '!Already for {}'.format( tab ), file = sys.stderr)
    return infor

def link_same_pos( chrom, gbin, key, already_infor, **kwargs ):
    span_bins = kwargs.get('span_bins')
    if span_bins in already_infor and chrom in already_infor[ span_bins ] and gbin in already_infor[ span_bins ][ chrom ]:
        tab = kwargs.get('tab') and kwargs.get('tab') or os.path.basename(already_infor[span_bins][chrom][gbin][0])
        pdf = fix.fix( tab ).change('pdf')
        key = legal( key )
        key_tab = fix.fix(key).append( '{}.spanBin{}.gbin{}.tab'.format( kwargs.get('rs_inKb'), span_bins, gbin) )
        key_pdf = fix.fix(key_tab).change('pdf')
        if key_tab != tab:
            if not os.path.exists( os.path.join( kwargs.get('output_root_dir'), key_tab)):
                print ('Same pos for {} with {}, build link only...\n'.format( key_tab, tab), file = sys.stderr )
                #print ('ignore {}, delete: {}'.format(line, kwargs.get('dchr')), file = sys.stderr)
                os.system('cd {}\nln -s {} {}'.format( os.path.abspath(kwargs.get('output_root_dir')), tab, key_tab))
                os.system('cd {}\nln -s {} {}'.format( os.path.abspath(kwargs.get('output_root_dir')), pdf, key_pdf))
                already_infor[span_bins][chrom][gbin].append( key_tab )
        return already_infor, False
    if kwargs.get('tab') :
        print ( kwargs, already_infor.get( span_bins ).get(chrom).get(gbin), kwargs.get('tab'), key )
        exit('test')
    return already_infor, True

def bed2bin( **kwargs ):
    already_genes, tmp_bed = defaultdict( list ), open( '{}.bed'.format(''.join(chr(random.randrange(65,90)) for i in range(6))), 'w')
    already_infor, span_bins = already( **kwargs ), str(kwargs.get('span_bins') )
    need_deal_num_line = 0
    with open( kwargs.get('bed') ) as fh:
        for line in fh:
            line_arr = line.rstrip().split('\t')
            chrom, start, end = line_arr[0:3]
            if chrom in kwargs.get('dchr'):
                continue
            chain = line_arr[-1]
            if chain == '+':
                gbin = int( start ) // kwargs.get('rs')
            else :
                gbin = int( end ) // kwargs.get('rs')
            line_arr.pop(4)
            key = '.'.join(line_arr)
            already_infor, action = link_same_pos( chrom, gbin, key, already_infor, **kwargs )
            if action:
                need_deal_num_line += 1
                tmp_bed.write( line )
    tmp_bed.close()
    print ('deal {} lines, after ignore, file: {}'.format(need_deal_num_line, tmp_bed.name), file = sys.stderr)
    if not need_deal_num_line:
        exit('!allDone')
    return already_infor, tmp_bed.name
def calready( chrom, start, end, name, ainfor ):
    if name in ainfor:
        for each in ainfor[name]:
            achrom, astart, aend = each
            astart,aend,start,end = [ int(i) for i in (astart,aend,start,end) ]
            if chrom == achrom and astart < start < aend  and astart < end < aend:
                return True 
    return False 

def need_chroms_get( **kwargs ):
    need_chroms = set()
    with open( kwargs.get('bed') ) as fh :
        for line in fh:
            chrom, start, end, name, place_hold, strand = line.strip().split('\t')
            if 'chrX' in chrom :
                continue
            need_chroms.add(chrom)
    return need_chroms 


def get_normal_ice_fls( **kwargs ):
    dit = defaultdict( lambda : defaultdict ( str ) )
    peirods = kwargs.get('peirod') and kwargs.get('peirod') or ['E50','E90','E120','0M','4M','45Y','20Y']
    for peirod in peirods:
        for line in system.run('project_beds_bws_bams_andSoOn_pos.py {} ICE {}'.format( peirod, kwargs.get('rs_inKb')), shell = True ):
            if '#' in line or 'ini' in line or not line.strip():
                continue
            chrom, fl = [ i for i in re.split( r'\s+', line.rstrip()) ]
            dit[peirod][chrom] = fl
    return dit


#def dir2dit(dirs, load = [], **kwargs):
def dir2dit( **kwargs ):
    ice_fls = kwargs.get('ice_fls')
    dit = defaultdict( lambda : defaultdict( list ) )
    for peirod in ice_fls :
        for chrom in ice_fls[peirod]:
            chrom_fl = ice_fls[peirod][chrom]
            #trick.set2dict(dit, peirod, 'work_dir', os.path.dirname(chrom_fl))
            if chrom in kwargs.get('need_chroms') and peirod in kwargs.get('peirod'):
                print('loading signal: {} to matrix for extract data'.format(chrom_fl), file = sys.stderr)
                arr = np.array(pd.read_csv( chrom_fl, header = None, sep = '\t'))
                #trick.set3dict(dit, peirod, 'chrom', chrom, arr)
                dit[peirod][chrom] = arr
            #else :
                #trick.set3dict(dit, peirod, 'chrom', chrom, abs_path)
    return dit

def legal( string ):
    lst = ['(',')','/']
    for each in lst:
        string = string.replace(each,'__name__')
    return string


if __name__ == '__main__':
    #tfh = NamedTemporaryFile('w+t', delete = False, dir = '.')
    kwargs = vars( args )
    kwargs.update({'span_bins': kwargs.get('span')//kwargs.get('rs'), 'rs_inKb': str(kwargs.get('rs')//1000) +'K'})
    output_root_dir = os.path.join( '.'.join(args.peirod), kwargs.get('rs_inKb'), 'spanBin' + str(kwargs.get('span_bins')) )
    if not os.path.exists( output_root_dir ):
        output_dir = os.makedirs( output_root_dir )
    kwargs.update( { 'output_root_dir': output_root_dir } )

    #check already for ignore and copy for same bin results
    already_infor, bed = bed2bin( **kwargs )
    print ( '{} for analysis ...'.format( bed ), file = sys.stderr )
    kwargs.update( {'bed': bed })
    #already signal chroms 
    kwargs.update({'ice_fls': get_normal_ice_fls( **kwargs )})
    #need chroms
    kwargs.update( {'need_chroms': need_chroms_get( **kwargs  ) })
    #assert os.path.getsize(tfh.name)
    #make signal file for speed
    kwargs.update( { 'make_signal_dit': make_signal( **kwargs  ) })
    kwargs.update({ 'signal': hic_dit( **kwargs )})
    kwargs.update({'array': dir2dit( **kwargs )})
    kwargs.update({'genome':  chromosome.chr('rh8').chr() } )

    fh = open( kwargs.get('bed') )
    for i,line in enumerate(fh):
        chrom, start, end, name, place_hold, strand = line.strip().split('\t')
        #ignore already bin matched results
        if strand == '+':
        	gbin = int( start ) // kwargs.get('rs')
        else :
        	gbin = int( end ) // kwargs.get('rs')
        tab = '%s.tab' % ('.'.join([ chrom, start, end, name, strand,kwargs.get('rs_inKb'),'spanBin{}.gbin{}'.format(kwargs.get('span_bins'), gbin)]))
        key = '.'.join( [ chrom, start, end, name, strand ] )
        already_infor, action = link_same_pos( chrom, gbin, key, already_infor, **kwargs)
        if not action :
            continue
        else :
            already_infor[ kwargs.get('span_bins') ][ chrom ][ gbin ].append(tab)
        #if chrom in dchr:
        #    continue
        coord, pvalues, postive_bins = bline_pick_pvalue( line, **kwargs )
        if coord == []:
            continue
        #write the result in the bedGraph
        pdf = tab.replace( '.tab', '.pdf' )
        tab, pdf = [ os.path.join( kwargs.get('output_root_dir'), legal(i)) for i in (tab, pdf) ]
        print ( 'Ouput {} and {} to {} ... '.format( tab, pdf, kwargs.get('output_root_dir')), file = sys.stderr )
        if 1 :
            bedGraph( chrom, coord, pvalues, kwargs.get('rs'), tab )
            df = pd.DataFrame.from_dict(pvalues)
            peirods = list( df )
            #df = dataframe.df(df).insert(49, [ 0 for i in peirods ])
            ax = df.plot(style='.-')
            #ax.set_ylim(0, 10)
            ax.axvline(x= kwargs.get('span_bins'), linestyle = '-.')
            ax.axhline(y=3, linestyle = '-.')
            plt.savefig( pdf, format='pdf')
            plt.close()
    fh.close()


















