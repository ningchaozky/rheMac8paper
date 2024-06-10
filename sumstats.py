#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import gzip
import argparse
from ningchao.nBio import gwas
from tempfile import NamedTemporaryFile
from ningchao.nSys import trick, system, fix, status, fileKit

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help ='dir for sumstats')
parser.add_argument('-p', nargs='?', help ='pattern for the file find', default = 'sum')
parser.add_argument('-fp', nargs='?', help ='fiter pattern for the file find', default = 'out')
parser.add_argument('-sp', nargs='?', help ='spliter for the line', default = '\t')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()




def prase( args ):
    fls = system.dir(args.dir).fls( pattern = args.p, level = 0, fpattern = args.fp )
    return fls

def check( fl ):
    f = NamedTemporaryFile(prefix='mytemp', suffix='.txt', dir='/tmp', delete = False )
    fh, pcol = fileKit.File( fl ).open(), 0
    header = fh.next().rstrip().split( args.sp )
    f.write( '\t'.join(header) + '\n' )
    default_cnames = gwas.sumstats().header()
    for i,v in enumerate( header ):
        if v in default_cnames and default_cnames[v] == 'P':
            pcol = i
    if not pcol :
        exit('Check the head, there is no pvalue\n')
    for line in fh:
        line_arr = line.rstrip('\n').split( args.sp )
        pval = line_arr[pcol]
        pval_arr = pval.split('-')
        if len(pval_arr) == 2:
            num = int( pval_arr[1] )
            if num > 30 :
                pval_arr[1] = '30'
                sys.stderr.write( 'More than 10^20, eq to 10^20 for %s:%s\n' % (fl, line))
        line_arr[pcol] = str('-'.join(pval_arr))
        f.write( '\t'.join(line_arr) + '\n' )
    f.close()
    return f.name

def header_len( fl ):
    fh = fileKit.File( fl ).open()
    out = len(fh.next().split( args.sp ))
    fh.close()
    return out

if __name__ == '__main__':
    fls = prase( args )
    for fl in fls:
        if header_len( fl ) >= 6:
            wd = os.getcwd() 
            out = fl.split('.')[0] + '.out'
            if '.txt' in out:
                out = out.replace('.txt','')
            stats = status.cmd( wd ).stats()
            if out not in stats:
                fl = check( fl )
                cmd = 'munge_sumstats.py --sumstats %s --merge-alleles /home/ningch/soft/gwas/ldsc/w_hm3.snplist --out %s && gunzip %s.sumstats.gz' % ( fl, out, out)
                cmd = status.cmd( wd, out).accessory( cmd, appendix = 1 )[0]
                print(cmd)




























