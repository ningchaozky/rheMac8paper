#!/usr/bin/env python3
import os
import sys
import argparse
import itertools
from ningchao.nSys import trick
example = ''' '''
gcs = ['hg19','hg38','rheMac8','rheMac3','rheMac10','mm10','mm9','dm6']
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('g', choices=gcs, help = 'genome you want to download')
parser.add_argument('-dir', nargs='?', help = 'dir want to put', default = '/var/lib/mysql/')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def parse( args ):
    g = args.g
    odir = args.dir and args.dir or '/home/mysql_data/mysql/mysql/'
    flst = ["all_est","chromInfo","gap","hgFindSpec","refFlat","refGene","refSeqAli","trackDb","xenoRefGene","xenoRefSeqAli","knownCanonical",'phyloP124way','ncbiRefSeq','ncbiRefSeqCurated','multiz124way','altSeqLiftOverPsl','gtexInfo']
    flst.extend(['rmsk'])
    typ = ['MYD','MYI','frm']
    dd = g not in odir and '/'.join( [ odir, g ] ).replace('//','/') or odir
    #print ('rsync -avzP rsync://hgdownload.cse.ucsc.edu/mysql/hgFixed/ %s' % os.path.dirname(dd) )
    print ("sudo rsync -avzP rsync://hgdownload.cse.ucsc.edu/mysql/hgFixed/gtexInfo.* /var/lib/mysql/hgFixed")
    for i in itertools.product( flst, typ ):
        link = 'sudo rsync -avzP rsync://hgdownload.cse.ucsc.edu/mysql/%s/%s %s' % ( g, '.'.join(i), dd )
        print ( link )
    print ( 'sudo chown -R mysql.mysql %s' % dd )
    dd = '/gbdb/%s' % g
    print ('sudo rsync -avzP rsync://hgdownload.cse.ucsc.edu/gbdb/%s/%s.2bit %s' % (g, g, dd ) )
    wibs = []
    if g == 'dm6':
        wibs = ['phyloP124way','phastCons124way']
    for wib in wibs:
        link = 'hgdownload.cse.ucsc.edu/gbdb/{}/multiz124way/{}.wib'.format(g, wib)
        fl = link.replace('hgdownload.cse.ucsc.edu', '')
        print ( 'sudo mkdir -p {}'.format( os.path.dirname(fl) ))
        print ('sudo rsync -avzP rsync://{} {}'.format( link, fl) )
if __name__ == '__main__':
    parse( args )






