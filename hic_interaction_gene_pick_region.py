#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,system,excel
from ningchao.nBio import bed,rheMac

example = '''txt one gene for one line. Can give the column for pick gene. The gene should upper'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='+', help = 'reference' )
parser.add_argument( '-dir', nargs='?', help = 'dir for interaction results', default = '/dataB/ftp/pub/rheMac3/prefrontal/hic/whole_genome/coding/index' )
parser.add_argument( '-pcut', nargs='?', help = 'p or q value cut off', default = 5, type = float )
parser.add_argument( '-output_inter_region', nargs = '?', help = 'output for inter_region', required = True )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def parse( args ):
    tab, oiregion = args.tab, args.output_inter_region
    if len( tab ) == 1 :
        tcol = 0
        tab = tab[0]
    elif len(tab) == 2 :
        tab, tcol = tab[0], int(tab[1]) - 1
    return tab, tcol, args.dir, oiregion

def tabs( directory ):
    dit, tabs = {}, system.dir(directory).fls(pattern = 'tab', fpattern ='pdf')
    for tab in tabs:
        key = os.path.basename(tab).strip().split('.')[0]
        trick.dit(dit).set(key,[])
        dit[key].append(tab)
    return dit

def tab_get_region( tabs, pcut = 5 ):
    out = {}
    for tab in tabs:
        fh = open( tab )
        header = fh.next().strip().split('\t')
        for line in fh:
            line_arr = line.strip().split('\t')
            dit = dict(list(zip(header,line_arr)))
            start = dit.pop('start')
            end = dit.pop('end')
            chrom = dit.pop('chr')
            for key in dit:
                val = float(dit[key])
                if val > pcut:
                    gkey = os.path.basename(tab).split('.')
                    gkey.pop(0)
                    gkey.pop(-1)
                    gkey = '.'.join(gkey)
                    trick.dit(out).set(key,gkey, [])
                    out[key][gkey].append( '\t'.join( [chrom, start, end] ) )
        fh.close()
    for key in out:
        for gkey in list(out[key].keys()):
            out[key][gkey] = ','.join(bed.bed( out[key][gkey] ).merge())
    return out

def write_tmp( oiregion, gene) :
    gchrom = rheMac.rheMac().gchrom()
    oirfh,fh = open(oiregion, 'w'), open( gene )
    fh = open( gene )
    header = excel.xls('/dataB/ftp/pub/rheMac3/prefrontal/hic/whole_genome/coding/index/PAX6.chr14.34086109.34088109.+.tab').header()
    header.pop( 0 ),header.pop( 0 ),header.pop( 0 ),header.insert(0, 'gene')
    keys = header[1:]
    oirfh.write( '\t'.join( header ) + '\n' )
    for line in fh:
        line_arr = line.strip().split('\t')
        gene = trick.lst( line_arr ).get(gcol)
        gene = gene.strip().split('.')[0].replace('/','--')
        out_line = [ gene ]
        if gene in infor:
            merge_dit = tab_get_region( infor[gene] )
            lst = trick.dit(merge_dit).get(keys)
            lst = [ str(i) for i in lst ]
            out_line.extend( lst )
        if len( out_line ) != 1 :
            oirfh.write( '\t'.join( out_line ) + '\n' )
        else :
            if gene in gchrom:
                chrom = str(gchrom[gene])
            else :
                chrom = gene
            sys.stderr.write('Ignore gene %s, chrom %s: %s\n' % (gene, chrom, line) )


if __name__ == '__main__':
    gene, gcol, directory, oiregion = parse( args )
    infor = tabs(directory)
    write_tmp( oiregion, gene )
