#!/usr/bin/env python3
import os
import sys
import argparse
import itertools  
from ningchao.nSys import trick, system

example = '''promoter like value pick up keep pace value'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help ='wgcna work dir')
parser.add_argument('newOrder', nargs='?', help ='new order for deal')
parser.add_argument('-RNA', nargs='?', help ='tpm file', default = '/dataB/ftp/pub/rheMac3/prefrontal/RNA/tpm/gene.tpm.sort.mini' )

if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    fls = list( system.dir( args.dir ).fls( 'txt', fpattern = 'Cytoscape'))
    moduleEigengenes = list( system.dir( args.dir ).fls( 'moduleEigengenes.pvalue.tab'))[0]
    eig_order = Eigengenes_get_order( moduleEigengenes )
    colors = [ i.split('.')[-2] for i in fls ]
    fls = dict(zip( colors, [ fl_pick_genes(fl) for fl in fls] )) 
    rna_order, rna_raw = order_RNA( args.RNA )
    return fls, rna_raw, rna_order, eig_order


def order_RNA( RNA ):
    fh, out, raw = open( RNA ), {}, {}
    header = next( fh ).strip()
    trick.dinit( raw, 'header', header )
    header = [i.split('.')[0] for i in  header.split('\t') ]
    for line in fh:
        line_arr = line.strip().split('\t')
        dit = dict(zip(header[1:], [ float(i) for i in line_arr[1:]]))
        order = sorted(dit.items(), key = lambda x: x[1], reverse = True )
        order = [i[0] for i in order ]
        trick.dinit( raw, line_arr[0], line.strip() )
        trick.dinit( out, line_arr[0], order )
    return out, raw

def fl_pick_genes( fl ):
    fh, lst = open( fl ), []
    for line in fh:
        line = line.strip()
        if 'modProbes' in line or not line:
            continue
        lst.append( line.split('\t')[1] )
    return lst

def Eigengenes_get_order( fl ):
    fh, out = open( fl ), {}
    header = next( fh ).strip().split('\t')
    for line in fh:
        line_arr = line.strip().split('\t')
        color = line_arr.pop(0)
        dit = dict(zip( header, [ float(i) for i in line_arr[1:]]) )
        order = sorted(dit.items(), key = lambda x:x[1], reverse = True )
        order = [ i[0] for i in order ]
        trick.dinit( out, color, order )
    return out

def ignore( gval, cutoff = 1 ):
    g_arr = [ float(i) for i in gval.split('\t')[1:] ]
    g_arr = sorted( g_arr )
    return g_arr[0] < cutoff

def main():
    fls, rna_raw, rna_order, eig_order = parse( args )
    num, peirod = range( 1000 ), ''
    print( rna_raw['header'] )
    for color in open( args.newOrder ):
        color = color.replace('ME','').strip()
        if color in fls:
            sys.stderr.write('Deal %s %s\n' % (color, peirod))
            i = 0
            fl_genes = fls[color]
            order = eig_order['ME' + color]
            for gene in fl_genes:
                if ignore( rna_raw[gene] ) :
                    continue
                if peirod == 'other':
                    print ( rna_raw[gene] )
                    continue
                pe = 2
                combins = list(itertools.combinations(order[0:pe],len(order[0:pe])))
                i += 1
                if peirod in rna_order[gene][0:pe]:
                    print( rna_raw[gene] )
                #else :
                if i not in num:
                    sys.stderr.write('%s num is %d\n' % ( color, i ))
                    break
                    
        else :
            peirod = color

if __name__ == '__main__':
    main()
    

























