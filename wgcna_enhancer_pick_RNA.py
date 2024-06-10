#!/usr/bin/env python3
import os
import sys
import importlib
import argparse
from ningchao.nSys import trick, system
from ningchao.nBio import rheMac

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'order', nargs='?', help = 'order for pub list')
parser.add_argument( '-match', nargs='?', help = 'match to gene', default = '/dataB/ftp/pub/rheMac3/prefrontal/enh/enhaner/match/enhancer.merge.bed.gene.addNum')
parser.add_argument( '-p', nargs='*', help = 'peirod you want to select', default = [ 'all' ])
parser.add_argument( '-tpm', nargs='?', help = 'RNA tpm list', default = '/dataB/ftp/pub/rheMac3/prefrontal/RNA/tpm/gene.tpm.sort.mini' )
parser.add_argument( '-w', nargs='?', help = 'weight for each module', default = 50, type = int )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def get_color_lst( order, peirod, weight):
    peirods = rheMac.rheMac().peirods() 
    peirods.append('other')
    fh, dit, key, ratio = open( order ), {}, '', {}
    for line in fh:
        line = line.strip('\n').replace('ME','')
        if line in peirods:
            key = line
            trick.dinit( dit, key, [] )
        else :
            dit[key].append( line )
    for key in dit:
        trick.dinit( ratio, key, len( dit[key] ) * weight )
    return trick.dit( dit ).get( peirod, dit = True ), ratio

def color_txt():
    fls = list(system.dir( '.' ).fls('txt', fpattern = 'Cy'))
    names = [ i.split('.')[1] for i in fls ]
    return dict(list(zip(names, fls)))

def fl_parse( fls ):
    genes = []
    for fl in fls:
        if not fls[fl].strip():
            continue
        fh = open( fls[fl] )
        next(fh)
        for line in fh:
            line_arr = line.strip('\n').split('\t')
            genes.append( line_arr[1] )
    return genes

def get_genes( colors, txts ):
    enhs = []
    for peirod in colors:
        pcolors = colors[ peirod ]
        fls = trick.dit( txts ).get( pcolors, dit = True )
        yield peirod, fl_parse( fls )
 
def enh_pick_genes( enhs, match ):
    mfh, dit = open( match ), {}
    for line in mfh:
        line_arr = line.strip('\n').split('\t')
        trick.dinit( dit, line_arr[0], 'distal', line_arr[3] )
        trick.dinit( dit, line_arr[0], 'local', line_arr[4] )
    for peirod, penhs in enhs:
        for penh in penhs:
            distal = dit[penh]['distal']
            local = dit[penh]['local']
            genes = distal.split(',')
            genes.extend( local.split(',') )
            for gene in genes:
                if not gene.strip():
                    continue
                gene_arr = gene.split('.')
                try :
                    gene = '.'.join( gene_arr[ 0: trick.lst( gene_arr ).index('copy', regular = True ) + 1 ] )
                except :
                    print(trick.lst( gene_arr ).index('copy', regular = True ), gene_arr)
                    exit()
                yield peirod, gene

def RNA_parse( exp, fc = 1.2 ):
    fh, out = open( exp ), {}
    header = next( fh ).strip('\n').split('\t')
    for line in fh:
        line_arr = line.strip('\n').split('\t')
        dit = dict( list(zip( header, line_arr)) )
        gene = dit.pop('gene')
        #trick.dinit( out, gene, [] )
        dit_sort = sorted(dit.items(), key = lambda x: float(x[1]), reverse = True )
        #dit_sort = sorted(list(dit.items()), lambda x, y: cmp( float(x[1]), float(y[1])), reverse = True )
        vmax, vmax2 = float( dit_sort[0][1] ), float( dit_sort[1][1] )
        if vmax == 0 :
            continue
        if vmax2 == 0:
            vmax2 = vmax2 + 0.1
        if vmax / vmax2 > fc:
            max_p = dit_sort[0][0].split('.')[0] 
            trick.dinit( out, gene, max_p )
            #out[ gene ].append([ max_p, dit_sort])
        else :
            trick.dinit( out, gene, 'other' )
            #out[ gene ].append( ['other', dit_sort ])
    return out
        

def main( args ):
    colors, ratio = get_color_lst( args.order, args.p, args.w)
    txts, tpm = color_txt(), args.tpm
    enhs = get_genes( colors, txts )
    rna_infor = RNA_parse( tpm )
    num, ugenes = 0, []
    for peirod, gene in enh_pick_genes( enhs, args.match ):
        if gene in rna_infor and peirod == rna_infor[gene]:
            if gene not in ugenes:
                num += 1
                if num < ratio[peirod]:
                    print('\t'.join( [ gene, str( ratio[peirod] ) ] ))
                ugenes.append( gene )
if __name__ == '__main__':
    main( args )




















