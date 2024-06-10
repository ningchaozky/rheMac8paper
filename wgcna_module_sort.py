#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, system
from ningchao.nBio import rheMac

example = '''select.moduleEigengenes.tab'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs = '?', help = 'tab to sort' )
parser.add_argument( '-cut', nargs = '?', help = 'cut for the pvalue', default = 0.7, type = float )
parser.add_argument( '-pf', nargs = '?', help = 'print order txt file' )
parser.add_argument( '-fp', nargs = '+', help = 'fiter periods', default = ['E50','20Y'] )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    tab = args.tab
    days = rheMac.rheMac().days()
    #peirods = rheMac.rheMac().peirods( dp = ['E80','liver'] )
    days.pop('E80')
    peirods = system.dir.periods( prefrontal_raw = True )
    print ( peirods, file = sys.stderr )
    return tab, args.cut, days, peirods


if __name__ == '__main__':
    tab, cut, days, peirods = parse( args )
    tfh, infor = open( tab ), {}
    header = next(tfh).rstrip().split('\t')
    header[0] = 'color'
    header = [ i.split('.')[0] for i in header ]
    dit = {}
    for line in tfh:
        line_arr = line.rstrip().split('\t')
        color = line_arr.pop( 0 )
        values = [ float(i) for i in line_arr ]
        tmp_dit = dict(list(zip(header[1:], values)))
        #p,v = sorted(list(tmp_dit.items()), lambda x, y: cmp(x[1], y[1]), reverse = True )[0]
        p,v = sorted(list(tmp_dit.items()), key = lambda x: float(x[1]), reverse = True )[0]
        print ( color, p, v, line_arr )
        trick.dinit( dit, p, [])
        if v > cut :
            dit[p].append( color )
    orderFl = open( args.tab + '.order', 'w' )
    All_Fls = open( args.tab + '.orderFls.sh', 'w' )
    print ( 'enrich_celltype_hypergeometric_distribution_v3.py ', end = '', file = All_Fls )
    for p in peirods:
        names = dit[p]
        #names = [ i.replace('ME','') for i in dit[p]]
        print (p, file = sys.stderr )
        Fls = open( p + '.'+ args.tab + '.orderFls.sh', 'w' )
        print ( 'enrich_celltype_hypergeometric_distribution_v3.py ', end = '', file = Fls )
        for fl in names:
                color = fl.replace('ME','')
                fls = system.dir('.').fls( '\.'+ fl.replace('ME','')+'\.' )
                for fl in fls:
                    if 'Cytoscape' in fl or not fl.endswith('.txt'):
                        continue
                    print ( fl, ',index_col,2', sep = '', end = ' ', file = Fls )
                    if p not in args.fp:
                        print ( fl, ',index_col,2', sep = '', end = ' ', file = All_Fls )
        print ( '-pt heatmap', end = '', file = Fls )
        print( *names, sep = '\n', file = orderFl)

    print ( '-pt heatmap', end = '', file = All_Fls )




