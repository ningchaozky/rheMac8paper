#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick
from ningchao.nBio import rheMac

example = '''select.moduleEigengenes.tab'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs = '?', help = 'tab to sort' )
parser.add_argument( '-cut', nargs = '?', help = 'cut for the pvalue', default = 0, type = float )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse( args ):
    tab = args.tab
    days = rheMac.rheMac().days()
    days.pop('E80')
    return tab, args.cut, days

def delta( pd, days ):
    dit = {}
    for p in days:
        pday = days[p]
        if days[pd] - pday :
            trick.dinit( dit, p, 1.0/abs((days[pd] - pday)) )
        else :
            trick.dinit( dit, p, 0.000035 )
    return dit



if __name__ == '__main__':
    tab, cut, days = parse( args )
    tfh, infor = open( tab ), {}
    header = tfh.next().rstrip().split('\t')
    header[0] = 'color'
    header = [ i.split('.')[0] for i in header ]
    for line in tfh:
        line_arr = line.rstrip().split('\t')
        color = line_arr.pop( 0 )
        values = [ float(i) for i in line_arr ]
        tmp_dit = dict(list(zip(header[1:], values)))
        p,v = sorted(list(tmp_dit.items()), lambda x, y: cmp(x[1], y[1]), reverse = True )[0]
        if v < 0.75 :
            continue
        #elif round(v, 4 ) == 0.8256 :
        #    continue
        delta_dit =  delta( p, days )
        sys.stderr.write(p + str(sorted(list(delta_dit.items()), lambda x, y: cmp(x[1], y[1]), reverse=False)))
        for p in delta_dit:
            dlta = delta_dit[p]
            #if dlta > 0.05 :
            #    continue
            trick.dinit( infor, dlta, [])
            if tmp_dit[p] > 0 :
                infor[ dlta ].append( (tmp_dit[p],p) )
    trick.write( 'index','Diff','cor','peirod' )
    i = 0
    for dlta in infor:
        for val in infor[dlta]:
            i += 1
            trick.write( i,dlta, val[0],val[1] )



























