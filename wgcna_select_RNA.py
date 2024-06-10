#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, system
from collections import defaultdict
desc = '''indexmatch: moduleText rowdata RNAexpression; no index match the two table'''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'Input', nargs = '+', help = 'module wgcnasignaldata RNA')
parser.add_argument( '-m', nargs = '?', help = 'sort math length', default = 1, type = int)
parser.add_argument( '-heatmap', nargs = '?', help = 'heatmap plot')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



def short_symbol( symbol ):
    return symbol.split('.')[0].upper()


def need_symbols():
    infor, symbols = defaultdict( list ), set()
    if len( args.Input ) == 3 :
        [kwargs.update({ v: args.Input[i]}) for i,v in enumerate(['module', 'signal', 'rna'])]
        files = [ args.Input[0] ]
    elif len( args.Input ) == 2 :
        [kwargs.update({ v: args.Input[i] }) for i,v in enumerate(['signal', 'rna'])]
        files = args.Input[1:]
    for fl in files:
        with open( fl ) as f :
            firstLine = next( f )
            if 'modProbes' not in firstLine :
                f.seek(0)
            for line in f :
                line_arr = line.rstrip().split('\t')
                if len( args.Input ) == 2 :
                    symbol = short_symbol( line_arr[0] )
                elif len( args.Input ) == 3 :
                    symbol = short_symbol( line_arr[1] )
                symbols.add( symbol )
    kwargs.update({'symbols': symbols})
    return symbols

def mean( dit ):
    tdit = defaultdict( list )
    for key in dit :
        tdit[ key.split('.')[0] ].append( float(dit[key]) )
    for key in tdit.keys():
        tdit[key] = round( sum(tdit[key]) / len(tdit[key]), 3 )
    tdit['sorted'] = [ i[0] for i in sorted ( tdit.items(), key = lambda x: x[1], reverse = True ) ]
    return tdit

def indexGet():
    infor = defaultdict( lambda : defaultdict( dict ) )
    for tp in ('signal','rna'):
        print ( '###Deal with {} tab ...'.format( kwargs.get(tp) ), file = sys.stderr )
        with open( kwargs.get( tp ) ) as f :
            header = next( f ).rstrip().split('\t')[ 1: ]
            kwargs.update( { 'header': [ short_symbol(i) for i in header if short_symbol(i) not in ['E80']] } )
            for line in f :
                line_arr = line.rstrip().split('\t')
                symbol = short_symbol( line_arr.pop( 0 ) )
                line_dit = mean(dict(zip(header, line_arr)))
                if symbol in kwargs.get('symbols'):
                    infor[tp][symbol] = line_dit
    kwargs.update( infor )


def matchOut():
    signal, rna = kwargs.get('signal'), kwargs.get('rna')
    if len( args.Input ) == 3 :
        ifiles = args.Input[1:]
    elif len( args.Input ) == 2 :
        ifiles = args.Input
    signal_out, rna_out = [ open('.'.join([ os.path.basename(i), 'matched']), 'w') for i in ifiles ]
    print ( 'symbol', *kwargs.get('header'), sep = '\t', file = signal_out)
    print ( 'symbol', *kwargs.get('header'), sep = '\t', file = rna_out )
    for key in signal:
        signalPeirods = set(signal[key].pop('sorted')[0:args.m])
        if key in rna :
            rnaPeirods = set(rna[key].pop('sorted')[0:args.m])
            if signalPeirods == rnaPeirods:
                if not max( signal[key].values() ) or not max( rna[key].values() ):
                    print ( key, 'signal:', *signal[key].values(), 'rna:', *rna[key].values(), file = sys.stderr)
                    continue
                print ( key, *[ signal[key][i] for i in kwargs.get('header') if i in signal[key] ], sep = '\t', file = signal_out)
                print ( key, *[ rna[key][i] for i in kwargs.get('header') if i in rna[key] ], sep = '\t', file = rna_out)
    signal_out.close()
    rna_out.close()
    kwargs.update( {'signal_out': signal_out.name })
    kwargs.update( {'rna_out': rna_out.name })




if __name__ == '__main__':
    kwargs = vars( args )
    need_symbols()
    indexGet()
    matchOut()
    if args.heatmap :
        cmd = 'heatmap_seaborn.py {} -z_score 0'.format( kwargs.get('signal_out') )
        system.run( cmd, shell = True, noReturn = True )
        cmd = 'heatmap_seaborn.py {} -z_score 0 -cmap YlOrBr'.format( kwargs.get('rna_out') )
        system.run( cmd, shell = True, noReturn = True )


















