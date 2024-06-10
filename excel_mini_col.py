#!/usr/bin/env python3
import os
import re
import sys
import argparse
import fileinput
from ningchao.nSys import trick, system
from ningchao.nBio import rheMac

parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='mini the excel by column', formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'tab', nargs='?', help ='excel you want to mini', default = sys.stdin, type = argparse.FileType('r'))
parser.add_argument( '-c', choices=['fcol','scol','skey','fkey', 'fitRep', 'join','copyAsReps', 'sortPeriod','uniq_tab_marker'], help = 'which method you want to mini you excel', default = 'skey')
parser.add_argument( '-v', nargs='*', help = 'values for the choices')
parser.add_argument( '-js', nargs='*', help = 'if args.c == join, join seperator', default = '_')
parser.add_argument( '-o', nargs='?', help = 'output file', default = sys.stdout, type = argparse.FileType('w') )
parser.add_argument( '-norep', action='store_true', help = 'not include the rep')
parser.add_argument( '-nr', action='store_false', help = 'no use regular select skey')
parser.add_argument( '-mergeColumns','-m', nargs='+', help = '. 1 2 3')
parser.add_argument( '-sp', nargs = '?', help = 'separator \t|,', default = '\t')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()


def column( values, **kwargs ):
    tfh = kwargs.get('tab')
    for line in kwargs.get('lines'):
        line_arr = line.strip('\n').split(args.sp)
        print ( *[ kwargs.get('header')[ i ] for i in values ], file = kwargs.get('o'), sep = args.sp )
    for line in tfh:
        line_arr = line.strip('\n').split(args.sp)
        if kwargs.get('c') == 'scol':
            print ( *[ line_arr[ i ] for i in values ], file = kwargs.get('o'), sep = args.sp )
        elif kwargs.get('c') == 'fcol':
            print ( * [ v for i,v in enumerate( line_arr ) if i not in values ], sep = args.sp, file = kwargs.get('o') )
        elif kwargs.get('c') == 'sortPeriod' :
            print ('yes')
        else :
            print ( kwargs )
            print('Please select the right [ scol, fcol ]\n')
            exit()

def get_cols( values, header, rep, copy):
    ints, strs, cols_int, cols_str = [],[],[],[]
    for each in values:
        try :
            ints.append( int(each) - 1 )
        except :
            strs.append( each )
    cols_int.extend( [ i for i in ints ] )
    cols_str.extend( trick.lst( header ).get( strs, regular = True ) )
    index = trick.lst(header).index( cols_str )
    keys = rheMac.trick().sort( trick.lst(header).get( index ), copy = copy )
    if rep:
        keys = [ i for i in keys if 'rep' not in i ]
    cols = cols_int
    cols.extend( trick.lst( header ).index( keys ) )
    return cols


def key( values, **kwargs):
    tfh = args.tab
    for j,line in enumerate(tfh):
        if not line.strip():
            continue
        line_arr = line.rstrip('\n').split( args.sp )
        if kwargs.get('mergeColumns') :
            osep = kwargs.get('mergeColumns')[0]
            pos = [ int(i) - 1 for i in kwargs.get('mergeColumns')[1:] ]
            values = [ i for i in values if i not in pos ]
            oline_arr = [ osep.join([ line_arr[i] for i in pos ]) ]
            oline_arr.extend( [ line_arr[ i ] for i in values ] )
        else :
            try :
                oline_arr = [ line_arr[i] for i in values ]
            except :
                print ( tfh.name, j, file = sys.stderr)
                print ( values, line_arr, line, file = sys.stderr)
                continue
        print ( *oline_arr, file = kwargs.get('o'), sep = args.sp )

def fitRep( values, **kwargs ):
    fh = kwargs.get('tab')
    fh.seek(0)
    for line in fh:
        line_arr = line.strip().split( args.sp )
        print ( *[ line_arr[i] for i in values ], sep = args.sp, file = kwargs.get('o'))

def expand_colon( string ):
    if string.startswith(':'):
        return '0-{}'.format(string.replace(':',''))
    elif string.endswith(':'):
        return '{}-{}'.format(string.replace(':',''), str(len(kwargs.get('header')) + 1))
    else :
        return string

def pos_deal_with():
    #check values is num or string
    header = kwargs.get('header')
    if args.c == 'fitRep':
        #print ( system.dir.str_map_rep(header[0]) )
        pos = [ i for i,v in enumerate( header ) if not re.search( r'rep', v) ]
        return pos
    if 'sortPeriod' in args.c:
        kwargs.update({'v': [ header[0], *system.dir.sort( header[1:] ) ] })
    varr = [ expand_colon(i) for i in kwargs.get('v') ]
    values_str = ' '.join( varr )
    values_arr = re.split(r'\s+|,', values_str)
    try :
        num = True
        int(values_arr[0])
    except :
        num = False
    if '-' in values_arr[0] :
        num = True
    #get raw 
    raw, pos = [], []
    header = kwargs.get('header')
    if num :
        for each in values_arr:
            if '-' in each:
                seRange = [ int(i) - 1 for i in each.split('-') ]
                if len( seRange ) == 2 :
                    seRange[-1] += 1
                    raw.extend( range( *seRange ) )
                else :
                    exit('- region only two numbers. Eg: 1-3')
            elif ',' in each:
                raw.extend( [ int(i) - 1 for i in each.split(',') ] )
            else :
                raw.append( int( each ) - 1 )
    else :
        raw_pos = trick.lst( header ).index( kwargs.get('v'), regular = args.nr, returnType = 'raw' )
        raw.extend( raw_pos )
    if 'f' in args.c:
        pos.extend( [ i for i,v in enumerate( header ) if i in raw ])
    else :
        pos.extend( raw )
    print ( 'Select pos is: ', pos, file = sys.stderr )
    kwargs.update({'pos': pos })
    return pos


def check():
    mfh = args.tab
    header = mfh.readline().strip('\n')
    while '#' in header :
        header = mfh.readline()
    header_arr = header.split( args.sp )
    kwargs.update({'lines': [ header ] }); kwargs.update({'header': header_arr })
    first_line = mfh.readline().strip('\n'); kwargs.get('lines').append( first_line )
    first_line_arr = first_line.split( args.sp )
    if kwargs.get('v') or args.c in [ 'fitRep','copyAsReps', 'sortPeriod']:
        print ('Pass check parameters...', file = sys.stderr )
        print ('header is {}...'.format( header ), file = sys.stderr)
        return True
    else :
        line_len, header_len = len( first_line_arr ), len( header_arr )
        if line_len != header_len :
            print( '!!!Warning, header lenght: {}, line lenght: {}'.format(header_len,line_len), 'not the same')
        if 1 :
            lst = sorted( [ '|'.join([ '|'.join(v), str(i+1)]) for i,v in enumerate( zip(header_arr, first_line_arr) )] )
            if len( lst ) <= 1 :
                print ( 'len header is <= 1, check the args.sp: {}'.format( args.sp ))
            print( *lst, sep = '\n' )
            exit()
def join_col( values, **kwargs):
    tfh = args.tab
    tfh.seek(0)
    for line in tfh:
        line_arr = line.strip().split( args.sp )
        key = kwargs.get('js').join([ line_arr[i] for i in values])
        line_arr = [ v for i,v in enumerate( line_arr ) if i not in values ]
        line_arr.insert( 0, key )
        print ( *line_arr, sep = args.sp, file = sys.stdout )


def copy_as_replicates( include_merge = True ):
    args.tab.seek(0)
    for line in args.tab:
        line_arr = line.rstrip().split( args.sp )
        line_output = []
        for i,v in enumerate(line_arr):
            if i in kwargs.get('pos'):
                line_output.extend([v,v])
            line_output.append(v)
        print ( *line_output, sep = args.sp, file = args.o )

def uniq_tab_marker( f ):
    if args.c == 'uniq_tab_marker':
        header = f.readline().rstrip().split('\t')
        print ( *[ '.'.join(['c',v, str(i)]) for i,v in enumerate( header ) ], sep = '\t' )
        for i,line in enumerate(f):
            line_arr = line.rstrip().split('\t')
            line_arr[0] = '.'.join([ 'r', line_arr[0], str(i) ])
            print ( *line_arr, sep = '\t' )
        exit()
def run( values ):
    if args.c == 'fitRep':
        fitRep( values, **kwargs )
    elif args.c in ['scol','fcol'] :
        column( values, **kwargs )
    elif args.c in ['skey','fkey'] :
        key( values )
    elif args.c == 'join':
        join_col( values, **kwargs )
    elif args.c == 'copyAsReps':
        copy_as_replicates()
    elif args.c == 'sortPeriod':
        key( values )


def main():
    uniq_tab_marker( args.tab )
    check()
    values = pos_deal_with()
    for line in kwargs.get('lines'):
        line_arr = line.rstrip('\n').split( args.sp )
        print ( *[line_arr[i] for i in values ], sep = args.sp, file = args.o )
    run( values )
if __name__ == '__main__':
    #tfh, ofh, choice, values, rep, copy = parse( args )
    kwargs = vars(args)
    main()







