#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, fix, system
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'exp', nargs = '?', help = 'reference')
parser.add_argument( '-fc', nargs = '?', help = 'fold change', default = 2, type = float)
parser.add_argument( '-cut', nargs = '?', help = 'cut for row max(lina_arr)', default = 1, type = float)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def get_dynamic_and_stable():
    data = []
    with open( args.exp ) as f:
        header = next(f).strip()
        stable, dynamic = open( fix.fix(args.exp).append('stable'), 'w'), open( fix.fix(args.exp).append('dynamic'), 'w' )
        print ( header, file = stable )
        print ( header, file = dynamic )
        kwargs.update( {'header': header.strip().split('\t') })
        for line in f:
            line_arr = line.strip().split('\t')
            numbers = [ float(i) for i in line_arr[1:] ]
            if max( numbers ) < args.cut:
                continue
            if [ i for i in numbers if not i ]:
                print( line.strip(), file = dynamic )
                continue
            compare = []
            for i,v in enumerate( numbers ):
                compare.append( [ v/j for j in numbers if v/j > args.fc] )
            compare = [ j for i in compare for j in i ]
            if compare :
                print ( line.strip(), file = dynamic )
                data.append( line.strip() )
            else :
                print ( line.strip(), file = stable )
    kwargs.update({'dynamic': dynamic.name})
    kwargs.update({'dat': data})

def data_sorted():
    output = open( fix.fix(args.exp).append('hardSorted'), 'w')
    print ( *kwargs.get('header'), sep = '\t', file = output )
    index_label = kwargs.get('header').pop( 0 )
    order = system.dir.sort( kwargs.get('header'), up = '', down = '',sort = 'peirod' )
    print ( order )
    for p in order :
        cmd = 'hard_pick.py {} -p {} -nfc {} -maxCut {}'.format( args.exp, p, args.fc, args.cut )
        out = system.run( cmd, shell = True)
        print ( out.pop(0), file = sys.stderr )
        for line in out:
            print ( line.strip(), file = output)
    output.close()

if __name__ == '__main__':
    kwargs = vars( args )
    get_dynamic_and_stable()
    data_sorted()

























