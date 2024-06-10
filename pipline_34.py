#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick,system
from ningchao.nBio import rheMac

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'ini', nargs = '?', help = 'ini file')
parser.add_argument( '-s', nargs = '?', help = 'specify', default = 'PFC' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



bw_dir = os.path.abspath( os.path.join('.', 'bw') )
infor = rheMac.infor( args.ini, args.s )
infor_period = ( infor.raw()['period'] )
for typ in infor_period:
    for region in infor_period[ typ ]:
        for deal in infor_period[ typ ][region]:
            for rep in infor_period[ typ ][ region ][ deal ]:
                name = infor_period[ typ ][ region ][ deal][ rep ] 
                wd = os.path.join('.', name )
                for fl in system.dir( wd ).fls('{}.*bw$'.format( name )):
                    cmd = 'ln -s {} {}/{}'.format( fl, bw_dir, '.'.join([deal,rep,'bw']))
                    print ( cmd )




























