#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help = 'dir for del')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def hic():
    dels = []
    dir_name = args.dir
    fls = os.listdir(dir_name)
    fiters = [ i for i in fls if 'hicup_filter_ditag_rejects' in i ]
    sams = [ i for i in fls if i.endswith('sam') ]
    bams = [ i for i in fls if i.endswith('bam') and i != '%s.hicup.bam' % dir_name ]
    fqs = [ i for i in fls if i.endswith('fq') ]
    svgs = [ i for i in fls if i.endswith('svg') ]
    dels.extend( fiters )
    dels.extend( sams )
    dels.extend( bams )
    dels.extend( fqs )
    return [ os.path.join(dir_name, i) for i in dels ]

print('\n'.join(hic()))

























