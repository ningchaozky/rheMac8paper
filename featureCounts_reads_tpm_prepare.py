#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
from ningchao.nSys import trick
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'mats', nargs = '+', help = 'featureCounts output file')
parser.add_argument( '-gene_mean_length','-gl', nargs = '?', help = 'gene mean lenght file output from gtf', default = '/home/soft/data/genome/mm10/mm10.ncbiRefSeq.gtf.mean' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def merge( mats ):
    dfs = []
    for mat in mats:
        df = pd.read_csv( mat, index_col = 0, comment='#', sep = '\t')
        dfs.append( df )
    return pd.concat( dfs, axis = 1 )



if __name__ == '__main__':
    df = merge( args.mats )
    len_df = pd.read_csv( args.gene_mean_length, index_col = 0, sep = '\t', )
    df = pd.concat( [df, len_df], axis = 1)
    #columns = [ i for i in df.columns if 'bam' in i ]
    columns = [ i for i in df.columns ]
    df = df.loc[ :, [ *columns, 'mean'] ]
    df.to_csv( sys.stdout, sep = '\t', index_label = 'gene' )

































