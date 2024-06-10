#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick, fix
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('tab', nargs='?', help = 'NMF the methylation data')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



def data_format( arr, fiter = 0.3 ):
    lst, nan_num = [], 0
    for each in arr:
        if not each.strip():
            print ('There is null string', file = sys.stderr)
            each = 2
        if each in ['NA', 'NaN'] :
            nan_num += 1
            each = 2
        elif float( each ) > 3 :
            each  = 2
        lst.append( each )
    if len(arr) - nan_num == 3 :
        return False
    else :
        return lst, nan_num
    #if 1 - nan_num/len(arr) > fiter:
    #    return lst, nan_num
    #else :
    #    print ( 'fiter line:', len(arr), nan_num, 1 - nan_num/len(arr), fiter )
    #    return False
def deal_matrix( mat ):
    with open( mat ) as f:
        ofh = open( fix.fix(mat).append('colMergeWithNum'), 'w')
        header = f.readline().rstrip().split('\t')
        if 'chr' in header[0] and '0000' in header[2] :
            header = [ 'cell_{}'.format(i+1) for i in range(len(header)-3) ]
            #print ( 'pos', *header, sep = '\t', file = ofh )
            f.seek(0)
        print ( 0, *header, sep = '\t', file = ofh )
        for i, line in enumerate(f):
            line_arr = line.rstrip().split('\t')
            #NA_num = sum([ 1 for i in line_arr if 'NA' in i ])
            #if NA_num + 3 == len( line_arr ):
            #    print ( '!!!fit line: {}'.format( '\t'.join(line_arr[:3]), file = sys.stderr ) )
            #    continue
            #line_arr[3:] = map( data_format, line_arr[3:] )
            line_format = data_format( line_arr[3:], fiter = 0.5 )
            if line_format:
                line_arr[3:], nan_num = line_format
                key = '.'.join( line_arr[:3] )
                print ( i, *line_arr[3:], sep = '\t', file = ofh )
        ofh.close()
        return ofh.name


if __name__ == '__main__':
    mat = deal_matrix( args.tab )
























