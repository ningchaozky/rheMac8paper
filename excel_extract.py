#!/usr/bin/env python3
import os
import sys
import argparse
from ningchao.nSys import trick
import ningchao.nSys.excel as xlsKit
import ningchao.nSys.env as envTookit
from collections import defaultdict
example = ''' excel extract by column or raw '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'e', nargs = '+', help = 'excel file and column' )
parser.add_argument( '-n', nargs = '+', help = 'file name and column' )
parser.add_argument( '-c', action = 'store_true', help = 'extract by column' )
parser.add_argument('-o', nargs='?', help ='output file')
parser.add_argument('-merge', action = 'store_true', help ='merge the find object')
parser.add_argument('-eH', action = 'store_true', help ='expression file include header? default no')
parser.add_argument('-igc', action = 'store_true', help = 'ignore case' )
parser.add_argument('-onlyName', nargs= '*', help ='.split pos for the name, eg: 0 0', type = int )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def parse():
    title, tindex = xlsKit.args_deal( args.n )
    exp, eindex = xlsKit.args_deal(args.e)
    onlyName = args.onlyName
    if onlyName == []:
        exit('Please set onlyName like 1 1 as input pos for n and exp\n')
    out = sys.stdout
    if args.o:
        out = open(args.o,'w')
    return title, tindex, exp, eindex, onlyName, out, args.eH, args.igc


def fiter_gene(lst):
    out = []
    lst = [ i for i in lst if i not in ['+','-','']]
    for each in lst:
        try :
            int(each)
        except :
            if not each.startswith('chr'):
                #if '.' in each:
                #    out.append(each.split('.')[0])
                #else :
                    out.append(each)
    return out


def dkey(line, index, igc, onlyName, typ = 'exp'):
    '''title|exp for typ'''
    lst = line.strip().split('\t')
    if onlyName and len(onlyName) == 2 :
        ecol = int( onlyName[1] ) - 1
        tcol = int( onlyName[0] ) - 1
        if typ == 'exp':
            col = ecol
        else :
            col = tcol
        key = '\t'.join([ i.split('.')[col] for i in trick.lst(lst).get(index, returnType = 'list') ] )
    else :
        key = '\t'.join(trick.lst(lst).get(index, returnType = 'list'))
    if igc :
        key = key.upper()
    return key

def deal( title, tindex, exp, eindex, eH, igc, onlyName):
    eInfor, out_lst, efh, tfh = defaultdict( list  ), [], open(exp), open(title)
    if onlyName:
        ecol = int( onlyName[1] ) - 1
        tcol = int( onlyName[0] ) - 1
    if eH:
        efirst = next( efh ).strip()
        out_lst.append( efirst )
    for line in efh:
        line = line.strip()
        key = dkey(line, eindex, igc, onlyName, typ = 'exp')
        eInfor[key].append( line )

    for line in tfh:
        key = dkey(line, tindex, igc, onlyName, typ = 'title')
        if key in eInfor:
            out_lst.extend( eInfor[key] )
        else :
            sys.stderr.write('#!!! %s: not in the %s\n' % (key,exp))
    tfh.close()
    efh.close()
    return out_lst
if __name__ == '__main__':
    title,tindex,exp,eindex,onlyName,out, eH, igc = parse()
    lst = deal( title, tindex, exp, eindex, eH, igc, onlyName )
    for line in lst:
        print ( line.strip('\n'), file = sys.stdout )




