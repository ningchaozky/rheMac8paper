#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from pybedtools.bedtool import BedTool
from ningchao.nSys import trick,fix
from ningchao.nBio import rheMac,bed
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help ='reference')
parser.add_argument('-tss', nargs='?', help ='reference', default = '/home/ningch/data/genome/rheMac8/exon/gene.sort.tss1K')
parser.add_argument('-qlst', nargs = '+', default = [ 10, 5, 3 ], help = 'qvalue for select the region. can use multi')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def check(wd, keys, qlst, tab):
    infor = {}
    if os.path.exists(os.path.join(wd,tab)):
        tfh = open(tab)
        header = tfh.next().strip().split('\t')
        for each in tfh:
            each_arr = each.strip().split('\t')
            dit = dict( list(zip(header, each_arr)) )
            for key in keys:
                if key in dit:
                    val = float(dit[key])
                    for cut_off in qlst:
                        if val >= cut_off :
                            trick.set2dict(infor, cut_off, key, [])
                            infor[cut_off][ key ].append('\t'.join([dit['chr'], dit['start'], dit['end']]))
        tfh.close()
    return merge(infor)


def merge(infor):
    dit = {}
    if infor :
        for cut_off in infor:
            for peirod_rep in infor[ cut_off ]:
                #bed = '\n'.join( infor[ cut_off ][ peirod_rep ] )
                #slow and use mine method
                #bed_obj = BedTool( bed, from_string=True )
                bed_merge = bed.merge(infor[ cut_off ][ peirod_rep ])
                trick.set2dict(dit, cut_off, peirod_rep, '\n'.join(bed_merge))
                
    return dit


def format(line, prks, dit):
    out = line.strip().split('\t')
    if dit:
        for prk in prks:
            peirod, rep, cut_off = prk.split('.')
            key = '.'.join([peirod, rep])
            cut_off = int(cut_off)
            if cut_off in dit:
                if key in dit[cut_off]:
                    out.append(dit[cut_off][key].replace('\t',',').replace('\n','|'))
            else :
                    out.append('')
        return '\t'.join(out)
    else :
        out.extend(['' for i in prks])
        return '\t'.join(out)

def title(qlst):
    header = ['chrom', 'start', 'end', 'name', 'place_hold', 'chain']
    keys = rheMac.rheMac().peirods_combine_reps()
    prks = []
    for cut_off in qlst:
        prks.extend(['.'.join([i,str(cut_off)]) for i in keys ])
    header.extend(prks)
    return header, prks, keys


if __name__ == '__main__':
    tss,wd,qlst = args.tss, args.dir, args.qlst
    header, prks, keys = title(qlst)
    tss_fh = open(tss)
    print('\t'.join(header))
    for line in tss_fh:
        chrom, start, end, name, place_hold, chain = line.strip().split('\t')
        prefix = fix.fix('.'.join([name, chrom, start, end, chain]))
        tab = prefix.append('tab')
        interaction = check(wd, keys, qlst, tab)
        print(format(line, prks, interaction))

























