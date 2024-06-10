#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nBio import chromosome,rheMac
from ningchao.nSys import trick
from collections import defaultdict
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'final_bam', nargs= '+', help ='final bam')
parser.add_argument( '-peirod','-p', choices = ['E50','E80','E90','E120','0M','4M','45Y','20Y'], help = 'peirod for choices' )
parser.add_argument( '-region', choices = ['PFC'], help = 'region', default = 'PFC' )
parser.add_argument( '-rep','-r', choices = ['merge'], help = 'rep for choices', default = 'merge')
parser.add_argument( '-res', choices = [ 40000, 10000 ], type = int, help = 'res for choices', default = 40000 )
parser.add_argument( '-rutine','-ru', choices = [ '01', '02', '03', '05', '06'], nargs = '*', help = 'for anasyis type', default = [ '01', '02', '03', '05', '06' ] )
parser.add_argument( '-wd', help = 'work_dir', default = '.' )
parser.add_argument( '-tmp', action = 'store_true' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def tmps(wd):
    out = [ os.path.join(wd,i) for i in os.listdir(wd) if 'Fassignment' in i ]
    out.extend( [ os.path.join(wd,i) for i in os.listdir(wd) if i.endswith('bedpe') ] )
    out.extend( [ os.path.join(wd,i) for i in os.listdir(wd) if i.endswith('assignment.txt') ] )
    return out


if __name__ == '__main__':
    cmds, bams, peirod, rep, resolution, rutines = defaultdict( str ), args.final_bam, args.peirod, args.rep, args.res, args.rutine
    ver, enzyme, res, chr_obj = 'rheMac8', 'Mbol', str(resolution/1000)+'k', chromosome.chr('rh8')
    #genome_file,region = chr_obj.genome(), args.region
    region = args.region
    perl_path = '/home/soft/soft/HiC-CXP-pipeline/pl'
    fragments_dir = '/home/soft/data/genome/rheMac8/hic'
    #fragments_dir = '/home/ningch/data/genome/rheMac8/Mbol/rheMac8_Mbol/rheMac8_Mbol_fragments/'
    #ref, work_dir = chr_obj.fa(), os.path.abspath(args.wd)
    work_dir = os.path.abspath( args.wd )
    peirods = rheMac.rheMac().peirods()
    if len(bams) == 1 :
        bam = bams[0]
        bam_dir = os.path.dirname(os.path.abspath(bam))
        if not peirod or not rep:
            peirod, rep = chromosome.trick( bam ).peirod_rep()
        if args.tmp:
            cis = os.path.join(bam_dir, rep, 'cis')
            print('\n'.join(tmps(cis)))
            exit( 1 )
        cmds['01'] = '01-HiC_pipe_CXP_matrix.py %s -peirod %s -rep %s -res %d' % (bam, peirod, rep, resolution)
    else :
        peirod = [ i for i in os.path.abspath('.').replace('4month','4M').split('/') if i in peirods ][0]
        trick.set1dict(cmds, '01', '01-HiC_pipe_CXP_matrix_combine.py rep1 rep2 -peirod %s -res %d' % (peirod, resolution))
    cmds['02'] = '02-HiC_pipe_CXP_ICE.py %s -res %d' % (rep, resolution)
    cmds['03'] = '03-HiC_pipe_CXP_HMM_TAD.py %s -res %d -s %s' % (rep, resolution, '.'.join([peirod,rep]))
    cmds['05'] = '05-HiC_pipe_CXP_hicFile.py %s -res %d -s %s' % (rep, resolution, '.'.join([peirod,rep]))
    cmds['06'] = '06-HiC_pipe_CXP_loop_compartment.py %s -s %s' % (rep, '.'.join([peirod,rep]))

    for each in rutines:
        print(cmds[each])
