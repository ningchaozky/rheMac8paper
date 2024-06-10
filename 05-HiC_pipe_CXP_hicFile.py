#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,system,status, fix
from ningchao.nBio import chromosome
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('work_dir', nargs='?', help ='reference')
parser.add_argument('-res', nargs='?', type = int, help ='reference', default = 40000 )
parser.add_argument('-sample', nargs='?', help ='sample for anasysis, 4M.PFC.hic.rep1', required = True )
parser.add_argument('-find_dir', nargs='+', help ='find dir append work_dir for files' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def txt(wd, chroms):
    #fls = [ i for i in os.listdir(os.path.join( wd, 'cis')) if 'FRMTdist.proper.read.txt' in i ]
    #dit = chromosome.trick(fls).lst2dict()
    outfls,cmds = [],[]
    FRMTdist_proper_already, FRMTdist_proper_need = system.dir(args.work_dir).file_check('FRMTdist.proper.read.txt'.format(args.res//1000), depth = 10)
    if args.find_dir:
        for fd in args.find_dir:
            FRMTdist_proper_already_2, FRMTdist_proper_need_2 = system.dir( fd ).file_check('FRMTdist.proper.read.txt'.format(args.res//1000), depth = 10)
            FRMTdist_proper_already.update( FRMTdist_proper_already_2 )
            FRMTdist_proper_need.update( FRMTdist_proper_need_2 )
    for chrom in FRMTdist_proper_already:
        work_dir = os.path.join(os.path.abspath(wd), 'juicer')
        #fl = system.dir(os.path.join( os.path.abspath(wd),'cis', dit[chrom])).check()
        #outfl = system.dir(os.path.join( os.path.abspath(wd),'juicer','cis.{}.{}.juicer.merged_nodups.txt'.format(chrom, sample))).check()
        fl = FRMTdist_proper_already[chrom]
        outfl_dir = system.dir(os.path.join( os.path.abspath(wd),'juicer')).check()
        outfl = os.path.join( outfl_dir, 'cis.{}.{}.juicer.merged_nodups.txt'.format(chrom, sample))
        outfls.append(outfl)
        cmd = '''awk '{split($1,a,"r");split($4,b,"r");split($17,fidx,"[_|F]");if($9=="+"){s1=0;}else{s1=16;};if($10=="+"){s2=0;}else{s2=16;};print $7" "s1" chr"a[2]" "$3" "fidx[2]-1" "s2" chr"b[2]" "$6" "fidx[4]-1" "$8" "$8}' %s > %s''' % (fl, outfl)
        cmd = status.cmd( wd =work_dir, flag='chrom.nodups.%s' % chrom).accessory( cmd )
        cmds.append(cmd)
    txt_ouput = 'cis.%s.chrall.juicer.merged_nodups.txt' % sample
    cmd = 'cat %s > %s' % (' '.join(outfls), txt_ouput)
    cmd = status.cmd( wd = work_dir, flag='cat.chrallnodups').accessory( cmd )
    cmds.append(cmd)
    return cmds,txt_ouput
def t2bin( wd, txt_output, q = 30 ):
    juicer = 'java -jar /home/soft/soft/juicer/CPU/scripts/juicer_tools.1.7.6_jcuda.0.8.jar'
    #juicer_Mbol = '/home/ningch/data/genome/rheMac8/Mbol/rh8_MboI.txt'
    juicer_Mbol = '/home/soft/data/genome/rheMac8/hic/rheMac8_Mbol.nochr.txt'
    work_dir = system.dir(os.path.join(os.path.abspath(wd), 'juicer')).check()
    genome = fix.fix( chromosome.chr('rh8').abspath ).append('nochr')
    hic = 'cis.%s.chrall.juicer.merged_nodups_%d.hic' % (sample, q)
    cmd = '''sed -i 's/chr//g' {}'''.format(txt_output)
    cmd1 = status.cmd( wd=work_dir, flag='{}'.format('chrdel.forHic')).accessory( cmd )
    cmd = '''%s pre -q %d -f %s %s  %s %s -r 25000000,1000000,500000,250000,100000,50000,40000,25000,20000,10000,5000''' % (juicer, q, juicer_Mbol, txt_output, hic, genome)
    cmd2 = status.cmd( wd=work_dir, flag='hic.%d.chrallnodups' % q).accessory( cmd )
    return '\n'.join( [ cmd1, cmd2] )
    


if __name__ == '__main__':
    chroms = chromosome.chr('rh8').chr
    sample, wd = args.sample, args.work_dir
    cmds, txt_output = txt(wd, chroms)
    print('\n'.join(cmds))
    print(t2bin( wd, txt_output))
























