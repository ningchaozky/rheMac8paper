#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,status,fix,system
from ningchao.nBio import chromosome
from collections import defaultdict

example = ''' xuepengs DI and TAD calling '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('work_dir', nargs='?', help =' dir for the first like rep1 rep2 merge')
parser.add_argument('-res','-r', nargs='?', help ='resultion for the domain calling', default = 40000, type = int )
parser.add_argument('-sample','-s', nargs='?', help ='sample name for di and tad calling. like, E50.PFC.hic.rep2', required = True )
parser.add_argument('-wss', action='store_true', help = 'wss results to continue')
parser.add_argument('-chrs', nargs='*')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def DI( **kwargs ):
    #pattern = 'clean.*\.{}[k|K]\..*ICE.scaleTo1.1w.txt'.format( kwargs.get('res')//1000)
    #scaleTo1_already, scaleTo1_need = system.dir( args.work_dir ).file_check( pattern, depth = 10 )
    window = kwargs.get('w') and kwargs.get('w') or 2000000
    chrom_obj, res = kwargs.get('chrom_obj'), args.res
    chrs, fai, genome, cmds = kwargs.get('chroms'),chrom_obj.fai,chrom_obj.abspath,[]
    gfiles = defaultdict( lambda : defaultdict ( str ) )
    for chrom in kwargs.get('chroms'):
        gfiles['raw_DI_name'][chrom] = 'DI.%s.clean.%s.%s.ICE.scaleTo1.1w.txt' % ( kwargs.get('sample'), chrom, kwargs.get('rsBinK'))
        gfiles['di_scale'][chrom] = '%s.%s.%s.norm_DI.scaleTo1.1w.txt' % ( kwargs.get('sample'), chrom, kwargs.get('rsBinK'))
        gfiles['raw_di_bdg'][chrom] = '%s.%s.%s.norm_DI.scaleTo1.1w.bdg' % ( kwargs.get('sample'), chrom, kwargs.get('rsBinK'))
    chrom_all_scale, chrom_all_bdg = [],[]
    #step 02 get resulst
    normal_scaleTo1_already, normal_scaleTo1_need = system.dir(args.work_dir, depth = 10).file_check('{}.*scaleTo1.1w.txt'.format(kwargs.get('rsBinK')))
    di_clean_already, di_clean_need = system.dir(args.work_dir, depth = 10).file_check('clean.*{}.ICE.scaleTo1.1w.txt'.format(kwargs.get('rsBinK')))
    for chrom in kwargs.get('chroms'):
        if chrom in ['chrY','chrMT','chrM']:
            continue
        cmd = '''awk -v bin=%s '{header=(NR-1)*bin;print header"\\t"$0}' %s | normdepth2DImatrix.pl - %s %s rheMac8 %s > %s''' % ( args.res, normal_scaleTo1_already[chrom],chrom, args.res, genome, gfiles['raw_DI_name'][chrom])
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.DI.unadj', force = ['all'])
        cmds.append(cmd)

    scaleTo1_1w_txt_already, scaleTo1_1w_txt_need = system.dir(args.work_dir).file_check('clean.*{}[k|K].*scaleTo1.1w.txt'.format(args.res//1000), depth = 10)
    scaleTo1_scale_already, scaleTo1_scale_need = system.dir(args.work_dir).file_check('{}[k|K].*norm_DI.scaleTo1.1w.txt'.format(args.res//1000), depth = 10)
    for chrom in kwargs.get('chroms'):
        if kwargs.get('chrs') and chrom not in kwargs.get('chrs'):
            continue
        if chrom not in scaleTo1_1w_txt_already:
            continue
        cmd = '''DI_from_matrix.pl %s %d %d %s |grep -v "inf" > %s ''' % ( gfiles['raw_DI_name'][chrom], res, window, fai, gfiles['di_scale'][chrom])
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + 'DI_from_matrix.pl', force = ['all'] )
        cmds.append(cmd)


    scaleTo1_1w_bed_already, scaleTo1_1w_bed_need = system.dir(args.work_dir).file_check('{}[k|K].*scaleTo1.1w.bdg'.format(args.res//1000), depth = 10)
    for chrom in chrs:
        #if kwargs.get('chrs') and chrom not in kwargs.get('chrs'):
        #    continue
        if chrom not in scaleTo1_1w_txt_already:
            continue
        raw_di_bdg = scaleTo1_1w_txt_already[chrom]
        cmd = '''awk '{print "chr"$0}' %s > %s''' % ( gfiles['di_scale'][chrom], gfiles['raw_di_bdg'][chrom])
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.DI.bdg', force = ['all'])
        cmds.append(cmd)

        cmds.append('\n\n\n')
    #chrom_all_scale_already, chrom_all_scale_need = system.dir( kwargs.get('work_dir'), depth = 10).file_check('norm_DI.scaleTo1.1w.txt')
    chrom_all_scale_file = '%s.chrall.DI.scaleTo1.1w.txt' % args.sample
    cmd = 'cat %s > %s' % (' '.join(gfiles['di_scale'].values()), chrom_all_scale_file)
    cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = 'DI.cat.txt', force = ['all'])
    cmds.append(cmd)
    
    chrom_all_scale_file_matlab = '%s.chrall.DI.scaleTo1.1w.matlab.txt' % args.sample
    cmd = '''perl -pe  's/^chr//; s/^X/23/; s/^Y/24/' %s > %s''' % (chrom_all_scale_file, chrom_all_scale_file_matlab)
    cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.matlabneedXto23', force = ['all'])
    cmds.append(cmd)

    chrom_all_bdg_file = chrom_all_scale_file.replace('txt','bdg')
    cmd = 'cat %s > %s' % (' '.join(gfiles['raw_di_bdg'].values()), chrom_all_bdg_file)
    cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.DI.cat.bdg', force = ['all'] )
    cmds.append(cmd)
    
    chrom_all_bdg_file_sort = fix.fix(chrom_all_bdg_file).insert('sort')
    cmd = 'sort -k1,1 -k2,2n %s -o %s' % (chrom_all_bdg_file, chrom_all_bdg_file_sort)
    cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.DI.sort.bdg', force = ['all'])
    cmds.append(cmd)
    
    di_bw = chrom_all_scale_file.replace('txt','bw')
    cmd = 'bedGraphToBigWig %s %s %s' % ( chrom_all_bdg_file_sort, genome, di_bw)
    cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('dir_unadj'), flag = chrom + '.DI.bw', force = ['all'])
    cmds.append(cmd)
    return cmds

def domain( **kwargs ):
    wd, sample, chrom_obj, tmin, prob, binsize = [ kwargs.get(i) for i in ['work_dir', 'sample', 'chrom_obj', 'tmin', 'prob','binsize'] ]
    hmm_dir = system.dir( kwargs.get('dir_unadj'), 'HMM%s' % tmin).check()
    hmm_chrom = system.dir( kwargs.get('dir_unadj'), 'HMM%s' % tmin, 'chr').check()
    fai = chrom_obj.fai
    chrom_all_scale_file_matlab = '%s.chrall.DI.scaleTo1.1w.matlab.txt' % sample
    chrom_all_scale_hmm_file = '%s.chrall.hmmout.scaleTo1.1w.txt' % sample
    fh = open('/home/soft/soft/HiC-CXP-pipeline/matlab/HMM_calls.m')
    mfile = os.path.join(hmm_dir,'HMM_calls.m')
    mfh = open(mfile,'w')
    for line in fh:
        if '/home/chenxp/3T3/mouse/Sperm/R1/domain/ICE/ScaleTo1/' in line:
            if line.startswith('filename'):
                line = '''filename = '%s';\n''' % os.path.join( os.path.abspath(os.path.dirname(hmm_dir)), chrom_all_scale_file_matlab)
            if line.startswith('file =  \'/home/chenxp/3T3/mouse/Sperm/'):
                line = '''file = '%s';\n''' % os.path.join( os.path.abspath(os.path.dirname(hmm_dir)), chrom_all_scale_hmm_file)
        mfh.write(line)
    mfh.close()
    cmd = '''matlab -nodesktop -nosplash -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r HMM_calls'''
    cmd = status.cmd_accessory( cmd, work_dir = hmm_dir, flag = 'runmatlab', force = ['all'])
    return cmd

def summary( **kwargs ):
    dir_unadj, chrom_unadj_dir, sample, chrom_obj, tmin = [ kwargs.get(i) for i in ['dir_unadj','chrom_unadj_dir','sample', 'chrom_obj','tmin'] ]
    tmin, prob, binsize = [ kwargs.get(i) for i in ['tmin','prob','binsize'] ]
    cmds, beds, fai = [],[], chrom_obj.fai
    chrom_all_scale_hmm_file = '%s.chrall.hmmout.scaleTo1.1w.txt' % sample
    chrom_all_scale_file = '%s.chrall.DI.scaleTo1.1w.txt' % sample
    norm_7colfile = '%s.ICE_7colfile.txt' % sample
    cmd = '''1file_ends_cleaner.pl %s %s | 2converter_7col.pl > %s''' % (chrom_all_scale_hmm_file,chrom_all_scale_file,norm_7colfile)
    cmd = status.cmd_accessory( cmd, work_dir = dir_unadj, flag = 'norm_7colfile.generate', appendix = True, force = ['all'])
    cmds.append(cmd)

    norm_7colfile = os.path.join(os.path.abspath(dir_unadj), norm_7colfile)

    for chrom in kwargs.get('chroms'):
        chrom_7colfile = '%s.%s.norm_7colfile.txt' % (sample, chrom)
        cmd = '''awk -v chr=%s '$1==chr{print $0}' %s > %s''' % (chrom, norm_7colfile, chrom_7colfile)
        cmd = status.cmd_accessory( cmd, work_dir = chrom_unadj_dir, flag = '%s.7colfile.generate' % chrom, force = ['all'])
        cmds.append(cmd)

        norm_domain_bed = '%s.%s.norm_domain.bed' % (sample,chrom)
        cmd = '''4hmm_probablity_correcter.pl %s %d %f %d | 5hmm-state_caller.pl %s %s | awk 'NR>=2{print}' - | 6hmm-state_domains.pl - | awk 'NF==3{print}' -> %s''' % (chrom_7colfile, tmin, prob, binsize, fai, chrom, norm_domain_bed)
        cmd = status.cmd_accessory( cmd, work_dir = chrom_unadj_dir, flag = '%s.norm_domain_bed' % chrom, force = ['all'])
        cmds.append(cmd)
        beds.append(norm_domain_bed)
    chrom_all_bed = '%s.chrall.ICE.domain.bed' % sample
    cmd = 'cat %s > %s' % (' '.join(beds), chrom_all_bed)
    cmd = status.cmd_accessory( cmd, work_dir = chrom_unadj_dir, flag = 'cat.chrom.7colfile', force = ['all'])
    cmds.append(cmd)
    
    cmd = 'wc -l %s > %s.TAD.summary.txt' % (chrom_all_bed, sample)
    cmd = status.cmd_accessory( cmd, work_dir = chrom_unadj_dir, flag = 'TAD.summary.txt.generate', force = ['all'])
    cmds.append(cmd)
    return cmds

if __name__ == '__main__':
    kwargs = vars(args)
    kwargs.update({ 'chrom_obj': chromosome.chr('rh8'), 'tmin': 2, 'prob': 0.99, 'binsize': args.res, 'rsBinK': str(args.res//1000) + 'K'} )
    #tmin,prob,binsize = 2,0.99,40000
    #chrs,fai,genome = chrom_obj.chr,chrom_obj.fai,chrom_obj.abspath
    if args.chrs :
        kwargs.update({'chroms': args.chrs})
    else :
        kwargs.update({'chroms': { i:v for i,v in chromosome.chr('rh8').chr.items() if i != 'chrX' }} )
    kwargs.update( {'scaleTo1_dir': system.dir( args.work_dir,'Norm_matrix/ICE/%s' % kwargs.get('rsBinK')).check()} )
    kwargs.update( {'domain_dir': system.dir( args.work_dir, 'domain/ICE/ScaleTo1/', kwargs.get('rsBinK')).check()} )
    kwargs.update( {'dir_unadj': system.dir( kwargs.get('domain_dir'), 'unadj').check()} )
    kwargs.update( {'chrom_unadj_dir': system.dir( kwargs.get('dir_unadj'), 'HMM%d' % kwargs.get('tmin'), 'chr').check()} )
    #domain_dir = dir_check(work_dir,'domain/ICE/ScaleTo1/', resk)
    #dir_unadj = dir_check(work_dir,'domain/ICE/ScaleTo1/%s/unadj' % resk )
    #chrom_unadj_dir = dir_check(work_dir,'domain/ICE/ScaleTo1/%s/unadj' % resk, 'HMM%d' % tmin, 'chr')
    #sys.stderr.write(scaleTo1_dir + '\n')

    #DI value get
    print('\n'.join(DI( **kwargs )) )

    #HMM calling
    print(domain( **kwargs ))

    #bed and TAD summary
    print('\n'.join(summary( **kwargs )))

























