#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import argparse
from ningchao.nSys import trick,path,status, system, fix
from collections import defaultdict
from ningchao.nBio import chromosome
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('rep_dirs', nargs='+', help ='replicate dirs')
parser.add_argument('-peirod','-p', choices = ['E50','E80','E90','E120','0M','4M','45Y','20Y'], help = 'peirod for choices', required = True )
parser.add_argument('-region', choices = ['PFC'], help = 'region', default = 'PFC' )
parser.add_argument('-res', choices = [ 40000, 10000 ], type = int, help = 'res for choices', default = 40000 )
parser.add_argument('-rutine', choices = [ 'bedpe', 'raw_matrix' ], nargs = '*', help = 'for anasyis type', default = [ 'bedpe', 'raw_matrix' ])
parser.add_argument('-wd', help = 'work_dir', default = '.')
parser.add_argument('-c', nargs = '*', help = 'chroms for analysis' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

######################################################################################
##  生成cis/trans BedPE文件  ##


def dirs( **kwargs ):
    work_dir, res, rep = [ kwargs.get(i) for i in ('work_dir','res', 'rep') ]
    dit = {}
    out_dir = path.mkdir( work_dir, rep, make = True )
    trick.dinit(dit, 'root', out_dir)
    trick.dinit(dit, 'work_dir', work_dir)
    trick.dinit(dit, 'cis', path.mkdir( out_dir, 'cis', make = True ))
    trick.dinit(dit, 'trans', path.mkdir( out_dir, 'trans', make = True ))
    trick.dinit(dit, 'raw', path.mkdir( out_dir, 'raw_matrix', make = True ))
    trick.dinit(dit, 'raw.ICE', path.mkdir( out_dir, 'raw_matrix','ICE', make = True ))
    trick.dinit(dit, 'raw.ICE.res', path.mkdir( out_dir, 'raw_matrix','ICE',res, make = True ))
    trick.dinit(dit, 'Norm', path.mkdir( out_dir, 'Norm_matrix', make = True ))
    trick.dinit(dit, 'Norm.ICE', path.mkdir( out_dir, 'Norm_matrix', 'ICE', make = True ))
    trick.dinit(dit, 'Norm.ICE.res', path.mkdir( out_dir, 'Norm_matrix', 'ICE', res, make = True ))
    trick.dinit(dit, 'domain', path.mkdir( out_dir, 'domain', make = True ))
    trick.dinit(dit, 'domain.ICE', path.mkdir( out_dir, 'domain', 'ICE', make = True ))
    trick.dinit(dit, 'domain.ICE.ScaleTo1', path.mkdir( out_dir, 'domain', 'ScaleTo1', make = True ))
    return dit


def get_rep_name( dirname ):
    if 'R1' in dirname or 'r1' in dirname or 'rep1' in dirname:
        return 'rep1'
    if 'R2' in dirname or 'r2' in dirname or 'rep2' in dirname:
        return 'rep2'

def dirs_chrom_find( dirs, **kwargs):
    fls = defaultdict( str )
    for fdir in dirs :
        for tfl in list( system.dir(fdir).fls('{}'.format( kwargs.get('pattern') ), depth = 1) ):
            chrom = [i for i in tfl.split('.') if 'chr' in i ][0]
            fls[chrom] = tfl
    return fls


def bedpe( **kwargs ):
#def bedpe(rep_dirs, peirod, region, rep, fragments_dir, **kwargs):
    rep_dirs, peirod, region, rep, fragments_dir = [ kwargs.get(i) for i in ('rep_dirs','peirod','region','rep','fragments_dir') ]
    cmds,rep_dirs,fls = [],[ os.path.abspath(i) for i in rep_dirs ],[]
    #for typ in 'cis','trans':
    #cmd = '''cat %s/cis.E90.PFC.rep1.hicup.dedup.sort.bedpe %s/cis.E90.PFC.rep2.hicup.dedup.sort.bedpe > %s''' % ()
    already_chroms = dirs_chrom_find( rep_dirs, pattern = 'FRMTdist.proper.read.txt' )
    for chrom in chr_obj.chr:
        if chrom in already_chroms :
            continue
        if 'chrX' in chrom :
            continue
        if kwargs.get('c') and chrom not in kwargs.get('c'):
            continue
        ifl = ' '.join( dirs_chrom_find( rep_dirs, pattern = 'hicup.dedup.sort.bedpe' ) )
        if not ifl:
            exit('check {}.hicup.dedup.sort.bedpe not find in {}'.format(chrom, rep_dirs))
        chrom_sort_bedpe = 'cis.%s.%s.%s.%s.hicup.dedup.sort.bedpe' % (peirod, region, rep, chrom)
        cmd = '''cat %s | awk -v chr="%s" '$1==chr{print}' - | sort -k1,1 -k2,2n - > %s''' % ( ifl, chrom, chrom_sort_bedpe)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = 'cis.%s.sort.bedpe' % chrom, echo = 'generate %s ' % chrom_sort_bedpe )
        cmds.append(cmd)
        
        reBed = os.path.join(fragments_dir, 'rheMac8.%s.Mbol.REfragment.bed' % chrom)
        end5 = 'cis.%s.%s.hicup.dedup.5end.sort.bedpe' % (chrom, peirod)
        assignment = 'cis.%s.%s.hicup.fragment.assignment.txt' % (chrom, peirod)
        read1_assignment = 'cis.%s.%s.hicup.read1.Fassignment.txt' % (chrom, peirod)
        read2_assignment = 'cis.%s.%s.hicup.read2.Fassignment.txt' % (chrom, peirod)
        FRMTdist_proper_read = 'cis.%s.%s.hicup.FRMTdist.proper.read.txt' % (chrom, peirod)
        cmd = '''awk 'BEGIN{OFS="\\t"}{if($9=="+"){R1e=$2+1}else{R1e=$3};R1s=R1e-1;if($10=="+"){R2e=$5+1}else{R2e=$6};R2s=R2e-1;print $1"\\t"R1s"\\t"R1e"\\t"$4"\\t"R2s"\\t"R2e"\\t"$7"\\t"$8"\\t"$9"\\t"$10}' %s > %s''' % ( chrom_sort_bedpe, end5)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = end5, echo = 'generate %s' % end5)
        cmds.append(cmd)

        cmd = '''bedtools pairtobed -f 1 -s -type both -a %s -b %s > %s''' % (end5, reBed, assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = assignment, echo = 'generate %s' % assignment)
        cmds.append(cmd)
        
        cmd = '''awk "%s == 1 {print}"  %s > %s''' % ('NR % 2', assignment, read1_assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = read1_assignment, echo = 'generate %s' % read1_assignment)
        cmds.append(cmd)

        cmd = '''awk '%s == 0 {print}' %s > %s''' % ('NR % 2', assignment, read2_assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = read2_assignment, echo = 'generate %s' % read2_assignment)
        cmds.append(cmd)
        
        cmd = '''paste %s %s | awk '$14!=$31{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13"\\t"$27"\\t"$28"\\t"$29"\\t"$14"_"$30}' - | awk '! ( $2 == 0 || $5 == 0) { print }' > %s''' % (read1_assignment, read2_assignment, FRMTdist_proper_read)
        cmd = status.cmd_accessory( cmd, work_dir = wds['cis'], flag = FRMTdist_proper_read, echo = 'generate %s' % FRMTdist_proper_read)
        cmds.append(cmd)
    return cmds


def raw_matrix( **kwargs ):
#def raw_matrix(fragments_dir, ver, dit):
    fragments_dir, ver, dit, peirod, res = [ kwargs.get(i) for i in ('fragments_dir','ver','dit','peirod', 'res') ]
    genome_file, resolution = [ kwargs.get(i) for i in ('genome_file','resolution') ]
    cmds = []
    for chrom in chr_obj.schr:
        if kwargs.get('c') and chrom not in kwargs.get('c'):
            continue
        #files generate 
        gbin = os.path.join(fragments_dir,'{}.{}.{}.bin.txt'.format( ver, chrom, res))
        reBed = os.path.join(fragments_dir, '%s.%s.%s.REfragment.bed' % ( ver, chrom, kwargs.get('enzyme')) )
        gdist = fix.fix(gbin).insert('distribution')
        #os.path.join(fragments_dir, '%s.%s.bin.distribution.txt' % ( ver, chrom))
        proper = os.path.join(dit['cis'], 'cis.%s.%s.hicup.FRMTdist.proper.read.txt' % (chrom, peirod))
        chr_bin = os.path.join(dit['root'], '{}.{}.txt'.format(chrom, res) )
        reads_count = os.path.join(dit['raw.ICE.res'], 'clean.%s.%s.%s.%s.matrix.txt' % (peirod, chrom, res, ver))
        ice_raw_matrix = os.path.join(dit['raw.ICE.res'], 'HTC.%s.%s.%s.%s.matrix.txt' % (peirod, chrom, res, ver))
        #### 生成 固定size bin ###
        cmd = '''bedtools makewindows -g %s -w %s | awk -v chr="%s" '$1==chr{print $0}' - | awk '{print $0"\\tbin"NR}' - > %s''' % (genome_file, resolution, chrom, os.path.basename(gbin))
        cmd = status.cmd_accessory( cmd, work_dir = fragments_dir, flag = gbin, echo = 'generate %s' % gbin)
        cmds.append(cmd)

        cmd = '''awk '$6=="+"{print $0}' %s | awk '{mid=int(($2+$3)/2);mids=mid-1;print $1"\\t"mids"\\t"mid"\\t"$4"\\t"$5"\\t"$6}' - | bedtools intersect -f 1 -wa -wb -a - -b %s |grep -v "-1" - > %s''' % (reBed, gbin, os.path.basename(gdist))
        cmd = status.cmd_accessory( cmd, work_dir = fragments_dir, flag = os.path.basename(gdist), echo = 'generate: {}'.format(gdist))
        cmds.append(cmd)

        ### 将reads分到指定bin里
        cmd = '''awk -F "\\t" 'NR==FNR{a[$4]=$10;next;}{split($17,b,"_");print a[b[1]]"\\t"a[b[2]]"\\t"$0}' %s %s | awk '{gsub("bin",""); print $0}' - | sort -k1,1n -k2,2n - | bedtools groupby -g 1,2 -c 9 -o count_distinct | awk '{print "bin"$1"\\tbin"$2"\\t"$3}' - > %s''' % (gdist, proper, chr_bin)
        cmd = status.cmd_accessory( cmd, work_dir = dit['root'], flag=os.path.basename(chr_bin), echo = 'generate: {}'.format(chr_bin))
        cmds.append(cmd)
        
        ### reads count计数 ###
        cmd = '''num=$(wc -l %s) && generate_plain.pl %s $num | awk '{ gsub("bin","");print }' - | sort -k1,1n -k2,2n - | bin2matrix.pl - >  %s''' % (gbin,chr_bin, os.path.basename(reads_count))
        cmd = status.cmd_accessory( cmd, work_dir = dit['raw.ICE.res'], flag=os.path.basename(reads_count), echo = 'generate %s' % reads_count)
        cmds.append(cmd)
        
        ### 生成HiTC兼容的ICE raw natrix矩阵  ###
        cmd = '''awk -v resolution=%d '{header=(NR-1)*resolution; print header"\\t"$0}' %s | HTCmatrix.pl - %s %d %s %s > %s''' % ( resolution, reads_count, chrom, resolution, ver, genome_file, ice_raw_matrix)
        cmd = status.cmd_accessory( cmd, work_dir = dit['raw.ICE.res'], flag=os.path.basename(ice_raw_matrix), echo = 'generate %s' % ice_raw_matrix)
        cmds.append(cmd)
    return cmds

if __name__ == '__main__':
    #rutines = args.rutine
    #rep_dirs, peirod, rep, resolution = args.rep_dirs, args.peirod, 'merge', args.res
    #ver, enzyme, res, chr_obj = 'rheMac8', 'Mbol', str(resolution//1000)+'K', chromosome.chr('rh8')
    #genome_file,region = chr_obj.abspath, args.region
    #perl_path = '/home/soft/soft/HiC-CXP-pipeline/pl'
    #fragments_dir = '/home/soft/data/genome/rheMac8/hic'
    #ref, work_dir = chr_obj.fa, os.path.abspath(args.wd)
    chr_obj = chromosome.chr('rh8')
    kwargs = { **vars(args), 'res': str( args.res//1000 )+'K', 'chr_obj': chromosome.chr('rh8'), 'ver': 'rheMac8', 'enzyme': 'Mbol'}
    kwargs.update({'genome_file': chr_obj.abspath, 'region': args.region, 'perl_path': '/home/soft/soft/HiC-CXP-pipeline/pl'})
    kwargs.update({'fragments_dir': '/home/soft/data/genome/rheMac8/hic', 'ref': chr_obj.fa, 'work_dir': os.path.abspath(args.wd)})
    kwargs.update({'rep': 'merge', 'resolution': args.res })
    #build dirs 
    #wds = dirs( kwargs.get('work_dir'), kwargs.get('res'), **kwargs)
    kwargs.update( {'dit': dirs( **kwargs ) } )
    
    #bedpe 
    if 'bedpe' in kwargs.get('rutine'):
        #print('\n'.join( bedpe( rep_dirs, peirod, region, rep, fragments_dir, **kwargs)))
        print('\n'.join( bedpe( **kwargs )))
    if 'raw_matrix' in kwargs.get('rutine'):
    #generate raw matrix
        #print('\n'.join(raw_matrix(fragments_dir, ver, wds, **kwargs)))
        print('\n'.join(raw_matrix(  **kwargs )))

