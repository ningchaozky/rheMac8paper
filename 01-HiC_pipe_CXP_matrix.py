#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import argparse
from ningchao.nSys import trick,path,status,fix
from ningchao.nBio import chromosome
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('final_bam', nargs='?', help ='final bam')
parser.add_argument('-peirod','-p', choices = ['E50','E80','E90','E120','0M','4M','45Y','20Y'], help = 'peirod for choices', required = True )
parser.add_argument('-region', choices = ['PFC'], help = 'region', default = 'PFC' )
parser.add_argument('-rep','-r', choices = ['rep1','rep2','merge'], help = 'rep for choices', required = True )
parser.add_argument('-res', choices = [ 40000, 10000 ], type = int, help = 'res for choices', default = 40000 )
parser.add_argument('-rutine', choices = [ 'bedpe', 'raw_matrix' ], nargs = '*', help = 'for anasyis type', default = [ 'bedpe', 'raw_matrix' ])
parser.add_argument('-wd', help = 'work_dir', default = '.')
parser.add_argument('-chrs','-c', nargs = '*', help = 'chroms for analysis')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

######################################################################################
##  生成cis/trans BedPE文件  ##


def dirs( **kwargs ):
    work_dir, res, rep = [ kwargs.get(i) for i in ('wd','res_dir','rep') ]
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


#def bedpe( bam, region, rep, fragments_dir, **kwargs):
def bedpe( **kwargs ):
    bam, cmds = os.path.abspath( kwargs.get('final_bam') ), []
    region, rep, fragments_dir = [ kwargs.get(i) for i in ('region','rep','fragments_dir') ]
    ref, peirod, wds = [ kwargs.get(i) for i in ('ref','peirod', 'wds') ]
    #for typ in 'cis','trans':
    for typ in 'cis',:
        output = '%s.%s.%s.%s.hicup.dedup.sort.bedpe' % (typ, peirod,region,rep)
        if typ == 'cis':
            cmd = '''samtools view %s | grep "CT:Z:FAR" | samtools view -bS -T %s - | bedtools bamtobed -bedpe -mate1 -i - |awk '$8>=10{print}' -|awk '$1==$4{print}' - > %s''' % (bam, ref, output)
        elif typ == 'trans':
            cmd = '''samtools view %s | grep "CT:Z:FAR" | samtools view -bS -T %s - | bedtools bamtobed -bedpe -mate1 -i - |awk '$8>=10{print}' -|awk '$1!=$4{print}' - > %s''' % (bam, ref, output)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = '{}.bedpe'.format(typ), echo = 'generate %s ' % (output), appendix = True )
        if typ == 'trans':
            cmd = cmd.replace('CT:Z:FAR','CT:Z:TRANS').replace('$1==$4','$1!=$4')
        cmds.append(cmd)
    output = '%s.%s.%s.%s.hicup.dedup.sort.bedpe' % ( typ, peirod,region, rep)
    for chrom in kwargs.get('chr_obj').chr:
        if kwargs.get('chrs') and chrom not in kwargs.get('chrs'):
            continue
        chrom_sort_bedpe = '%s.%s.%s.%s.%s.hicup.dedup.sort.bedpe' % ( typ, peirod, region, rep, chrom)
        cmd = '''awk -v chr="%s" '$1==chr{print}' %s | sort -k1,1 -k2,2n - > %s''' % ( chrom, output, chrom_sort_bedpe)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = '{}.{}.sort.bedpe'.format( typ, chrom), echo = 'generate %s ' % chrom_sort_bedpe, appendix = True )
        cmds.append(cmd)
        
        reBed = os.path.join(fragments_dir, 'rheMac8.%s.Mbol.REfragment.bed' % chrom)
        end5 = '%s.%s.%s.hicup.dedup.5end.sort.bedpe' % ( typ, chrom, peirod)
        assignment = '%s.%s.%s.hicup.fragment.assignment.txt' % ( typ, chrom, peirod)
        read1_assignment = '%s.%s.%s.hicup.read1.Fassignment.txt' % ( typ, chrom, peirod)
        read2_assignment = '%s.%s.%s.hicup.read2.Fassignment.txt' % ( typ, chrom, peirod)
        FRMTdist_proper_read = '%s.%s.%s.hicup.FRMTdist.proper.read.txt' % ( typ, chrom, peirod)
        cmd = '''awk 'BEGIN{OFS="\\t"}{if($9=="+"){R1e=$2+1}else{R1e=$3};R1s=R1e-1;if($10=="+"){R2e=$5+1}else{R2e=$6};R2s=R2e-1;print $1"\\t"R1s"\\t"R1e"\\t"$4"\\t"R2s"\\t"R2e"\\t"$7"\\t"$8"\\t"$9"\\t"$10}' %s > %s''' % ( chrom_sort_bedpe, end5)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = end5, echo = 'generate %s' % end5, appendix = True )
        cmds.append(cmd)

        cmd = '''bedtools pairtobed -f 1 -s -type both -a %s -b %s > %s''' % (end5, reBed, assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = assignment, echo = 'generate %s' % assignment, appendix = True )
        cmds.append(cmd)
        
        cmd = '''awk "%s == 1 {print}"  %s > %s''' % ('NR % 2', assignment, read1_assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = read1_assignment, echo = 'generate %s' % read1_assignment, appendix = True )
        cmds.append(cmd)

        cmd = '''awk '%s == 0 {print}' %s > %s''' % ('NR % 2', assignment, read2_assignment)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = read2_assignment, echo = 'generate %s' % read2_assignment, appendix = True )
        cmds.append(cmd)
        
        cmd = '''paste %s %s | awk '$14!=$31{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13"\\t"$27"\\t"$28"\\t"$29"\\t"$14"_"$30}' - | awk '! ( $2 == 0 || $5 == 0) { print }' > %s''' % (read1_assignment, read2_assignment, FRMTdist_proper_read)
        cmd = status.cmd_accessory( cmd, work_dir = wds[typ], flag = FRMTdist_proper_read, echo = 'generate %s' % FRMTdist_proper_read, appendix = True )
        cmds.append(cmd)
    return cmds


def raw_matrix( **kwargs ):
    fragments_dir, ver, dit, enzyme, peirod = [ kwargs.get(i) for i in ('fragments_dir','ver', 'wds', 'enzyme', 'peirod') ]
    resolution, binKb,genome_file = [ kwargs.get(i) for i in ('res','res_dir','genome_file') ]
    cmds = []
    for chrom in kwargs.get('chr_obj').chr:
        if kwargs.get('chrs') and chrom not in kwargs.get('chrs'):
            continue
        binKb = kwargs.get('res_dir')
        gbin = os.path.join(fragments_dir,'%s.%s.%s.bin.txt' % ( ver, chrom, binKb))
        reBed = os.path.join(fragments_dir, '%s.%s.%s.REfragment.bed' % ( ver, chrom, enzyme))
        gdist = os.path.join(fragments_dir, '%s.%s.%s.bin.distribution.txt' % ( ver, chrom, binKb))
        proper = os.path.join(dit['cis'], 'cis.%s.%s.hicup.FRMTdist.proper.read.txt' % (chrom, peirod))
        chr_bin = os.path.join(dit['root'], '%s.txt' % chrom)
        reads_count = os.path.join(dit['raw.ICE.res'], 'clean.%s.%s.%s.%s.matrix.txt' % (peirod, chrom, ver, binKb))
        ice_raw_matrix = os.path.join(dit['raw.ICE.res'], 'HTC.%s.%s.%s.%s.matrix.txt' % (peirod, chrom, ver, binKb))
        #### 生成 固定size bin ###
        cmd = '''bedtools makewindows -g %s -w %s | awk -v chr="%s" '$1==chr{print $0}' - | awk '{print $0"\\tbin"NR}' - > %s''' % (genome_file, resolution, chrom, os.path.basename(gbin))
        cmd = status.cmd_accessory( cmd, work_dir = fragments_dir, flag = gbin, echo = 'generate %s' % gbin, appendix = True )
        cmds.append(cmd)
        cmd = '''awk '$6=="+"{print $0}' %s | awk '{mid=int(($2+$3)/2);mids=mid-1;print $1"\\t"mids"\\t"mid"\\t"$4"\\t"$5"\\t"$6}' - | bedtools intersect -f 1 -wa -wb -a - -b %s |grep -v "-1" - > %s''' % (reBed, gbin, os.path.basename(gdist))
        cmd = status.cmd_accessory( cmd, work_dir = fragments_dir, flag = os.path.basename(gdist), echo = 'generate %s' % os.path.basename(gdist), appendix = True )
        cmds.append(cmd)

        ### 将reads分到指定bin里
        cmd = '''awk -F "\\t" 'NR==FNR{a[$4]=$10;next;}{split($17,b,"_");print a[b[1]]"\\t"a[b[2]]"\\t"$0}' %s %s | awk '{gsub("bin",""); print $0}' - | sort -k1,1n -k2,2n - | bedtools groupby -g 1,2 -c 9 -o count_distinct | awk '{print "bin"$1"\\tbin"$2"\\t"$3}' - > %s''' % (gdist, proper, chr_bin)
        cmd = status.cmd_accessory( cmd, work_dir = dit['root'], flag = os.path.basename(chr_bin), echo = 'generate %s' % chr_bin, appendix = True )
        cmds.append(cmd)
        
        ### reads count计数 ###
        cmd = '''num=$(wc -l %s) && generate_plain.pl %s $num | awk '{ gsub("bin","");print }' - | sort -k1,1n -k2,2n - | bin2matrix.pl - >  %s''' % (gbin,chr_bin, os.path.basename(reads_count))
        cmd = status.cmd_accessory( cmd, work_dir = dit['raw.ICE.res'], flag = os.path.basename(reads_count), echo = 'generate %s' % reads_count, appendix = True)
        cmds.append(cmd)
        
        ### 生成HiTC兼容的ICE raw natrix矩阵  ###
        cmd = '''awk -v resolution=%d '{header=(NR-1)*resolution; print header"\\t"$0}' %s | HTCmatrix.pl - %s %d %s %s > %s''' % ( resolution, reads_count, chrom, resolution, ver, genome_file, ice_raw_matrix)
        cmd = status.cmd_accessory( cmd, work_dir = dit['raw.ICE.res'], flag = os.path.basename(ice_raw_matrix), echo = 'generate %s' % ice_raw_matrix, appendix = True )
        cmds.append(cmd)
    return cmds

if __name__ == '__main__':
    kwargs = {'perl_path':'/home/soft/soft/HiC-CXP-pipeline/pl', 'fragments_dir': '/home/soft/data/genome/rheMac8/hic'}
    kwargs.update( ** vars( args ) )
    kwargs.update({ 'wd': os.path.abspath(args.wd), 'chr_obj': chromosome.chr('rh8') })
    kwargs.update({'res_dir': str(args.res//1000)+'K', 'ver': 'rheMac8', 'enzyme': 'Mbol', 'genome_file': kwargs.get('chr_obj').abspath})
    kwargs.update({'ref': kwargs.get('chr_obj').fa, 'wds': dirs( **kwargs )})
    #bedpe 
    print ( kwargs, file = sys.stderr )
    if 'bedpe' in kwargs.get('rutine'):
        print('\n'.join(bedpe( **kwargs ) ) )
    if 'raw_matrix' in kwargs.get('rutine'):
    #generate raw matrix
        print('\n'.join(raw_matrix( **kwargs ) ) )

    exit()
    #rutines = args.rutine
    #bam,peirod,rep,resolution = args.final_bam, args.peirod, args.rep, args.res
    #ver, enzyme, res, chr_obj, kwargs = 'rheMac8', 'Mbol', str(resolution/1000).replace('.0','')+'K', chromosome.chr('rh8')
    #perl_path = '/home/ningch/soft/HiC-CXP-pipeline/pl'
    #perl_path = '/home/soft/soft/HiC-CXP-pipeline/pl'
    #fragments_dir = '/home/soft/data/genome/rheMac8/hic'
    #fragments_dir = '/home/ningch/data/genome/rheMac8/Mbol/rheMac8_Mbol/rheMac8_Mbol_fragments/'
    #ref, work_dir = chr_obj.fa, os.path.abspath(args.wd)
    #build dirs 
    #wds = dirs(work_dir, res)
