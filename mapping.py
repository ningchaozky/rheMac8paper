#!/usr/bin/env python3
import os
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
from ningchao.nSys import trick, system, status, fix
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'dir', nargs = '?', help = 'dir for file trim output file')
parser.add_argument( '-i', nargs = '?', help = 'index for mapping', default = 'hg38')
parser.add_argument( '-d', nargs = '?', help = 'depth for find reads pair', default = 0, type = int)
parser.add_argument( '-s', choices = ['STAR','bwa','bowtie'], help = 'mapping software', required = True)
parser.add_argument( '-f', nargs = '*', help = 'fiter out dir', default = [])
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def CONF( *args, **kwargs):
    index = {'mm10': '/home/soft/data/genome/mm10/Mus_musculus.GRCm38.chr.fa', 'hg19': '/home/soft/data/genome/human/hg19/human.hg19.chrall.fa',
    'hg38': '/dataE/smb_share/data/genome/human/hg38/ncbi_dataset/data/GCF_000001405.26/hg38.fa'}
    print ('Choose genome:', kwargs.get('i'), index.get(kwargs.get('i')), file = sys.stderr )
    conf = {'index': index.get(kwargs.get('i')) }
    return conf

def reads_pair( args ):
    infor = defaultdict( lambda : defaultdict( list ) )
    reads = system.dir( args.dir ).fls(r'clean.fq.gz$|P\.fq\.gz$', depth = args.d, abspath = True)
    for read in reads:
        rdir = os.path.dirname( read )
        infor[rdir]['reads'].append( read )
        infor[rdir]['reads'].append( next(reads) )
        #infor[rdir]['prefix'] = os.path.commonprefix([ infor[rdir]['reads'][0], infor[rdir]['reads'][1]] ).strip('.').strip('_').replace('.fq.gz','')
        print (infor[rdir]['reads'][0], infor[rdir]['reads'][1] )
        infor[rdir]['prefix'] = trick.lst([ infor[rdir]['reads'][0], infor[rdir]['reads'][1]]).commonprefix()
    for wd in infor:
        infor[wd]['reads'] = sorted( infor[wd]['reads'] )
    return infor
def single_read():
    infor = defaultdict( lambda : defaultdict( list ) )
    reads = system.dir( args.dir ).fls('U.fq.gz$', depth = args.d, abspath = True)
    return reads
def srun( cmd ):
    for line in system.run( cmd, shell = True, noReturn = True):
        print ( line )

def multi_run( cmds ):
    with mp.Pool(3) as p:
        p.map( srun, cmds )
    p.close()
    p.join()

def check_sam( wd, sam):
    index = True
    sam = os.path.basename( sam )
    for e in os.listdir(wd ):
        if e.endswith('sam') and os.path.getsize(os.path.join(wd,e)) > 100:
            print ( e, ' file\n' )
        if e == sam:
            print ( 'exist: for already', sam )
            break
            exit()
    return index

if __name__ == '__main__':
    reads = reads_pair( args  )
    sreads = single_read()
    #fls = list(system.dir( args.dir ).fls('P\..*fq.gz$|P\.fq\.gz$', depth = args.d))
    #print ( 'Find reads:', list(fls), file = sys.stderr )
    conf = CONF( ** vars(args) )
    cmds = []
    if args.s == 'STAR':
        for wd in reads:
            prefix = os.path.basename(wd)
            #cmd = '''STAR --readFilesCommand zcat --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 6 --genomeDir /home/soft/data/genome/human/hg38/STAR --alignIntronMin 20 --alignIntronMax 500000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outSAMattrRGline ID:plus1 SM:plus2 PL:ILLUMINA --outFilterMismatchNmax 2 --outSAMmultNmax 2 --outSAMmapqUnique 60  --outFileNamePrefix {}_result --alignMatesGapMax 1000000 --outSJfilterReads Unique --readFilesIn {} {}'''.format( prefix, *reads[wd] )
            cmd = '''STAR --readFilesCommand zcat --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 6 --genomeDir /home/soft/data/genome/human/hg38/STAR --alignIntronMin 20 --alignIntronMax 500000 --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:plus1 SM:plus2 PL:ILLUMINA --outFilterMismatchNmax 2 --outSAMmultNmax 2 --outSAMmapqUnique 60  --outFileNamePrefix {}_result --alignMatesGapMax 1000000 --outSJfilterReads Unique --readFilesIn {} {}'''.format( prefix, *reads[wd]['reads'] )
            print ( status.cmd_accessory( cmd, wd = wd, flag = 'mapping', appendix = True ) )
    elif args.s == 'bwa':
        #fls = system.dir( args.dir ).fls('P\.|\.fastq')
        for sread in sreads:
            wd = os.path.dirname( sread )
            SAM_prefix = os.path.basename(sread).replace('U.fq.gz','')
            sam = SAM_prefix + '.sam'
            if check_sam( wd, sam ):
                cmd1 = '''/home/soft/soft/bwa/v0.7.17/bwa mem -M -R "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:ILLUMINA" {} {} -t 8  > {}\n'''.format( SAM_prefix, SAM_prefix, SAM_prefix, conf.get('index'), sread, sam )
                cmd1 = status.cmd( cmd1, wd = wd).accessory( flag = 'bwamapping_{}'.format(SAM_prefix) )
                cmds.append( cmd1 )
        for wd in reads:
            read_1, read_2 = reads[wd]['reads']
            prefix = reads[wd]['prefix']
            #prefix = os.path.basename( read_1 ).replace( '_1P.fq.gz', '')
            #/home/soft/data/genome/mm10/Mus_musculus.GRCm38.chr.fa
            SAM = os.path.basename(prefix)
            LB, ID = SAM, SAM.split('_')[0]
            #cmd1 = '''/home/soft/soft/bwa/v0.7.17/bwa mem -M -R "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:ILLUMINA" {} {} -t 8  > {}_1.sam\n'''.format( ID, SAM, LB, conf.get('index'), read_1, prefix )
            #cmd2 = '''/home/soft/soft/bwa/v0.7.17/bwa mem -M -R "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:ILLUMINA" {} {} -t 8  > {}_2.sam\n'''.format( ID, SAM, LB, conf.get('index'), read_2, prefix )
            sam_out = prefix + '.bwa.sam'
            index = check_sam( wd, sam_out )
            print ( index )
            cmd = '''/home/soft/soft/bwa/v0.7.17/bwa mem -M -R "@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:ILLUMINA" {} {} {} -t 8  > {}\n'''.format( ID, SAM, LB, conf.get('index'), read_1, read_2, sam_out )
            for fit in args.f:
                if fit in cmd:
                    index = False
                    break 
            if index :
                cmd3 = status.cmd( cmd, work_dir = wd).accessory( flag = 'bwamapping' )
                cmds.append( cmd3 )
                #cmd1 = status.cmd( cmd1, wd = wd).accessory( flag = 'bwamapping_read_1' )
                #cmds.append( cmd1 )
                #cmd2 = status.cmd( cmd2, wd = wd).accessory( flag = 'bwamapping_read_2' )
                #cmds.append( cmd2 )

    elif args.s == 'bowtie':
        for wd in reads:
            read_1, read_2 = reads[wd]['reads']
            cmd = 'bowtie2 -p 6 -q --phred33 --very-sensitive --end-to-end --no-unal --no-mixed --no-discordant -N 1 -X 3000000 -x {} -1 {} -2 {} -S {}.sam'.format( conf.get('index'), read_1, read_2, reads[wd]['prefix'] )
            cmd = status.cmd_accessory( cmd, wd = wd, flag = 'snpmapping' )
            cmds.append( cmd )
    
    multi_run(cmds)
































