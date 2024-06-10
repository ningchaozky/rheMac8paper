#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick,status,system
from ningchao.nBio import chromosome

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'wd', nargs='?', help = ' FRMTdist.proper.read.txt ' )
parser.add_argument( '-sample','-s', nargs='?', help = 'sample name' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

def summary( cis, hdir, sample ):
    dit = chromosome.trick([ i for i in fls if 'FRMTdist.proper.read.txt' in i ]).lst2dict()
    chroms = chromosome.chr( 'rh8' ).size()
    tmps,cmds = [],[]
    for chrom in chroms:
        work_dir = cis
        fl = dit[ chrom ]
        out = 'homer.%s.%s.FRMT.comb.summary.txt' % (sample, chrom)
        cmd = '''awk '{if($9=="+"){r1e=$2+1}else{r1e=$3};if($10=="+"){r2e=$5+1}else{r2e=$6};print $7"\t"$1"\t"r1e"\t"$9"\t"$4"\t"r2e"\t"$10}' %s > %s''' % (fl, out)
        tmps.append(os.path.join(work_dir,out))
        cmd = status.cmd( wd = work_dir, flag = '%s.FRMT.comb.summary.txt' % chrom).accessory( cmd, appendix = True )
        cmds.append(cmd)
    out = 'homer.%s.FRMT.comb.summary.txt' % sample
    cmd = 'cat %s > %s' % (' '.join(tmps), out)
    cmd = status.cmd( wd = hdir, flag = 'cat.FRMT.comb.summary.txt').accessory( cmd, appendix = True )
    cmds.append( cmd )
    cmds.extend(system.nSys(tmps).rm(mark = '-rf'))
    return cmds, out


def hommer(hdir, sout, sample):
    cmds = []
    prodir = system.dir(os.path.abspath(args.wd), hdir, '%s_pro' % sample).check()
    cmd = '''makeTagDirectory %s -format HiCsummary %s''' % ( prodir, sout)
    cmd = status.cmd( wd = hdir, flag = 'makeTagDirectory.un.FRMT.comb.summary.txt').accessory( cmd, appendix = True )
    cmds.append( cmd )

    cmd = 'makeTagDirectory %s_pro -update -genome rheMac8 -removePEbg -restrictionSite GATC -both -removeSelfLigation' % sample
    cmd = status.cmd( wd=hdir, flag='makeTagDirectory.pro').accessory( cmd, appendix = True )
    cmds.append( cmd )

    #compartment 
    compart_dir = system.dir(os.path.abspath(args.wd), 'compartment').check()
    cmd = 'runHiCpca.pl %s.100k.400k.comb.CompartmentHomer %s -res 100000 -window 400000 -genome rheMac8 -pc 3 -cpu 8 && mv %s/%s.* %s' % ( sample, prodir, hdir, sample, compart_dir)
    cmd = status.cmd( wd=hdir, flag='runHiCpca.pl.pro').accessory( cmd, appendix = True )
    cmds.append( cmd )
    
    cmd = 'runHiCpca.pl %s.500k.500k.comb.CompartmentHomer %s -res 500000 -window 500000 -genome rheMac8 -pc 3 -cpu 8 && mv %s/%s.* %s' % ( sample, prodir, hdir, sample, compart_dir)
    cmd = status.cmd( wd=hdir, flag='runHiCpca.pl.pro').accessory( cmd, appendix = True )
    cmds.append( cmd )
    
    #loop 
    loop_dir = system.dir(os.path.abspath(args.wd), 'loop').check()
    cmd = 'findTADsAndLoops.pl find %s -res 10000 -window 20000 -genome rheMac8 -minDist 500000 -maxDist 3000000 -cpu 8 && mv %s/%s* %s' % (prodir, prodir, sample, loop_dir)
    cmd = status.cmd( wd=hdir, flag='findTADsAndLoops.pl.pro.10K').accessory( cmd, appendix = True )
    cmds.append( cmd )
    return cmds


if __name__ == '__main__':
    cis = os.path.join(os.path.abspath(args.wd), 'cis')
    hdir = system.dir(os.path.abspath(args.wd), 'homer').check()
    fls, sample = os.listdir( cis ), args.sample
    cmds, sout =  summary( cis, hdir, sample)
    cmds.extend( hommer(hdir, sout, sample) )
    print('\n'.join( cmds ))




















