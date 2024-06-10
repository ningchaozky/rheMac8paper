#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
from ningchao.nSys import trick, status, system
from ningchao.nBio import chromosome,rheMac
example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('sample_work_dir', nargs = '?', help ='rep1')
parser.add_argument('-res', nargs = '?', help ='resolution, default 40000', type = int, default = 40000 )
parser.add_argument('-chrs', nargs = '*' )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def outdir_get_name(dir_name):
    dirs = dir_name.split('/')
    name = '.'.join([dirs[5], dirs[6], dirs[8]]).replace('prefrontal','PFC')
    return name
def ICE_Normal( need, **kwargs ):
    cmds = []
    HTC_already, HTC_need = system.dir( args.sample_work_dir, depth = 10).file_check( 'HTC.*{}'.format(kwargs.get('res_dir_basename')), echo = '01-HiC_pipe_CXP_matrix.py results: HTC need')
    name_prefix = rheMac.trick().name( kwargs.get('normal_dir') )
    name_prefix = '.'.join([name_prefix, kwargs.get('res_dir_basename')])
    for chrom in need:
        if kwargs.get('chrs') and chrom not in kwargs.get('chrs'):
            continue
        if chrom not in HTC_already:
            print ( '#Warning: {} not find 01 step'.format( chrom, HTC_already ) )
            continue
        dirname, basenmae = os.path.split( HTC_already[chrom] )
        cmd = 'HiCT_ICE_Normalization.r %s --output_dir %s -p %s --name %s' % ( dirname, kwargs.get('normal_dir'), basenmae, name_prefix)
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('normal_dir'), flag = 'HiCT_ICE_Normalization.r.{}'.format(chrom), force = ['all'])
        cmds.append( cmd )
    return cmds

def clean_ice( normal_matrix_already, need, **kwargs):
    raw, cmds, fai = normal_matrix_already, [], chromosome.chr('rh8').fai
    #fls = [ i for i in os.listdir(normal_dir) if name in i and i.endswith('ICE_matrix.rout.mat') ]
    cmds.append('#clean for {}'.format( need ) )
    for chrom in raw:
        if kwargs.get('chrx') and chrom not in kwargs.get('chrx'):
            continue
        if chrom not in need:
            continue
        #if chrom in kwargs.get('clean_already'):
        #    continue
        name = rheMac.trick().name( raw[chrom] )
        normal_matrix = '%s.%s.%s.ICE_matrix' % ( name, chrom, str( kwargs.get('res')//1000  ) + 'K')
        clean_matrix = 'clean.%s.%s.%s.ICE_matrix' % (name, chrom, str( kwargs.get('res')//1000 ) + 'K')
        cmd = 'norm2DImatrix.pl %s %s %s > %s' % ( raw[chrom], chrom, fai, normal_matrix)
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('normal_dir'), flag='norm2DImatrix.pl.{}'.format( chrom ), force = ['all'])
        cmds.append(cmd)
        cmd = 'cut -f 4- %s > %s' % (normal_matrix, clean_matrix)
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('normal_dir'), flag = 'cutnormal2clean.%s' % chrom, force = ['all'])
        cmds.append(cmd)
    return cmds

def normal_depth( clean_already, depth_normal_matrix_need, **kwargs ):
    cmds = []
    for chrom in clean_already:
        if chrom not in kwargs.get('chrs'):
            continue
        name = rheMac.trick().name( clean_already[chrom] )
        cmd = 'HiCT_depth_Normalization.r %s --sample %s --res %s' % ( os.path.dirname(clean_already[chrom]), chrom, str( kwargs.get('res')//1000 ) + 'K')
        cmd = status.cmd_accessory( cmd, work_dir = kwargs.get('normal_dir'), flag = 'HiCT_depth_Normalization.r.{}'.format(chrom), force = ['all'])
        cmds.append(cmd)
    return cmds

def check( pattern, **kwargs):
    all_chroms = chromosome.chr('rh8').chr
    find_alreay = list( system.dir(args.sample_work_dir).fls( pattern, depth = 10, abspath = True ) )
    find_already_chroms = trick.string(find_alreay).get_chrom( basename = True )
    need_chroms = set(all_chroms.keys()) - set( find_already_chroms )
    if need_chroms:
        print ('#Warning!!!: pattern {} find, {} {}'.format( pattern, kwargs.get('echo'), need_chroms), file = sys.stderr )
    return dict(zip(find_already_chroms, find_alreay)), need_chroms

if __name__ == '__main__':
    kwargs, work_dir, res_dir_basename = vars(args), os.path.abspath(args.sample_work_dir), str( args.res//1000   ) + 'K'
    kwargs.update({'res_dir_basename': res_dir_basename, 'work_dir': work_dir})
    kwargs.update({'normal_dir': os.path.join( work_dir,'Norm_matrix/ICE/', res_dir_basename)})
    normal_matrix_already, normal_matrix_need = system.dir( args.sample_work_dir, depth = 10 ).file_check( '.*{}.*ICE_matrix.rout.mat'.format(kwargs.get('res_dir_basename')), echo = 'Find ICE_matrix.rout.mat')
    clean_already, clean_need = system.dir( args.sample_work_dir, depth = 10  ).file_check( '{}.*ICE_matrix$|{}.*ICE_matrix\.txt$$'.format( res_dir_basename, res_dir_basename), echo = '02: clean need')
    depth_normal_matrix_already, depth_normal_matrix_need = system.dir(args.sample_work_dir, depth = 10).file_check('{}.*scaleTo1'.format(res_dir_basename), echo = '03: depth normal to 1w scaleTo1.1w.txt')
    #find already normal
    #HiTC need chroms normal the data
    #HiTC normal the data
    normal_matrix_need.add('chrX')
    cmds = ICE_Normal( normal_matrix_need, **kwargs )
    print( '\n'.join( cmds ) )

    #clean ice need chroms
    print('\n'.join(clean_ice( normal_matrix_already, clean_need | depth_normal_matrix_need | set(['chrX']), clean_already = clean_already, **kwargs)))


    print('\n'.join(normal_depth( clean_already, depth_normal_matrix_need | set(['chrX']), **vars(args))))




















