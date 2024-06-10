#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
from ningchao.nSys import trick, system, fix, parse, status
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('dir', nargs='?', help = ' MarkDuplicates picard')
parser.add_argument('-d', nargs='?', help = 'depth for sam find', default = 10, type = int )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

fls = system.dir( os.getcwd() ).fls('sam$', depth = args.d)
for fl in fls:
    wd = os.path.dirname(fl)
    stats = parse.ini( os.path.join( wd, 'readme.status.txt' ) ).to_dict()
    sbam = fix.fix(fl).change('sorted.bam')
    mbam = fix.fix(sbam).change('marked_duplicates.bam')
    bw = fix.fix(mbam).change('bw')
    print ( 'ulimit -u 50000')
    cmd = '''/home/soft/soft/jdk/v21/bin/java -jar /home/soft/soft/picard.jar SortSam I={} O={} SORT_ORDER=coordinate'''.format(fl, sbam)
    cmd = status.cmd_accessory( cmd, wd = wd, flag = 'sortCoordinate', appendix = True )
    print ( cmd )
    cmd = '''/home/soft/soft/jdk/v21/bin/java -jar /home/soft/soft/picard.jar MarkDuplicates I={} O={} M=marked_dup_metrics.txt REMOVE_DUPLICATES=True'''.format( sbam, mbam )
    cmd = status.cmd_accessory( cmd, wd = wd, flag = 'MarkDuplicates', appendix = True )
    print ( cmd )
    cmd = '''samtools index {}'''.format(mbam)
    cmd = status.cmd_accessory( cmd, wd = wd, flag = 'bamindex', appendix = True )
    print ( cmd )
    cmd = '''bamCoverage -b {} -bs 10 -o {} --normalizeUsing CPM'''.format(mbam, bw)
    cmd = status.cmd_accessory( cmd, wd = wd, flag = 'bam2bw', appendix = True )
    print ( cmd )



























