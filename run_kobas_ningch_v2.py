#!/usr/bin/env python3
import os
import re
import sys
import sqlite3
import argparse
from ningchao.nSys import trick, excel, fix

example = ''' '''
parser = argparse.ArgumentParser( prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter )
parser.add_argument( 'fl', nargs = '?', help = 'file input excel,index_col,4' )
parser.add_argument( '-t', choices = [ 'name', 'transcripts'], help = 'gene name or transcripts name', default = 'name')
parser.add_argument( '-index_split_name','-isn', nargs = '*', help = 'indexl split and get name . 0 mean use . split and the first place', default = ['.', 0 ])
parser.add_argument( '-cluster','-c', nargs = '?', help = 'key words for select gene. Eg: c1$', default = '.*')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


hsa_db_file = '/home/soft/soft/kobas/v3.03/sqlite3/hsa.db'
con = sqlite3.connect( hsa_db_file )
cur = con.cursor()


matrix, eargs = excel.parse( args.fl )
index_col = eargs['index_col']


with open( matrix ) as f :
    uniq_genes = []
    trans_fl_name = '{}.kobas.id'.format( os.path.basename(matrix) )
    cluster_file_name = 'all' if args.cluster == '.*' else args.cluster.replace('$','').replace('^','')
    trans_fl_name_raw = open( '{}.kobas.{}.id.raw'.format( os.path.basename(matrix), cluster_file_name ), 'w')
    if args.cluster:
        trans_fl_name = fix.fix(trans_fl_name).insert( args.cluster.replace('$','').replace('.','' ).replace('*','all'))
    idfh, uniq = open( trans_fl_name, 'w'), []
    for gene in f:
        if args.cluster:
            p = re.compile(r'{}'.format(args.cluster))
            if not p.search(gene):
                continue
        line_arr = gene.strip().split('\t')
        if len( line_arr ) <= index_col:
            print ( '{} small than {}'.format(line_arr, index_col), file = sys.stderr )
            continue
        gene = line_arr[index_col]
        other_cols = [ v for i,v in enumerate(line_arr) if i != index_col ]
        if args.index_split_name:
            gene = gene.split(args.index_split_name[0])[args.index_split_name[1]]
            if gene not in uniq_genes:
                print ( 'Match {} to /home/soft/soft/kobas/v3.03/sqlite3/hsa.db hsa database'.format(gene), file = sys.stderr )
            uniq_genes.append(gene)
        gene = gene.strip()
        for hsaid, genes in cur.execute('''SELECT * FROM Genes where name like ? ''', tuple( [ '%'+gene+'%'] ) ) :
            gos = list(cur.execute( '''select * from GeneGos where gid = ?''', tuple([hsaid])))
            print ( hsaid.replace('hsa:',''), file = idfh )
            print ( gene, *other_cols, hsaid.replace('hsa:',''), file = trans_fl_name_raw, sep = '\t' )
            #for gid,goid in gos:
            #    desc = list(cur.execute( '''select * from Gos where goid = ?''', tuple([goid])))
            #    for each in desc:
                    #print ( hsaid, gene, genes, each, file = sys.stderr )
            #        if hsaid not in uniq :
            #            print ( hsaid.replace('hsa:',''), file = idfh )
            #        uniq.append( hsaid )
    idfh.close()
cmd = '''run_kobas.py -i {} -o {}.go -s hsa -t id:ncbigene'''.format( idfh.name, idfh.name )
print ( cmd )
cmd = '''run_kobas_ningch_return.py {}'''.format( idfh.name.replace('.kobas.all.id','' ) )
print ( cmd )
cmd = '''goBar.py {}.go -n 1000 -k '.*' -s 8 150'''.format( idfh.name )
print ( cmd )
print ('run_kobas_ningch_go_kegg_pathway.py . ')




























