#!/usr/bin/env python
import os
import sys
import argparse
import ningchao.nBio.fasta as faTookit
import ningchao.nSys.env as envTookit
import ningchao.nSys.trick as trickTookit

parser = argparse.ArgumentParser(prog=sys.argv[0],description='')
parser.add_argument('-gtf', nargs='?', help ='gtf download from public', default = 'Macaca_mulatta.Mmul_8.0.1.90.gtf')
parser.add_argument('-cds', nargs='?', default = 'Macaca_mulatta.Mmul_8.0.1.cdna.all.fa', help ='all cds download from public database')
parser.add_argument('-w', nargs='?', default = '.', help ='key word')
parser.add_argument('-o', nargs='?', help ='output')
if len(sys.argv) == 1:
	parser.print_help().__str__
	sys.exit(2)
args = parser.parse_args()
out = open(args.o,'w')


Infor = {}
if 1 :
	fh = open(args.gtf)
	for line in fh:
		line_arr = line.split('\t')
		if len(line_arr) > 7:
			if line_arr[2] == 'transcript':
				transcript_id = trickTookit.list2dict(line_arr[8].split(';')[2].split(' '))['transcript_id'].replace("\"","")
				gene_id = trickTookit.list2dict(line_arr[8].split(';')[0].split(' '))['gene_id'].replace("\"","")
				trickTookit.set1dict(Infor,transcript_id,gene_id)


cds = faTookit.fa2dict(args.cds)
gene_infor = {}
for each in cds:
	tid = (cds[each].name).split('.')[0]
	ver = each
	trickTookit.set2dict(gene_infor,Infor[tid],tid,{})
	gene_infor[Infor[tid]][tid].update({'seq':str(cds[each].seq)})
	gene_infor[Infor[tid]][tid].update({'len':len(cds[each].seq)})
	gene_infor[Infor[tid]][tid].update({'ver':ver})

for gene in gene_infor:
	seq,seq_len,tid,ver = '',0,'',''
	for each in gene_infor[gene]:
		if gene_infor[gene][each]['len'] > seq_len:
			seq_len = gene_infor[gene][each]['len']
			seq = gene_infor[gene][each]['seq']
			ver = gene_infor[gene][each]['ver']
			tid = each
	out.write('>%s %s %s %s\n' % (ver,gene,seq_len,tid))
	out.write(faTookit.format(seq)+'\n')







