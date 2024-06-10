#!/usr/bin/env python
import sys
import os

def usage():
	if len(sys.argv) <= 1:
		print(sys.argv[0],'file:gtf','length')
		exit()

if __name__ == '__main__':
	usage()
	input_gtf = open(sys.argv[1])
	tss = []
	for line in input_gtf:
		line_arr = line.strip().split('\t')
		if len(line_arr) >= 6 and line_arr[2].startswith('gene'):
			arttributes = line_arr[-1].split(';')
			refseq_gene,gene_biotype,exon_number = '','',''
			for each in arttributes:
				each = each.strip()
				each_arr = iter(each.split("\""))
				for each in each_arr:
					if each.startswith('gene_name'):
						refseq_gene = next(each_arr)
					if each.startswith('gene_biotype'):
						gene_biotype = next(each_arr)
			tss_start = int(line_arr[3]) - int(sys.argv[2])/2
			tss_end = int(line_arr[3]) + int(sys.argv[2])/2
			if tss_start < 0 or tss_end < 0:
				continue
#			if gene_biotype == 'protein_coding':
			value = line_arr[0]+'\t'+str(tss_start)
			if value not in tss:
				print(line_arr[0]+'\t'+str(tss_start)+'\t'+str(tss_end)+'\t'+refseq_gene+'\t'+gene_biotype)
			tss.append(value)








