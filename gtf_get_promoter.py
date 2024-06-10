#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import itertools
#from gtfparse import read_gtf
from ningchao.nSys import trick
from collections import defaultdict

example = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description='%s %s' % (os.path.basename(sys.argv[0]), example), formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'gtf', nargs = '?', help = 'gtf file' )
parser.add_argument( '-s', nargs = '?', help = 'span for start end end', default = 3000, type = int )
parser.add_argument( '-bu', nargs = '?', help = 'span for gene Body up span ', default = 5000, type = int )
parser.add_argument( '-bd', nargs = '?', help = 'span for gene Body down span ', default = 5000, type = int )
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()

mean, saf = open(args.gtf +'.mean','w'), open(args.gtf +'.SAF','w')
promoter, gene2transcript, geneTransPromoter, geneBody = open(args.gtf +'.promoter','w'), open(args.gtf +'.gene2transcripts','w'), open(args.gtf +'.geneTransPromoterMean','w'), open(args.gtf +'.geneBody.updown5K','w')
chroms = [ *[ str(i) for i in range( 1, 25 )], 'X', 'Y' ]
chroms = [ *[ 'chr{}'.format(i) for i in  chroms ], * chroms ]
chroms.extend(['X','3L','2L','3R','2R','4','Y'])

def gtf_transcripts_parse( gtf ):
    infor = defaultdict( lambda : defaultdict( list  )  )
    with open( args.gtf ) as f :
        for line in f :
            line_arr = line.strip().split('\t')
            if line_arr[0] not in chroms :
                continue
            if 'exon' != line_arr[2] :
                continue
            tmp = [[ j for j in i.strip().split(' ') if j.strip() ] for i in line_arr[-1].split(';') if i.strip() ]
            tmp = { k:v.strip('\"') for k,v in tmp }
            line_arr[3], line_arr[4] = int( line_arr[3]), int( line_arr[4])
            tmp.update( {'strand': line_arr[6], 'chrom': line_arr[0], 'start': line_arr[3], 'end': line_arr[4], 'len': abs( line_arr[4] - line_arr[3]) })
            if 'transcript_id' in tmp:
                infor[tmp['gene_id']][tmp['transcript_id']].append( tmp )
    return infor
if __name__ == '__main__':
    infor = gtf_transcripts_parse( args.gtf )
    for gene_id in infor:
        length, trans, symbols = [],[], []
        for transcript_id in infor[gene_id]:
            transcript_len = 0
            starts_ends, chrom, strand = [], '', ''
            for exon in infor[gene_id][transcript_id]:
                #exon_infor = infor[gene_id][transcript_id][exon]
                transcript_len += abs(exon['start'] - exon['end'])
                print ( gene_id, exon['chrom'].replace('chr',''), exon['start'], exon['end'], exon['strand'], sep = '\t', file = saf )
                starts_ends.extend([exon['start'], exon['end']])
                strand, chrom = exon['strand'], exon['chrom']
                if 'gene_name' in exon:
                    symbols.append(exon['gene_name'])
                #if 'transcript_name' in exon:
                #    symbols.append(exon['transcript_name'])
            starts_ends = sorted( starts_ends )
            if strand == '-' :
                pstart = starts_ends[-1] - args.s
                pend = starts_ends[-1] + args.s
            elif strand == '+':
                pstart = starts_ends[0] - args.s
                pend = starts_ends[0] + args.s
            
            body_start = starts_ends[0] - args.bu
            body_end = starts_ends[-1] + args.bd
            print ( chrom.replace('chr',''), pstart, pend, transcript_id, '0', strand, sep = '\t', file = promoter )
            print ( chrom.replace('chr',''), starts_ends[0] - 5000, starts_ends[-1] + 5000, transcript_id, '0', strand, sep = '\t', file = geneBody )
            print ( gene_id, transcript_id, strand, sep = '\t', file = gene2transcript )
            chrom = chrom.replace('chr','')
            #tss start end, body +up +down, and body starts ends
            trans.append ( [ chrom, pstart, pend, transcript_id, '0', strand, gene_id, chrom, starts_ends[0], starts_ends[-1], chrom, body_start, body_end ] )
            length.append( transcript_len )
        for etrans in trans :
            if symbols:
                symbol = list( set(symbols) )[0]
            else :
                symbol = etrans[6]
            print ( * etrans,  symbol, int(sum(length) / len(length)), sep = '\t', file = geneTransPromoter)
        print ( gene_id, int(sum(length) / len(length)), sep = '\t', file = mean)

mean.close()
saf.close()
promoter.close()
















exit()
raw_df = 'ncbiRefseq.xls'
if not os.path.exists( raw_df ):
    df = read_gtf( args.gtf )
    df_genes = df[df["feature"] == "transcript"]
    df.to_csv( raw_df, sep = '\t', index = None)
else :
    with open( raw_df ) as f :
        header = next(f).strip().split('\t')
        for line in f:
            line_arr = line.strip().split('\t')
            tmp = dict(zip(header,line_arr))
            if tmp['feature'] == 'transcripts':
                print ( tmp )


transcript = 'transcript.xls'
if not os.path.exists( transcript ):
    df = read_gtf( args.gtf )
    df_genes = df[df["feature"] == "transcript"]
    df_transcript = df[df["feature"] == "transcript"]
    df_transcript.to_csv( transcript, sep = '\t', index = None)


if 0 :
    promoter = 'promoter.{}.xls'.format(args.s)
    if not os.path.exists( promoter ):
        df = read_gtf( args.gtf )
        df_genes = df[df["feature"] == "transcript"]
        df.loc[df["strand"] == '-', 'promoter_start'] = df['end'] - 1000
        df.loc[df["strand"] == '+', 'promoter_start'] = df['start'] - 1000 
        df.loc[df["strand"] == '-', 'promoter_end'] = df['end'] + 1000 
        df.loc[df["strand"] == '+', 'promoter_end'] = df['start'] + 1000
        df = df.astype({'promoter_start': 'int32', 'promoter_end': 'int32'})
        df.to_csv( promoter, sep = '\t', index = None)




















