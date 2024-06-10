#!/usr/bin/env Rscript
#first
alen <- commandArgs()
suppressMessages(library(argparser))
suppressMessages(library(TFBSTools))
suppressMessages(library(motifmatchr))
suppressMessages(library(regioneR))
suppressMessages(library(Matrix))
p <- arg_parser('script description')
p <- add_argument(p, "motif", help = "motif path for jaspar, must tranlate by yourself")
p <- add_argument(p, "bed", help = "bed")
if ( is.null(p$help) || length(alen) < 5) {
    print(p)
    quit(status=1)
}
args <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))


peaks <- toGRanges( args$bed )
motifs <- readJASPARMatrix( args$motif, matrixClass="PFM")
motifs_ix <- matchMotifs(motifs, peaks, genome = "BSgenome.Mmusculus.UCSC.mm10", out = "scores")
#pfm <- getMatrixByName(JASPAR2020, name = "MYC")
#MycMotifs <- matchMotifs( pfm, peaks, genome = "BSgenome.Mmusculus.UCSC.mm10", out = "positions")
#export.bed(MycMotifs[[1]], con = "MycMotifs.bed")
df <- motifCounts(motifs_ix)
writeMM(df, file= paste( basename(args$bed), basename(args$motif),'motifs', sep = '.') )

#write.table( df,  file = paste( basename(args$bed), basename(args$motif),'motifs', sep = '.'), quote = F, sep="\t", row.names = F, col.names = T )



























