#!/usr/bin/env Rscript
suppressMessages(library(argparser))
p <- arg_parser('script description')
p <- add_argument(p, "-gtf", short = '-g', help="gtf file", default= "/home/soft/data/genome/rheMac8/Macaca_mulatta.Mmul_8.0.1.97.gtf")
p <- add_argument(p, "-genome", help="gtf file", default= "/home/soft/data/genome/rheMac8/rheMac8.fa")
p <- add_argument(p, "--sam", short = '-s', nargs = Inf, help="sam file")
p <- add_argument(p, "--out", short = '-o', help="out file")
if ( is.null(p$help) || length( commandArgs() ) < 5 ) {
        print(p)
        quit(status=1)
}

args <- parse_args(p)
library(Rsubread)
feature = featureCounts(files=args$sam,annot.ext=args$gtf,isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="transcript_id", isPairedEnd=TRUE,minFragLengt=20,maxFragLength=100000,nthreads=6,strandSpecific=2,reportReads='BAM', juncCounts=TRUE, genome=args$genome)
#feature = featureCounts(files=args$sam,annot.ext=args$gtf, isPairedEnd=TRUE,minFragLengt=20,maxFragLength=100000,nthreads=6,strandSpecific=0,reportReads='BAM')
fc = feature
write.table(
	x=data.frame(fc$annotation[,c("GeneID","Length")],
    fc$counts,
    stringsAsFactors=FALSE),
  	file=args$out,
  	quote=FALSE,
  	sep="\t",
  	row.names=FALSE)

#sink(file=args$out)
#feature$count
#sink()




