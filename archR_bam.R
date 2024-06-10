#!/usr/bin/env Rscript
#first
alen <- commandArgs()
suppressMessages(library(argparser))
p <- arg_parser('bam')
p <- add_argument(p, "bam", help="bam input")
p <- add_argument(p, "name", help="output dir name")
if ( is.null(p$help) || length(alen) < 5) {
    print(p)
    quit(status=1)
}
args <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))

suppressMessages(library(ArchR))
addArchRGenome('mm10')
inputFiles = args$bam
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = args$bam,
  minTSS = 0,
  minFrags = 0,
  maxFrags = 1e+15,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE, force=T
)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = args$name,
  minTSS = 6,
  minFrags = 0,
  maxFrags = 1e+05,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE, force=T
)
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = args$name,
  copyArrows = FALSE )
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

























