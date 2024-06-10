#!/usr/bin/env Rscript
#first
alen <- commandArgs()
library(prodlim)
library(data.table)
library(e1071)
suppressMessages( library(pheatmap) )
suppressMessages( library(dplyr) )
suppressMessages( library(viridis) )
suppressMessages( library(Rmpfr) )
suppressMessages( library(circlize) )
suppressMessages( library(seriation) )
suppressMessages( library(stringr) )
suppressMessages( library(TCseq) )
suppressMessages(library(argparser))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(amap)) # for Dist
suppressMessages(library(dendsort)) #sort
suppressMessages(library(cluster)) #sort
suppressMessages(library(matrixStats)) #sort
suppressMessages(library(Rclusterpp)) #sort
p <- arg_parser('script description')
p <- add_argument(p, "mat", help="number of significant digits to print")
p <- add_argument(p, "k", help="kmean cluster number")
p <- add_argument(p, "-c", nargs = '*', help="culster manual define include all for sort the position like: 6,2,1,4,3,5")
p <- add_argument(p, "-ymax", nargs = '?', help = "heatmap max set" )
p <- add_argument(p, "-ymin", nargs = '?', help = "heatmap min set" )
p <- add_argument(p, "-z", help = "zscore", flag = TRUE)
p <- add_argument(p, "-l", help = "+1log2", flag = TRUE)
if ( is.null(p$help) || length(alen) < 5) {
    print(p)
    quit(status=1)
}
args <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
   )
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))
#mat = t(mat)
#cluster columns
#group = kmeans(t(mat), centers = args$k)$cluster
#if args$c check the unique
if ( is.character(args$c)  ){
    manual_order = as.integer( unlist(strsplit(args$c,','))  )
    if ( length(unique(manual_order)) != length(manual_order)){
        #print ( sort.int(manual_order)
        print ( unique(manual_order) )
        stop('#args$c must unique')
    }
}

mat <- read.table( args$mat, sep = '\t', header = 1, check.names = 0, row.names = 1 )
#replace 0 by 0.01

#kmean_row = kmeans(mat, centers = args$k, iter.max=50)
kmean_row = cmeans( mat, centers = as.numeric( args$k ), iter.max=500 )
group = kmean_row$cluster
raw <- mat
print ('#Warning: Not normal anything you should normal by yourself')
#RNA-seq +1log(2)
#mat[ mat < 0 ] <- 0
if ( args$l ){
    mat <- log( mat + 1, base = 2)
}
#zscore
if ( args$z ){
    mat <- (mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
}
#center
#mat <- t(scale( t(mat) ))
#dat2 <- dat %>% mutate_at(c("y", "z"), ~(scale(.) %>% as.vector))
#mat <- mat %>% mutate_all(~(scale(.) %>% as.vector))
#hclust <- Rclusterpp.hclust
#k split
#Heatmap(mat, name = "mat", cluster_columns = cluster_within_group(mat, group), row_km = args$k)
#dendrogram split
#Heatmap(mat, name = "mat", cluster_columns = cluster_within_group(mat, group), row_split = args$k)
#Heatmap(mat, name = "mat", cluster_rows = cluster_within_group(mat, group))
#Heatmap(mat, name = "mat", cluster_rows = cluster_within_group(t(mat), group))
# "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall"
#use fast_hclust method
pdf( paste( args$mat, "kMeansGroup.cluster.pdf", sep = '.'))
ht_opt$fast_hclust = TRUE
#fh = function(x) fastcluster::hclust(Dist(x, method = "pearson"), method="ward.D2")
#row_dend = dendsort(fh(mat))

tca <- timeclust( as.matrix( mat ), algo = "cm", k = as.numeric(args$k), standardize = FALSE)
group = tca@cluster
row_dend = cluster_within_group( t(mat), group )
outmat <- mat[order.dendrogram(row_dend),]
#ht_row_cmeans_order <- unique(group[rownames(outmat)])
#df <- cbind( outmat, group = paste0('c',group[rownames(outmat)]))
df <- cbind( outmat, group = group[rownames(outmat)])
df = data.frame( 'gene_id' = rownames(df), df )
write.table( df, file = paste( args$mat, "ClustOrder.Normal", sep = '.'), quote = F, sep="\t",row.names = F, col.names = T)
df <- cbind( raw[rownames(outmat),], group = paste0('c',group[rownames(outmat)]))
df <- data.frame('gene_id' = rownames(df), df)
write.table( df, file = paste( args$mat, "ClustOrder.raw", sep = '.'),quote = F,sep="\t",row.names = F, col.names = T)
#Heatmap( mat, name = "mat", cluster_rows = row_dend, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D2", cluster_columns = FALSE, row_dend_reorder = TRUE)
#Heatmap( mat, name = "mat", cluster_rows = row_dend, cluster_columns = FALSE, row_dend_reorder = TRUE)
#split = data.frame(cutree(hclust(dist(mat)), k = 2), rep(c("A", "B"), 9))
#ht <- Heatmap( mat, name = "mat", cluster_rows = cluster_within_group(t(mat), group), cluster_columns = FALSE, row_dend_reorder = TRUE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), clustering_distance_rows = function(x, y) 1 - cor(x, y), row_split = 16)
#if define cluster will use the order show heatmap
CreateAdjacencyMatrix <- function(x) {
  s <- gsub("\\.", "", x)
  m <- matrix(0, 10, 10)
  for (i in 1:(nchar(s)-1)) {
    m[as.numeric(substr(s, i, i))+1,
      as.numeric(substr(s, i+1, i+1))+1] <-
      m[as.numeric(substr(s, i, i))+1,
        as.numeric(substr(s, i+1, i+1))+1]+1
  }
  rownames(m) = 0:9
  colnames(m) = 0:9
  m
}
m1 = CreateAdjacencyMatrix(formatMpfr(Const("pi",2000)))
col = colorRamp2(quantile(m1, seq(0, 1, by = 0.1)), magma(11))

#col_fun = colorRamp2(c(0.2, 10), c("white",'red'))
#print ( str(col_fun) )
#ht <- Heatmap( mat, name = "mat", col = col, cluster_rows = cluster_within_group(t(mat), group), cluster_columns = FALSE, row_dend_reorder = TRUE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), clustering_distance_rows = function(x, y) 1 - cor(x, y), row_split = as.numeric(args$k), show_row_names = F, row_title = 1:as.numeric(args$k), column_title_rot = 0 )
row_title = 1:as.numeric(args$k)
#row_title = letters[1:as.numeric(args$k)]
row_title = unique(unlist(df$group))
print ( row_title )
ht <- Heatmap( mat, name = "mat", cluster_rows = cluster_within_group(t(mat), group), cluster_columns = FALSE, row_dend_reorder = TRUE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), clustering_distance_rows = function(x, y) 1 - cor(x, y), row_split = as.numeric(args$k), show_row_names = F, row_title = row_title, column_title_rot = 0 )
draw(ht)
#test <- cluster_within_group( t(mat), group )
#print ( group )
#mat <- as.matrix(mat)
#TSP GW BEA_TSP
print ( str(list_seriation_methods("matrix")) )
#o = seriate(max(mat) - mat, method = "BEA_TSP")
#ht <- Heatmap( mat, name = "mat", row_order = get_order(o, 1), cluster_columns = FALSE, row_dend_reorder = FALSE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), show_row_names = F, column_title_rot = 0 )
#draw(ht)
#ht <- Heatmap( mat, name = "mat", cluster_rows = cluster_within_group(t(mat), group), cluster_columns = FALSE, row_dend_reorder = FALSE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), clustering_distance_rows = function(x, y) 1 - cor(x, y), row_split = manual_order, show_row_names = F, row_title = paste0('c',ht_row_cmeans_order), column_title_rot = 0)
#ht = ht + rowAnnotation(rn = anno_text(rownames(mat)))


#write cluster within group resulst
if ( is.character(args$c) ) {
    #manual_order = as.integer( unlist(strsplit(args$c,',')) )
    print ( paste(c('Length is: ',length(manual_order), 'Order is: ', manual_order, 'Uniqe order is', unique(manual_order)), collapse = " "))
    #group = group[order(match(group, manual_order))]
    #mat = mat[names(group),]
    need_raw_order <- row_order(ht)
    names(need_raw_order) = row_title
    out = c()
    #print ( row_order(ht) )
    for (i in manual_order){
        i <- paste0('c',i)
        #print ( mat[row_order(ht)[[i]], ] )
        clu <- t(t(row.names(mat[need_raw_order[[i]],])))
        #clu <- t(t(row.names(mat[row_order(ht)[[i]],])))
        #clu <- t(t(colnames(mat[, column_order(ht)[[i]]])))
        #out <- cbind(c1, paste("cluster", i, sep=""))
        out <- c(out, clu)
    }
	df = mat[out,]
    out_df = data.frame("symbol"=rownames(df), df)
    colnames(out_df) <- str_replace( colnames(out_df), "X", "")
    write.table( out_df,  file = paste( args$mat, "ClustOrder.raw.withInGroup", sep = '.'), quote = F,sep="\t",row.names = F, col.names = T )
    ht <- Heatmap( mat[out,], name = "mat", cluster_rows = FALSE, cluster_columns = FALSE, row_dend_reorder = FALSE, heatmap_width = unit(16, "cm"), heatmap_height = unit(18, "cm"), show_row_names = F, column_title_rot = 0 )
    draw(ht)
}

#capture.output( summary(row_order(ht)), file = "MyNewFile.txt" )
#htorder <- data.frame(row_order(ht))
#write.table( htorder, file = paste( args$mat, "ClustOrder2", sep = '.'),quote = F,sep="\t",row.names = T,col.names = T)
pdf(paste(args$mat, 'timeclust','pdf', sep="."), width=16, height=16)
#tca@membership <- kmean_row$membership
#tca@cluster <- group
#tca@data <- obj@matrix
p <- timeclustplot(tca, categories = "timepoint", value = "expression", cols = 3, cl.color = "gray50", membership.color = rainbow(30, s = 3/4, v = 1, start = 1/6), title.size=8, axis.title.size=8, axis.text.size=8, axis.line.size = 1, legend.text.size=8, legend.title.size=8)

dev.off()

























