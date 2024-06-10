#!/usr/bin/env sh
if [ $# -lt 1 ]; then
echo $0 'sample'
exit
fi

##  HiTC ICE Normalization  ## 输入文件 HTCrawMatrix 
## R里执行
#if (( 0 )); then
#library("HiTC")
#for(sm in c("Sperm")){
#    a<-paste("HTC.",sm,sep="")
#    for(res in c("40k")){
#        list<-sapply(list.files(paste("/home/chenxp/3T3/mouse/",sm,"/R1/raw_matrix/ICE/40k",sep=""),pattern=a,full.names=TRUE),import.my5C)
#        hic<- HTClist(list)
#        hic<- hic[isIntraChrom(hic)]
#        hiC_iced <- HTClist(lapply(hic,normICE,max_iter=1500))
#        hiC_iced<- HTClist(lapply(hiC_iced,forceSymmetric))
#        c<-seqlevels(hiC_iced)
#        for( m in c){n<-paste(m,m,sep="");export.my5C(hiC_iced[[n]],file=paste("/home/chenxp/3T3/mouse/",sm,"/R1/Norm_matrix/ICE/40k/R1.",m,".",res,".",sm,".ICE_matrix.rout",sep=""),genome="mm10")}
#    }
#}
#fi

#  生成 clean ICE matrix
## shell 执行
species=rheMac8
sample=$1
type=autosome
version=rheMac8
genome=rheMac8
norm=ICE
windows=2000000
work_dir=`pwd`
outdir=$work_dir/../../../domain
script_path=/home/ningch/soft/HiC-CXP-pipeline/domaincall_software/perl_scripts
fai=/home/ningch/data/genome/rheMac8/rheMac8.fai
echo $outdir
exit 
for bin in 40000
    do
    res=`expr $bin / 1000`k
    DI=DI_from_matrix.pl
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 
        do
        for EXP in  Sperm
            do
            matrix_path=$work_dir
            outdir=$work_dir/../../../domain
            routmat=$matrix_path/$sample.$chr.$res.$EXP.ICE_matrix.rout.mat
            norm_matrix=$matrix_path/$EXP.$sample.$chr.$res.$version.$norm\_matrix.txt
            clean_norm_matrix=$matrix_path/clean.$EXP.$sample.$chr.$res.$version.$norm\_matrix.txt
            perl $script_path/norm2DImatrix.pl $routmat $chr $fai > $norm_matrix
            less $norm_matrix | cut -f 4- - > $clean_norm_matrix
            rm $norm_matrix
            done
        done
    done
exit 
#  Normailize Sequence Depth
## R里执行
Norm=10000
for(sm in c("Sperm")){
dir=paste("/home/chenxp/3T3/mouse/",sm,"/R1",sep="")
for(res in c("40k")){
tmpdir<- paste(dir,"/Norm_matrix/ICE/40k/",sep="")
list<-sapply(list.files(path = tmpdir,pattern = paste("clean.",sm,sep=""),full.names = T),read.table,header=F)
list<-sapply(list,as.matrix)
names(list)<- paste("chr",c(10:13,1,14:19,2:9),sep="")  ###### 染色体顺序要注意
for(chr in names(list)){
a<-rowSums(list[[chr]])
depth<-a[a!=0][1]
NormDepth<- Norm*list[[chr]]/depth
write.table(NormDepth,file=paste("/home/chenxp/3T3/mouse/",sm,"/R1/Norm_matrix/ICE/40k/R1.clean.",chr,".",res,".",sm,".ICE.scaleTo1.1w.txt",sep=""),col.names = F,sep="\t",row.names = F,quote = F)
}
}
}

