#!/usr/bin/env sh
if [ $# -lt 1 ]; then
echo $0 ''
exit
fi

## Dekker Lab ## 安装详见 https://github.com/dekkerlab/cworld-dekker
## 利用dekkerlab包进行画图 heatmap.pl命令
## Matrix
species=mouse
resolution=40000
res=40k
genome=mm10
norm=ICE
perl_path=/home/chenxp/Hi-C/bin/domaincall_software/perl_scripts
dekkerlab=/home/chenxp/software/giorgetti-nature-2016-master
fai=/home/chenxp/Reference/$species/$genome/$species.$genome.chrall.fa.fai
for EXP in Sperm
do
for rep in R1
do
outdir=/home/chenxp/3T3/$species/$EXP/$rep
mkdir $outdir/Norm_matrix/ICE/${res}/dekkerlab
matrix_path=$outdir/Norm_matrix/ICE/${res}
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 
do
less $outdir/Norm_matrix/ICE/${res}/$rep.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt|awk -v resolution=$resolution '{header=(NR-1)*resolution;print header"\t"$0}' -|perl $perl_path/addHeader.pl - $chr $resolution $genome /home/chenxp/Reference/$species/$genome/$genome.genome > $outdir/Norm_matrix/ICE/${res}/dekkerlab/dekker.$rep.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt
gzip $outdir/Norm_matrix/ICE/${res}/dekkerlab/dekker.$rep.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt
done
done
done
## Plot
species=mouse
resolution=40000
genome=mm10
res=40k
R=R1
perl_path=/home/chenxp/Hi-C/bin/domaincall_software/perl_scripts
st=70842690
end=82779146
for EXP in Sperm
do
cd /home/chenxp/3T3/$species/$EXP/R1/Norm_matrix/ICE/$res/dekkerlab
for chr in chr9 
do
perl /home/chenxp/software/giorgetti-nature-2016-master/scripts/perl/subsetMatrix.pl -i dekker.R1.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt.gz --minDist 40000 --maxDist 4520000 --z ${chr}:${st}-${end}
perl /home/chenxp/software/giorgetti-nature-2016-master/scripts/perl/heatmap.pl -i dekker.R1.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt---${chr}-${st}-${end}.subset.matrix.gz  --dt  --iq 0  --pc white,red --start 0.5 --end 350 --mc grey --bg ## 可以改变颜色和color范围 详见heatmap.pl命令 
mv dekker.R1.clean.*.scaleTo1.1w.txt*.png /home/chenxp/Hi-C/figure/mouseEmbryo
rm dekker.R1.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt---${chr}*.subset.matrix.gz
done
done
