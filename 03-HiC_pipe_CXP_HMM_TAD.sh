#!/usr/bin/env sh
if [ $# -lt 1 ]; then
echo $0 ''
exit
fi
##  DI ##
species=mouse
sample=R1
type=autosome
version=mm10
genome=mm10
norm=ICE
windows=2000000

for bin in 40000
do

res=`expr $bin / 1000`k

script_path=/home/chenxp/Hi-C/bin/domaincall_software/perl_scripts
DI=DI_from_matrix.pl
fai=/home/chenxp/Reference/$species/$version/$species.$version.chrall.fa.fai

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
do

for EXP in Sperm
do

matrix_path=/home/chenxp/3T3/$species/$EXP/$sample/Norm_matrix/ICE/${res}
mkdir /home/chenxp/3T3/$species/$EXP/$sample/domain/ICE/ScaleTo1/$res
mkdir /home/chenxp/3T3/$species/$EXP/$sample/domain/ICE/ScaleTo1/$res/unadj
outdir=/home/chenxp/3T3/$species/$EXP/$sample/domain/ICE/ScaleTo1/$res/unadj


less $matrix_path/$sample.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt|awk -v bin=$bin '{header=(NR-1)*bin;print header"\t"$0}' -|perl $script_path/normdepth2DImatrix.pl - $chr $bin $version /home/chenxp/Reference/$species/$genome/$genome.genome > $matrix_path/DI.$sample.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt
perl $script_path/$DI $matrix_path/DI.$sample.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt $bin $windows $fai |grep -v "inf" >$outdir/$EXP.$sample.$chr.$res.$version.$norm\_DI.scaleTo1.1w.txt 
less $outdir/$EXP.$sample.$chr.$res.$version.$norm\_DI.scaleTo1.1w.txt |awk '{print "chr"$0}' ->$outdir/$EXP.$sample.$chr.$res.$version.$norm\_DI.scaleTo1.1w.bedGraph

rm $matrix_path/DI.$sample.clean.$chr.$res.$EXP.ICE.scaleTo1.1w.txt
cat $outdir/$EXP.$sample.$chr.$res.$version.$norm\_DI.scaleTo1.1w.txt >>$outdir/$EXP.$sample.chrall.DI.scaleTo1.1w.txt
cat $outdir/$EXP.$sample.$chr.$res.$version.$norm\_DI.scaleTo1.1w.bedGraph >>$outdir/$EXP.$sample.chrall.DI.scaleTo1.1w.bedGraph
done

done
done

for EXP in Sperm
do
Rep=R1
cd /home/chenxp/3T3/mouse/$EXP/$Rep/domain/ICE/ScaleTo1/40k/unadj
sort -k1,1 -k2,2n $EXP.$Rep.chrall.DI.scaleTo1.1w.bedGraph -o $EXP.$Rep.chrall.DI.scaleTo1.1w.bedGraph
bedGraphToBigWig $EXP.$Rep.chrall.DI.scaleTo1.1w.bedGraph ~/Reference/mouse/mouse/mm10.genome $EXP.$Rep.chrall.DI.scaleTo1.1w.bw
mv $EXP.$Rep.chrall.DI.scaleTo1.1w.bw /home/chenxp/3T3/mouse/Epi/mouseHiC/DI
done

##
matlab   # 执行HMM_call.m
##
species=mouse
sample=R1
map=bwa
version=mm10
norm=ICE
res=40k

perl_path=~/Hi-C/bin/domaincall_software/perl_scripts
fai=~/Reference/$species/$version/$species.$version.chrall.fa.fai


min=2
prob=0.99
binsize=40000

for EXP in Sperm
do
dir=/home/chenxp/3T3/$species/$EXP/$sample/domain/ICE/ScaleTo1/$res/unadj
mkdir $dir/HMM$min
mkdir $dir/HMM$min/chr
rm $dir/HMM$min/chr/$EXP.$sample.chrall.$res.$version.$norm\_domain.bed
perl $perl_path/1file_ends_cleaner.pl $dir/$EXP.$sample.$res.$version.${norm}_hmmout.scaleTo1.1w.txt $dir/$EXP.$sample.chrall.DI.scaleTo1.1w.txt |perl $perl_path/2converter_7col.pl >$dir/$EXP.$sample.$res.$version.$norm\_7colfile.txt

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
do
awk -v chr=$chr '$1==chr{print $0}' $dir/$EXP.$sample.$res.$version.$norm\_7colfile.txt >$dir/HMM$min/chr/$EXP.$sample.$chr.$res.$version.$norm\_7colfile.txt

perl $perl_path/4hmm_probablity_correcter.pl $dir/HMM$min/chr/$EXP.$sample.$chr.$res.$version.$norm\_7colfile.txt $min $prob $binsize|perl $perl_path/5hmm-state_caller.pl $fai $chr|awk 'NR>=2{print}' -|perl $perl_path/6hmm-state_domains.pl -|awk 'NF==3{print}' ->$dir/HMM$min/chr/$EXP.$sample.$chr.$res.$version.$norm\_domain.bed

cat $dir/HMM$min/chr/$EXP.$sample.$chr.$res.$version.$norm\_domain.bed >>$dir/HMM$min/chr/$EXP.$sample.chrall.$res.$version.$norm\_domain.bed
done
wc -l $dir/HMM$min/chr/$EXP.$sample.chrall.$res.$version.$norm\_domain.bed >>$dir/$EXP.$sample.TAD.summary.txt
done

