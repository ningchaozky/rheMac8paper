#!/usr/bin/env sh
if [ $# -lt 1 ]; then
    echo $0 'final.bam' 'peirod' 'rep' 'res'
    exit
fi

######################################################################################
##  生成cis/trans BedPE文件  ##
bam=$1
peirod=$2
rep=$3
EXP=PFC
genome=rheMac8
enzyme=Mbol
resolution=40000
res=40k
ref=/home/ningch/data/genome/rheMac8/rheMac8.fa
outdir=`pwd`/$rep
perl_path=/home/ningch/soft/HiC-CXP-pipeline/pl
genome_file=/home/ningch/data/genome/rheMac8/rheMac8.genome
fragments_dir=/home/ningch/data/genome/rheMac8/Mbol/rheMac8_Mbol/rheMac8_Mbol_fragments/
if [ ! -e $outdir/raw_matrix/ICE/$res ]; then
echo $outdir/raw_matrix/ICE/$res
mkdir -p $outdir/cis
mkdir -p $outdir/trans
mkdir -p $outdir/raw_matrix
mkdir -p $outdir/raw_matrix/ICE
mkdir -p $outdir/raw_matrix/ICE/$res
mkdir -p $outdir/Norm_matrix
mkdir -p $outdir/Norm_matrix/ICE
mkdir -p $outdir/Norm_matrix/ICE/$res
mkdir -p $outdir/domain
mkdir -p $outdir/domain/ICE
mkdir -p $outdir/domain/ICE/ScaleTo1
fi


echo bam2bedpe
cis_output=$outdir/cis/cis.$peirod\_hicup_dedup.sort.bedpe
if [ ! -e $cis_output ]; then
    echo generate file $cis_output
    samtools view $bam | grep "CT:Z:FAR" |samtools view -bS -T $ref - | bedtools bamtobed -bedpe -mate1 -i - |awk '$8>=10{print $0}' -|awk '$1==$4{print $0}' - > $cis_output
fi


if (( 0 )); then
trans_output=$outdir/trans/trans.$peirod\_hicup_dedup.bedpe
if [ ! -e $trans_output ]; then
    echo gnerate file: $trans_output
    samtools view $bam | grep "CT:Z:TRANS"|samtools view -bS -T $ref - | bedtools bamtobed -bedpe -mate1 -i - |awk '$8>=10{print $0}' -|awk '$1!=$4{print $0}' - > $trans_output
fi
fi


echo sort bedpe by chr
if [ ! -e $outdir/cis/cis.chr1.$peirod\_hicup_dedup.sort.bedpe ]; then
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX 
do
    sort_bedpe=$outdir/cis/cis.$chr.$peirod\_hicup_dedup.sort.bedpe
    echo sort $cis_output to $sort_bedpe
    awk -v chr="$chr" '$1==chr{print $0}' $cis_output | sort -k1,1 -k2,2n - > $sort_bedpe
done
fi

echo build FRMTdist
if [ ! -e $outdir/cis/cis.chr1.$peirod.hicup.FRMTdist.proper.read.txt ]; then
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX 
do
    sort_bedpe=$outdir/cis/cis.$chr.$peirod\_hicup_dedup.sort.bedpe
    reBed=$fragments_dir/$genome.$chr.$enzyme.REfragment.bed
    end5=$outdir/cis/cis.$chr.$peirod\_hicup_dedup.5end.sort.bedpe
    assignment=$outdir/cis/cis.$chr.$peirod.hicup.fragment.assignment.txt
    read1_assignment=$outdir/cis/cis.$chr.$peirod.hicup.read1.Fassignment.txt
    read2_assignment=$outdir/cis/cis.$chr.$peirod.hicup.read2.Fassignment.txt
    FRMTdist_proper_read=$outdir/cis/cis.$chr.$peirod.hicup.FRMTdist.proper.read.txt
    echo generate $end5
    awk 'BEGIN{OFS="\t"}{if($9=="+"){R1e=$2+1}else{R1e=$3};R1s=R1e-1;if($10=="+"){R2e=$5+1}else{R2e=$6};R2s=R2e-1;print $1"\t"R1s"\t"R1e"\t"$4"\t"R2s"\t"R2e"\t"$7"\t"$8"\t"$9"\t"$10}' $sort_bedpe > $end5
    echo generate $assignment
    bedtools pairtobed -f 1 -s -type both -a $end5 -b $reBed > $assignment
    echo generate $read1_assignment
    awk 'NR%2==1{print $0}' $assignment > $read1_assignment
    echo generate $read2_assignment
    awk 'NR%2==0{print $0}' $assignment > $read2_assignment
    echo generate $FRMTdist_proper_read
    paste $read1_assignment $read2_assignment | awk '$14!=$31{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$27"\t"$28"\t"$29"\t"$14"_"$30}' - | awk '! ( $2 == 0 || $5 == 0) { print }' > $FRMTdist_proper_read
done
#rm  $outdir/cis/*assig*.txt $outdir/cis/*5end* $outdir/cis/*chr*_hicup_dedup.sort.bedpe
fi
#######################################################################################
##  生成 ICE Raw Matrix  ##
echo generate raw matrix
if [ ! -e $outdir/raw_matrix/ICE/$res/HTC.$rep.chr1\_$peirod\_$res\_$genome\_matrix.txt ]; then
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX
do 
    gbin=$fragments_dir/$genome.$chr.bin.txt
    reBed=$fragments_dir/$genome.$chr.$enzyme.REfragment.bed
    gdist=$fragments_dir/$genome.$chr.bin.distribution.txt
    proper=$outdir/cis/cis.$chr.$peirod.hicup.FRMTdist.proper.read.txt
    chr_bin=$outdir/$chr.txt
    reads_count=$outdir/raw_matrix/ICE/$res/clean.$chr\_$peirod\_$res\_$genome\_matrix.txt
    ice_raw_matrix=$outdir/raw_matrix/ICE/$res/HTC.$rep.$chr\_$peirod\_$res\_$genome\_matrix.txt
    echo generate $gbin
    bedtools makewindows -g $genome_file -w $resolution | awk -v chr="$chr" '$1==chr{print $0}' - | awk '{print $0"\tbin"NR}' - > $gbin  #### 生成 固定size bin ###
    
    echo generate $gdist
    awk '$6=="+"{print $0}' $reBed |awk '{mid=int(($2+$3)/2);mids=mid-1;print $1"\t"mids"\t"mid"\t"$4"\t"$5"\t"$6}' - | bedtools intersect -f 1 -wa -wb -a - -b $gbin |grep -v "-1" - > $gdist
    
    echo generate $chr_bin
    awk -F "\t" 'NR==FNR{a[$4]=$10;next;}{split($17,b,"_");print a[b[1]]"\t"a[b[2]]"\t"$0}' $gdist $proper | awk '{gsub("bin","");print $0}' - | sort -k1,1n -k2,2n - | bedtools groupby -g 1,2 -c 9 -o count_distinct | awk '{print "bin"$1"\tbin"$2"\t"$3}' - > $chr_bin  ### 将reads分到指定bin里
    num=$(wc -l $gbin)
    
    echo generate $reads_count
    generate_plain.pl $chr_bin $num | awk '{gsub("bin","");print $0}' - | sort -k1,1n -k2,2n - | bin2matrix.pl - >  $reads_count  ### reads count计数 ###
    
    echo generate $ice_raw_matrix
    awk -v resolution=$resolution '{header=(NR-1)*resolution;print header"\t"$0}' $reads_count | HTCmatrix.pl - $chr $resolution $genome $genome_file > $ice_raw_matrix  ### 生成HiTC兼容的ICE raw natrix矩阵  ###

#rm $hind3/$species\_$enzyme\_fragments/$genome.$chr.bin.txt $hind3/$species\_$enzyme\_fragments/$genome.$chr.bin.distribution.txt $outdir/$chr.txt
done
fi






