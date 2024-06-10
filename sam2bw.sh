#!/usr/bin/env bash
help() {
    echo $0 sams
}
if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
        help
        exit 1
fi


if [[ $# -gt 2 ]]; then 
    bam=merge.bam
    if [[ ! -e merge.bam ]]; then 
        show samtools merge $bam $@
    fi
else 
    bam=$1
fi

sbam=`change $bam sorted.bam`
if [ ! -e $sbam ]; then
    show "samtools sort $bam > $sbam "
fi

if [ ! -e $sbam.bai ]; then
    show "samtools index $sbam "
fi

sbw=`change $bam sorted.bw`
if [ ! -e $sbw ]; then
    show "bamCoverage -b $sbam -bs 1 --maxFragmentLength 1000000 -o $sbw -p 6"
fi


markdup=`change $sbam marked_duplicates.bam`
if [ ! -e $markdup ]; then
    show "sambamba markdup -r -t 4 --overflow-list-size 6000000 $sbam $markdup"
fi

if [ ! -e $markdup.bai ]; then
    show "samtools index $markdup "
fi
mbw=`change $markdup bw`
#show "java -jar /home/soft/soft/picard.jar  MarkDuplicates I=$bam O=$bam.redup.bam M=marked_dup_metrics.txt"
if [ ! -e $mbw ]; then
    show "bamCoverage -b $markdup -bs 1 --maxFragmentLength 1000000 -o $mbw -p 6"
fi

























