#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: human_rrna_check.sh <read_1.fastq.gz> <read_2.fastq.gz> <number_of_read_pairs> <output_file>"
    exit 1
fi

R1=$1
R2=$2
N=$(($3*4))
OUTFILE=$4

if [ ! -f $R1 ]; then
    echo "ERROR: Read 1 file not found!"
    exit 1
fi

if [ ! -f $R2 ]; then
    echo "ERROR: Read 2 file not found!"
    exit 1
fi

TEMP=$HOME/.rrna_tmp
mkdir -p $TEMP

zcat $R1 | head -n $N | gzip > $TEMP/rrna_r1.fastq.gz
zcat $R2 | head -n $N | gzip > $TEMP/rrna_r2.fastq.gz

bwa mem -t 12 /opt/references/bwa_rRNA/human_rRNA.fa $TEMP/rrna_r1.fastq.gz $TEMP/rrna_r2.fastq.gz | samtools view -@ 12 -bS -F 2304 -o $TEMP/rrna.bam -

samtools flagstat $TEMP/rrna.bam > $TEMP/rrna_stat
RRNA_COUNTED=`head -n1 $TEMP/rrna_stat | cut -f1 -d\ `
RRNA_PERCENT=`head -n3 $TEMP/rrna_stat | tail -n1 | cut -f1 -d\: | cut -f2 -d\(`

echo -e "Read_1\tRead_2\tReads_checked\tHuman_rRNA_percent" 
echo -e "${R1}\t${R2}\t${3}\t${RRNA_PERCENT}"

echo -e "${R1},${R2},${3},${RRNA_PERCENT}" >> $OUTFILE

rm -rf $TEMP/*

echo "Deleting temporary files"

echo "Done!"