#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: human_globinrna_check.sh <read_1.fastq.gz> <number_of_reads> <output_file>"
    exit 1
fi

R1=$1

N=$(($2*4))
OUTFILE=$3

if [ ! -f $R1 ]; then
    echo "ERROR: Read 1 file not found!"
    exit 1
fi



TEMP=$HOME/.rrna_tmp
mkdir -p $TEMP

zcat $R1 | head -n $N | gzip > $TEMP/rrna_r1.fastq.gz
 

bwa mem -t 12 ~/lib/globin_fasta/Human_globin.fa $TEMP/rrna_r1.fastq.gz | samtools view -@ 12 -bS -F 2304 -o $TEMP/rrna.bam -

samtools flagstat $TEMP/rrna.bam > $TEMP/rrna_stat
RRNA_COUNTED=`head -n1 $TEMP/rrna_stat | cut -f1 -d\ `
RRNA_PERCENT=`head -n3 $TEMP/rrna_stat | tail -n1 | cut -f1 -d\: | cut -f2 -d\(`

echo -e "Read_1\tReads_checked\tHuman_globinRNA_percent"
echo -e "${R1}\t${2}\t${RRNA_PERCENT}"

echo -e "${R1}\t${2}\t${RRNA_PERCENT}" >> $OUTFILE

rm -rf $TEMP/*

echo "Deleting temporary files"

echo "Done!"
                                