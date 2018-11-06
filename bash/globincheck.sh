#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Run in folder with .fastq.gz files. Usage: human_globinrna_check.sh <project_id>"
    exit 1
fi

PROJ=$1

#R1=$1

N=$((1000000*4))
FNAME="_globin.txt"
OUTFILE=$1$FNAME

for R1 in *fastq.gz; do

echo $R1
TEMP=$HOME/.rrna_tmp
mkdir -p $TEMP
zcat $R1 | head -n $N | gzip > $TEMP/rrna_r1.fastq.gz
bwa mem -t 12 ~/lib/globin_fasta/Human_globin.fa $TEMP/rrna_r1.fastq.gz \
| samtools view -@ 12 -bS -F 2304 -o $TEMP/rrna.bam -

samtools flagstat $TEMP/rrna.bam > $TEMP/rrna_stat
RRNA_COUNTED=`head -n1 $TEMP/rrna_stat | cut -f1 -d\ `
RRNA_PERCENT=`head -n3 $TEMP/rrna_stat | tail -n1 | cut -f1 -d\: | cut -f2 -d\(`

echo -e "Read_1\tReads_checked\tHuman_globinRNA_percent"
echo -e "${R1}\t${2}\t${RRNA_PERCENT}"

echo -e "${R1}\t${2}\t${RRNA_PERCENT}" >> $OUTFILE

rm -rf $TEMP/*

echo "Deleting temporary files"

echo "Done!"
done