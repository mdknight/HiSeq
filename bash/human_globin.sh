#!/bin/bash

#updated script to calculate % of reads in .fastq align to globin reference for SR (or R1 for PE) data
#will perform calculations on all .fastq files in current directory and save to one .csv
if [ $# -ne 1 ]; then
    echo "Run in folder with .fastq.gz files. Usage: human_globin.sh <project_id>"
    exit 1
fi

#used to identify run or project
PROJ=$1
#number of reads to use from fastq
READS=1000000
#number of lines to read from fastq
N=$((READS*4))
REPDIR="/data/fastq/Reports/"
FNAME="_globin.txt"
OUTFILE=$1$FNAME

#output header to file, saved in reports folder
#append to prevent overwriting
echo -e "Read_1,Reads_checked,Human_globin_percent" >> $REPDIR$OUTFILE

#will find all directories with run/proj ID from input
#exclude report folders
PROJDIR=$(find /data/fastq -name ${PROJ}* -type d -not -path "*/Reports/*")
echo "Run/Project ID was found in the following directories:"
echo ${PROJDIR}

#go through each directory to find .fastq files to process
for CURPROJ in $PROJDIR;do

#search for .fastq.gz files for read 1 and save to list
echo "Searching for .fastq files in the following directory:"
echo $CURPROJ

PROJFASTQ=$(find ${CURPROJ} -name *R1_001.fastq.gz)

#for each fastq, calculate globin and output to file
for R1 in $PROJFASTQ; do
echo $R1
TEMP=$HOME/.rrna_tmp
mkdir -p $TEMP
#output first N lines to temp file 
zcat $R1 | head -n $N | gzip > $TEMP/rrna_r1.fastq.gz
#align temp file to globin fasta save as temp .bam
bwa mem -t 12 ~/lib/globin_fasta/Human_globin.fa $TEMP/rrna_r1.fastq.gz \
| samtools view -@ 12 -bS -F 2304 -o $TEMP/rrna.bam -
#calculate globin content 
samtools flagstat $TEMP/rrna.bam > $TEMP/rrna_stat
RRNA_COUNTED=`head -n1 $TEMP/rrna_stat | cut -f1 -d\ `
RRNA_PERCENT=`head -n3 $TEMP/rrna_stat | tail -n1 | cut -f1 -d\: | cut -f2 -d\(`
#output results to terminal
echo -e "Read_1\tReads_checked\tHuman_globinRNA_percent"
echo -e "${R1}\t${READS}\t${RRNA_PERCENT}"
#save results to file
echo -e "${R1},${READS},${RRNA_PERCENT}" >> $REPDIR$OUTFILE
#remove temporary files
rm -rf $TEMP/*
echo "Deleting temporary files"
echo "Done!"
done

done