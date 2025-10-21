#!/bin/bash
# Example BWA mapping command
GENOME="./data/genome/refgenome.fasta"
READ1="./data/reads/sample01_1.fq.gz"
READ2="./data/reads/sample01_2.fq.gz"
SAMPLE="sample01"
BAM_OUT="./data/bam/sample01.bam"

bwa mem -R "@RG\tID=$SAMPLE\tSM=$SAMPLE\tPL=ILLUMINA" $GENOME $READ1 $READ2 |     samtools view -bS - |     samtools sort -o $BAM_OUT
