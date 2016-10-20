#!/bin/bash
module load samtools

bam=$1

samtools sort -n $bam | samtools fastq - -1 ${bam%.bam}_1.fastq -2 ${bam%.bam}_2.fastq
