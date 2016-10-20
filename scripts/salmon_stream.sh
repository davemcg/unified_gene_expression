#!/bin/bash
module load samtools/1.3.1
module load sratoolkit/2.8.0

srr_id=$1
library=$2
output=$3


if [[ "$library" == "paired" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
		--libType A --seqBias --gcBias \
		-1 <(sam-dump $srr_id | samtools view -bu - | samtools fastq - | grep '^@.*/1$' -A 3 --no-group-separator) \
		-2 <(sam-dump $srr_id | samtools view -bu - | samtools fastq - | grep '^@.*/2$' -A 3 --no-group-separator)  \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
elif [[ "$library" == "single" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
        --libType A --seqBias --gcBias \
		-r <(sam-dump $srr_id | samtools view -bu - | samtools fastq -) \
		-o $output \
		-p $SLURM_CPUS_PER_TASK  
else
	echo $library "not 'paired' or 'single'"
fi
