#!/bin/bash
module load samtools/1.3.1
module load sratoolkit/2.8.0

srr_id=$1
library=$2
output=$3
file_fate=$4

if [[ "$library" == "paired" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25_commonTx/ \
		--libType A --seqBias --gcBias --numBootstraps 100 \
		-1 <(fastq-dump -I --split-files --stdout $srr_id | grep '^@.*\..*\.1 ' -A 3 --no-group-separator) \
		-2 <(fastq-dump -I --split-files --stdout $srr_id | grep '^@.*\..*\.2 ' -A 3 --no-group-separator)  \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
elif [[ "$library" == "single" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25_commonTx/ \
        --libType A --seqBias --gcBias --numBootstraps 100 \
		-r <(fastq-dump --stdout $srr_id) \
		-o $output \
		-p $SLURM_CPUS_PER_TASK  
else
	echo $library "not 'paired' or 'single'"
fi

if [[ "$file_fate" == "delete" ]]; then
	rm $srr_id
fi
