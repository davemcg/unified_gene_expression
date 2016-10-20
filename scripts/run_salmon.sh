#!/bin/bash

fastq_1=$1
fastq_2=$2
library=$3
output=$4

echo $1 $2 $3 $4

if [[ "$library" == "paired" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
		--libType A --seqBias --gcBias \
		-1 <(zcat -f $fastq_1) \
		-2 <(zcat -f $fastq_2)  \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
elif [[ "$library" == "single" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
        --libType A --seqBias --gcBias \
		-r <(zcat -f $fastq_1) \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
else
	echo $library "not 'paired' or 'single'"
fi
