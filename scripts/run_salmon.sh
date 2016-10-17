#!/bin/bash

fastq=$1
paired=$2
output=$3

echo $1 $2 $3

if [["$paired" == "yes"]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCm38/salmon_0.7.2_gencode_v10 \
		--libType A --seqBias --gcBias \
		-1 <($fastq | zcat ) \
		-2 <(${fastq%_1.fastq.gz}_2.fastq.gz | zcat )  \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
elif [["$paired" == "no"]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCm38/salmon_0.7.2_gencode_v10 \
        --libType A --seqBias --gcBias \
		-r <($fastq | zcat)
		-o $output
		-p $SLURM_CPUS_PER_TASK
else
	echo $paired "not 'yes' or 'no'"
fi
