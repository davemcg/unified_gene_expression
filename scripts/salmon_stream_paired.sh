#!/bin/bash

# example input: full ftp url
# will stream to salmon, not saving fastq
# ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030887/ERR030887_1.fastq.gz
fastq1_url=$1
output=$2 
~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCm38/salmon_0.7.2_gencode_v10 \
	--libType A --seqBias --gcBias \
	-1 <(curl -s $fastq_url | zcat ) \
	-2 <(curl -s ${fastq_url%_1.fastq.gz}_2.fastq.gz | zcat )  \
	-o $output \
	-p $SLURM_CPUS_PER_TASK
