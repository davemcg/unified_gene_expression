#!/bin/bash
module load samtools

srr_id=$1
library=$2
output=$3

id=$srr_id

echo $id $2 $3
if [[ $srr_id =~ .*\.sra$ ]]; then id=${srr_id%.sra}; fi
echo $id $2 $3

if [[ "$library" == "paired" ]]; then
	mkdir /scratch/mcgaugheyd/$id
	sam-dump $srr_id | samtools view -bu - | \
		samtools bam2fq - -1 /scratch/mcgaugheyd/$id/${id}_1.fastq -2 /scratch/mcgaugheyd/$id/${id}_2.fastq
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
		--libType A --seqBias --gcBias \
		-1 <(cat /scratch/mcgaugheyd/$id/${id}_1.fastq) \
		-2 <(cat /scratch/mcgaugheyd/$id/${id}_2.fastq)  \
		-o $output \
		-p $SLURM_CPUS_PER_TASK
elif [[ "$library" == "single" ]]; then
	~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon quant --index /data/mcgaugheyd/genomes/GRCh38/salmon_0.7.2_gencode_v25/ \
        --libType A --seqBias --gcBias \
		-r <(sam-dump $srr_id | samtools view -bu - | samtools bam2fq -) \
		-o $output \
		-p $SLURM_CPUS_PER_TASK  
else
	echo $library "not 'paired' or 'single'"
fi
