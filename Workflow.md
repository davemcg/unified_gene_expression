**Workflow**

1. Genome and annotation:
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz
 
2. Create STAR indices, with different indices for different read lengths 

```
[mcgaugheyd@biowulf GRCh38]$ cat create_STAR_index.sh 
#!/bin/bash

module load STAR/2.5.1b 


index_size=$1
# Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.

STAR \
	--runThreadN $SLURM_CPUS_PER_TASK \
	--runMode genomeGenerate \
	--genomeDir /data/mcgaugheyd/genomes/GRCh38/STAR_indices/$1 \
	--genomeFastaFiles /data/mcgaugheyd/genomes/GRCh38/GRCh38.primary_assembly.genome.fa \
	--sjdbGTFfile /data/mcgaugheyd/genomes/GRCh38/gencode.v24.primary_assembly.annotation.gtf \
	--sjdbOverhang $1
```

```
 [mcgaugheyd@biowulf GRCh38]$ cat STAR_index_creation.sh 
 sbatch --cpus-per-task 10 --mem=30G create_STAR_index.sh 35
 sbatch --cpus-per-task 10 --mem=30G create_STAR_index.sh 49
 sbatch --cpus-per-task 10 --mem=30G create_STAR_index.sh 100
 sbatch --cpus-per-task 10 --mem=30G create_STAR_index.sh 129
```
3. Align with STAR 2.5.1b, with appropriate sjdbOverhang for different read lengths
