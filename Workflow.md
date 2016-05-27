**Workflow**

1. Genomes, transcripts, and annotation:
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.metadata.HGNC.gz
 
2. Going cutting edge and using kallisto for alignment. Why? Speed.
```
[mcgaugheyd@biowulf public_RNA-seq]$ cat kallisto_quant.sh 
#!/bin/bash

module load kallisto/0.42.4

kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
			   -o ${1%_*}_kallisto \
               -b 100 \
               -t $SLURM_CPUS_PER_TASK \
			   $1 $2
```
Different script for single end. Have to explicitly give insert size (usually 200) and SD (usually 50)
```
[mcgaugheyd@biowulf public_RNA-seq]$ cat kallisto_quant_single.sh 
#!/bin/bash

module load kallisto/0.42.4

fastq1=$1
frag_length=$2
frag_sd=$3

kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
			   -o ${1%_*}_kallisto \
               -b 100 \
               -t $SLURM_CPUS_PER_TASK \
			   --single -l $frag_length -s $frag_sd  \
			   $fastq1
```
 
 
 
 
