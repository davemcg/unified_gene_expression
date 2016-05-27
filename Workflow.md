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
3. Next is calculating some kind of RPKM/FPKM/TPM score that will semi-accurately reflect gene expression across a wide variety of library sizes. Another significant issue is reducing multiple transcripts into each gene. After lots of Googling, I've found the following useful readings:
 - http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
 - http://f1000research.com/articles/4-1521/v1
 - https://www.biostars.org/p/143458/
 - https://benchtobioinformatics.wordpress.com/2015/07/10/using-kallisto-for-gene-expression-analysis-of-published-rnaseq-data/
 - https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
I ended up going with the Mike Love tximport tool, which deals with both core issues: aggregating transcripts to gene and performing a TPM-like calculation to normalize counts/expression by transcript length and library size (lengthScaledTPM). My implementation is done in R with calculate_lengthScaledTPM.R (https://github.com/davemcg/unified_gene_expression/blob/master/calculate_lengthScaledTPM.R).
 
