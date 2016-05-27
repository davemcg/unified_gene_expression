**Workflow**

1. Genomes, transcripts, and annotation:
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz
 - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.metadata.HGNC.gz
 
2. Going cutting edge and using kallisto for reference-free transcrdipt alignment. Why? Speed. I'd also really like to try Salmon, which looks fabulous.
	https://github.com/davemcg/biowulf2-bin/blob/master/kallisto_quant.sh
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
	Different script for single end. Have to explicitly give insert size (usually 200) and SD (usually 50).
	https://github.com/davemcg/biowulf2-bin/blob/master/kallisto_quant_single.sh
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
 - https://liorpachter.wordpress.com/2014/04/30/estimating-number-of-transcripts-from-rna-seq-measurements-and-why-i-believe-in-paywall
 - https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

	I ended up going with the Mike Love tximport tool, which deals with both core issues: aggregating transcripts to gene and performing a TPM-like calculation to normalize counts/expression by transcript length and library size (lengthScaledTPM). My implementation is done in R with calculate_lengthScaledTPM.R (https://github.com/davemcg/unified_gene_expression/blob/master/calculate_lengthScaledTPM.R).

4. Next is massaging the SRA metadata to better name the files, since they are currently labeled with "SRRlotsofdigits." For now I'll just re-name the files to reflect the cell type and source (something like "RPE_fetal" or "RPE_iPSC"). It's going to be ugly. I was hoping to use some kind of R package (SRAdb?) to pull the info I needed automatically, but I can't get anything to work. Since the plan is to hand-curate this data source, for now it is reasonable to just hand-pull the info. The two sources I'm using to get meta-data are https://www.ebi.ac.uk and http://www.ncbi.nlm.nih.gov/Traces/study/ . The former contains the fastq ftp links with some basic info on the experiment. The latter has more granular info on the biological info on the experiment.
 - The code that adds the meta-data and gives example plots:
 	- https://github.com/davemcg/unified_gene_expression/blob/master/plot_by_gene.R
5. Then I'll need to figure out some way of displaying/sharing the data with OGVFB. ~~This will likely be very tricky~~. The current plan is to run a Shiny web-app on Cyclops that will allow users to type a gene name in and get the lengthScaledTPM scores across the tissues (medians for replicates?) along with some basic stats (rank/decile expression scoring?). Put ui.R and server.R on Cyclops along with data files metaData_with_lengthScaledTPM.Rdata and all_human_genes.Rdata (both derived from plot_by_gene.R)
6. http://cyclops.nei.nih.gov:3838/unified_gene_expression/
	- only works in NIH/NEI

