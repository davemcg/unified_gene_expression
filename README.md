# unified_gene_expression
Meta-analyze and collate genome-wide gene expression from <b>normal human tissues and popular cell lines</b> (GTEx and ENCODE, respectively) and eye tissues.

These large repositories, unfortunately, largely lack eye tissues or cell lines (RPE/choroid, retina, cornea). Several studies have done RNA-seq in (relatively to GTEx) small numbers of samples in RPE/choroid (both derived from embryonic stem cells and straight tissue),  retina (tissue, some have done macula vs peripheral), or cornea. However they are all disjointed at the bioinformatic level. There would be tremendous value in unifying analysis across all of these tissues. 

GTEx is a tremendous resource, sequencing dozens to hundreds of different people across many human organs. 

Integration with this GTEx would be valuable for at least two reasons: 

1. Are RPE or retina global gene expression approximately similar to any of the GTEx tissues? 
- then it would be possible to use GTEx's allele-specific expression set to look whether common variation influences genes known to be important in eye function (e.g. ABCA4)
2. Which genes are differentially expressed in eye relative to other human tissues?

Unfortunately GTEx uses a very boutique (some would say controversial) method of RNA-seq quantification. 
- https://liorpachter.wordpress.com/2013/10/21/gtex/

ENCODE is also tremendous, integrating RNA-seq/methyl-seq/ChIP-seq across many cell lines. While they have a more standard pipeline, it would be expensive to replicate with the GTEx data. 

GTEx's software is used by no one but themselves. To make this a more useful resouce, I think it's important to re-process <b> everything </b> in a similar manner. 



This will be computationally expensive, especially with a 'classic' workflow (aligment to reference genome). This project has two things in its favor:
1. hpc.nih.gov
2. salmon/kallisto

biowulf2 is a huge computer cluster. Coupled with the very fast and accurate salmon/kallisto for gene/transcript quantification, this is tractable issue. 

## Querying the SRA for relevant eye RNA-seq datasets
https://github.com/davemcg/unified_gene_expression/blob/master/scripts/sraDB_search_select.R leverages the R package sraDB to download the entire sequence read archive metadata. Then the tidyverse is used to ID relevant run accessions. As of 2016-10-12, I've identified 179 normal samples, covering RPE/choroid, retina, and cornea tissue. There are another 50 retina samples from another study (E-MTAB-4377) that are not in the SRA. So, 229 tissues total. 

This is fairly streamlined, but needs to be manually re-run periodically to identify new studies. https://github.com/davemcg/unified_gene_expression/blob/master/data/eye_studies_hand_checked_2016-10-12.Rdata contains the studies that I've hand checked on the search parameters in https://github.com/davemcg/unified_gene_expression/blob/master/scripts/sraDB_search_select.R. In the future I can simply rerun the search and check studies that aren't in the Rdata file.

# Bioinformatic workflow
## SRA to transcript-level counts
Very roughly, SRA acquired via wget or sratoolkit's sam-dump, converted to fastq, gene/transcript levels quantified with Salmon against Gencode v25 protein coding transcripts. This initial effort was broken into three categories:

1. E-MTAB-4377, which is only available with bam files from EBI/ENA
2. Directly from SRA to counts, not retaining .sra or derived .fastq files
3. dbGap access for some samples - download and hold .sra, then process

Scripts that largely cover these cases:

- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/create_salmon_index.sh (1,2,3)
- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/E_MTAB_4377_bam_handling.sh (1)
- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/E_MTAB_4377_salmon_call.sh (1)
- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/sra_to_salmon.py (2)
- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/run_salmon.sh (1,2)
- https://github.com/davemcg/unified_gene_expression/blob/master/scripts/salmon_stream.sh (1,3)

## 

