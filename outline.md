# Abstract
GTEx provides gene/transcript level expression datasets across dozens of people over dozens of human tissues. However, the eye was not used. Over a dozen disparate studies have looked at the transcriptome in different portions of the human eye, more commonly in the RPE, retina, and cornea. Here we collate and re-process all of the available human eye-tissue RNA-seq data. We also process GTEx and ENCODE tissue and cell line data, respectively, in an identical manner. This allows for a comparison of gene expression signatures of eye tissues relative to GTEx and ENCODE data sets. This large dataset allows us to more confidently define retina, RPE, and cornea-specific gene expression patterns, look at correlated gene networks between the functionally disparate but physically adjacent RPE and retina, and identify eQTL and ASE. Data and plots are shared in an interactive web application (link here).

# Intro
- RNA-seq overview
- GTEx and ENCODE
  * overview
  * value
- What does eye have?
  * no inclusion in GTEx and ENCODE
  * disparate kind-of-small n studies
- RNA-seq processing is expensive computationally
  * alignment and gene-level counting for each sample takes many (6+) hours
- New algorithms (kallisto/sleuth) allow far faster computation
  * directly get gene counts in less than 2 hours on a 8-core node
- Collating publically available RNA-seq sets gives over 200 eye-specific human tissues
  * allows for higher confidence DE (retina <-> RPE <-> cornea)
  * ID eye-specific genes (not or very lowly expressed in other tissues)
  * ID genes NOT expressed in eye tissues
- Integration of GTEx and ENCODE data allows for high-value utilization of public data
  * if any GTEx tissues have grossly similar expression to eye tissues, then you can leverage eQTL, ASE from that tissue
  * if ENCODE matches by RNA-seq, then you can leverage ChIP-seq data
  
# Results
- How samples were selected
  * Queried SRA database for *retina*, *RPE*, *cornea*, *eye* with *RNA-seq*
- Workflow
  * Stage 1: ID samples, pull metadata, pull SRA file, convert to fastq, quantitate with Salmon
  * Stage 2: tximport to get gene-level TPM, remove genes with low counts across all samples, remove Samples with low median gene counts > 1, t-sne cluster to ID outliers, qsmooth to quantile normalize TPM
- Accounting
  * Eye, GTEx, ENCODE 
  * table: tissue, study, number of samples (pre and post filter), reads processed, base pairs sequenced, alignment percentage
- Clustering
  * recapitulate organ grouping in t-sne parsing
  * see fetal (either tissue or cell line) RPE and retina cluster away from adult RPE and retina
  * Retina is a unique tissue
  * RPE/cornea tends to stay together
    * closest GTEx tissue is transformed fibroblasts to RPE/cornea
  * ENCODE cell lines
    * cluster far apart from everything else
- Basic gene list parsing
  * what is unique (only expressed in eye) to retina, rpe, cornea?
  * what is not expressed in eye tissues?
- DE with DESeq2
  * RPE vs Retina
  * 'fetal' vs adult (both RPE and Retina)
- DTE with Sleuth?
- Transcriptome
  * cufflinks
  * or try to get pre-published transcriptomes for retina/rpe and add to gencode v25 transcript fasta file then re-run salmon?
- Network analysis (WGCNA, GSEA)?
- Pseudo ASE/eQTL
  * Use RNA-seq data to call SNPs, then can do ASE/eQTL analysis
  * https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0152-4
  * Would have to call RNA-seq with STAR.....
- Loop back to known eye biology
  * AMD loci?
  * GTEx fibroblast ASE/eQTL have any eye biology links?
  * Get Rob/Brian/OGVFB feedback/idea!!!
