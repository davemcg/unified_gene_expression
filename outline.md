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
  
