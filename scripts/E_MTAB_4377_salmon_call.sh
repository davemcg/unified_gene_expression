#!/bin/bash

cd /data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/E-MTAB-4377/

for i in *_1.fastq; do sbatch --partition=quick --cpus-per-task 16 --mem=16G ~/git/unified_gene_expression/scripts/run_salmon.sh $i ${i%_1.fastq}_2.fastq  paired ${i%_1.fastq}; done
