

#!/bin/bash

# take folder containing SRA reads, convert to fastq, send to salmon for quantification

module load sratoolkit/2.8.0


sample_accession=$1
sra_ftp_links=$2 #comma separated
paired_or_single=$3

# workflow
# download sra (or sras)
# use fastq-dump (--split for paired) to convert to fastq
# feed fastqs into salmon to generate counts

wget $sra_ftp_links && \
if paired=='yes':
	wget

#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description=\
	'Wrapper script that downloads sra, converts to fastq, and sends to salmon for quantification')
parser.add_argument('sample_accession',\
	help="SRA sample accession")
parser.add_argument('sra_ftp_links',\
	help="SRA ftp links, comma separated if there are multiple per sample")
parser.add_argument('library',\
	choices=['paired','single'],
	help="Specify \"paired\" or \"single\" ended library")

args = parser.parse_args()
sample_accession = args.sample_accession
ftp_links = args.sra_ftp_links
library = args.library

# create folder for sample_accession
try:
	mkdir_call = 'mkdir ' + sample_accession
	subprocess.check_call(mkdir_call, shell=True)
except:
	print('mkdir ' + sample_accession + ' fail!')
	sys.exit(1)

# download each sra file to this folder
try:
	for i in ftp_links.split(','):
		wget_call = 'wget --tries=50 -P ' + sample_accession + ' ' + i 
		subprocess.check_call(wget_call, shell=True) 
except:
	print('wget fail for ' + ftp_links)
	sys.exit(1)

# convert to fastq
try:	
	if library.lower == 'paired':
		for i in ftp_links.split(','):
			run_accession = i.split('/')[-1]
			fastq_dump_call = 'fastq-dump --gzip --split-files --outdir ' + 
				sample_accession + 
				sample_accession + '/' + run_accession
			subprocess.check_call(fastq_dump_call)  
		if library.lower == 'single':
			run_accession = i.split('/')[-1]
			fastq_dump_call = 'fastq-dump --gzip --outdir ' + 
				sample_accession + 
				sample_accession + '/' + run_accession
			subprocess.check_call(fastq_dump_call)  
except:
	print('fastq-dump fail!')

# salmon quantification time
 
