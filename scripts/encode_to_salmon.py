#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import sys
import glob
import re

parser = argparse.ArgumentParser(description=\
	'Wrapper script that downloads fastqs from encode and sends to salmon for quantification.')
parser.add_argument('encode_sample_ID',\
	help="ENCODE ID")
parser.add_argument('encode_fastq_urls',\
	help="ENCODE fastq urls. Comma separated. If multiple sets, put in order, splitting each set by ||. E.g. forward_rep1,reverse_rep1||forward_rep2, reverse_rep2")
parser.add_argument('library',\
	choices=['paired','single'],
	help="Specify \"paired\" or \"single\" ended library")

args = parser.parse_args()
encode_sample_ID = args.encode_sample_ID
encode_fastq_urls = args.encode_fastq_urls
library = args.library

# create main folder in /scratch
subprocess.call('mkdir /scratch/mcgaugheyd', shell=True)
# let's download fastq files to /scratch
temp_folder = '/scratch/mcgaugheyd/' + encode_sample_ID
# create folder (don't care if it already exists) for counts
salmon_main_dir = '/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/salmon_counts_bootstrap50_txUsed'
mkdir_salmon_call = 'mkdir ' + salmon_main_dir + '/' + encode_sample_ID
subprocess.call(mkdir_salmon_call, shell=True)


# download each fastq file 
#try:
#	for i in re.split('__|,',encode_fastq_urls):
#		wget_call = 'wget --directory-prefix=' + temp_folder + ' ' + i 
#		subprocess.check_call(wget_call, shell=True) 
#except:
#	print('wget command failed: ' + wget_call)
#	sys.exit(1)

# salmon quantification time
if library == 'paired':
	forward_reads = [x.split(',')[0].split('/')[-1] for x in encode_fastq_urls.split('__')]
	reverse_reads = [x.split(',')[1].split('/')[-1] for x in encode_fastq_urls.split('__')]
  
	forward_reads = ' '.join([temp_folder + '/' + x for x in forward_reads])
	reverse_reads = ' '.join([temp_folder + '/' + x for x in reverse_reads])
	salmon_call = 'sbatch --time=6:00:00 --cpus-per-task 16 ~/git/unified_gene_expression/scripts/run_salmon2.sh ' + \
				  '\\"' + forward_reads + '\\"\ ' + '\\"' + reverse_reads + '\\"' + \
           ' paired ' + salmon_main_dir + '/' + encode_sample_ID
	subprocess.check_call(salmon_call, shell=True)
if library == 'single':
	reads = [x.split(',')[0].split('/')[-1] for x in encode_fastq_urls.split('__')]
	reads = ' '.join([temp_folder + '/' + x for x in reads])
	salmon_call = 'sbatch --time=6:00:00 --cpus-per-task 16 ~/git/unified_gene_expression/scripts/run_salmon2.sh ' + \
				  '\"' + reads + '\"\ ' + '\"' + \
           ' single ' + salmon_main_dir + '/' + encode_sample_ID
	subprocess.check_call(salmon_call, shell=True)
