#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import sys
import glob

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

# let's keep sra and fastq files in /scratch
sra_path = '/scratch/mcgaugheyd/' + sample_accession
# create folder (don't care if it already exists) for counts
salmon_main_dir = '/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/salmon_counts_bootstrap50_txUsed'
mkdir_salmon_call = 'mkdir ' + salmon_main_dir + '/' + sample_accession
subprocess.call(mkdir_salmon_call, shell=True)
# create main folder in /scratch
subprocess.call('mkdir /scratch/mcgaugheyd', shell=True)

# create folder for sample_accession
try:
	mkdir_call = 'mkdir ' + sra_path 
	subprocess.check_call(mkdir_call, shell=True)
except:
	print('mkdir ' + sra_path + ' fail!')
	sys.exit(1)

# download each sra file to this folder
try:
	for i in ftp_links.split(','):
		wget_call = 'wget --tries=50 -P ' + sra_path + ' ' + i 
		subprocess.check_call(wget_call, shell=True) 
except:
	print('wget fail for ' + ftp_links)
	sys.exit(1)

# convert to fastq
try:	
	if library == 'paired':
		for i in ftp_links.split(','):
			run_accession = i.split('/')[-1]
			fastq_dump_call = 'fastq-dump --clip -I --gzip --split-files --outdir ' + \
				sra_path + ' ' + \
				sra_path + '/' + run_accession
			print(fastq_dump_call)
			subprocess.check_call(fastq_dump_call, shell=True)  
	if library == 'single':
		for i in ftp_links.split(','):
			run_accession = i.split('/')[-1]
			fastq_dump_call = 'fastq-dump --clip --gzip --outdir ' + \
				sra_path + ' ' + \
				sra_path+ '/' + run_accession
			subprocess.check_call(fastq_dump_call, shell=True)
except:
	print('fastq-dump fail!')

# remove sra files
try:
	rm_call = 'rm ' + sra_path + '/*sra'
	subprocess.check_call(rm_call, shell=True)
except:
	print('rm fail')

# salmon quantification time
if library == 'paired':
	f_reads = glob.glob(sra_path + '/*_1.fastq.gz')
	r_reads = glob.glob(sra_path + '/*_2.fastq.gz')
	f_reads = str(' '.join(f_reads))
	r_reads = str(' '.join(r_reads))
	salmon_call = 'sbatch  --time=6:00:00  --cpus-per-task 16 ~/git/unified_gene_expression/scripts/run_salmon2.sh' + \
					r' \"' + f_reads + r'\" \"' + r_reads + r'\" paired ' + salmon_main_dir + '/' + sample_accession
	subprocess.check_call(salmon_call, shell=True)
if library == 'single': 
	reads = glob.glob(sra_path + '/*fastq.gz')
	reads = str(' '.join(reads))
	salmon_call = 'sbatch --time=6:00:00 --cpus-per-task 16 ~/git/unified_gene_expression/scripts/run_salmon2.sh' + \
					r' \"' + reads + r'\" \"' + 'INTENTIONAL_BLANK' + r'\" single ' + salmon_main_dir + '/' + sample_accession
	subprocess.check_call(salmon_call, shell=True)


