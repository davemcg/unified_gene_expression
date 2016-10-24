#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import sys
import glob

parser = argparse.ArgumentParser(description=\
	'Wrapper script that downloads dbGaP protected sra, converts to fastq, and sends to salmon for quantification. MUST RUN IN /data/mcgaugheyd/dbGaP/11588.')
parser.add_argument('sample_accession',\
	help="SRA sample accession")
parser.add_argument('run_accession',\
	help="SRA run accession, comma separated if there are multiple per sample")
parser.add_argument('library',\
	choices=['paired','single'],
	help="Specify \"paired\" or \"single\" ended library")

args = parser.parse_args()
sample_accession = args.sample_accession
run_accession = args.run_accession
library = args.library

# let's keep sra and fastq files in /scratch
sra_path = '/data/mcgaugheyd/dbGaP/11588/sra'
# create folder (don't care if it already exists) for counts
salmon_main_dir = '/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/salmon_counts'
mkdir_salmon_call = 'mkdir ' + salmon_main_dir + '/' + sample_accession
subprocess.call(mkdir_salmon_call, shell=True)
# create main folder in /scratch
subprocess.call('mkdir /scratch/mcgaugheyd', shell=True)

# download each sra file 
try:
	for i in run_accession.split(','):
		prefetch_call = 'prefetch ' +  i 
		subprocess.check_call(prefetch_call, shell=True) 
except:
	print('prefetch fail for ' + i)
	sys.exit(1)

# salmon quantification time
if library == 'paired':
	salmon_call = 'sbatch --cpus-per-task 8 ~/git/unified_gene_expression/scripts/salmon_stream.sh ' + \
				  sra_path + '/' + run_accession + '.sra paired ' + salmon_main_dir + '/' + sample_accession + ' delete'
	subprocess.check_call(salmon_call, shell=True)
if library == 'single': 
	salmon_call = 'sbatch --cpus-per-task 8 ~/git/unified_gene_expression/scripts/salmon_stream.sh ' + \
				  sra_path + '/' + run_accession + '.sra single ' +  salmon_main_dir + '/' + sample_accession + ' delete'
	subprocess.check_call(salmon_call, shell=True)

