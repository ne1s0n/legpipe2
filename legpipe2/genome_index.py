#index the speficied genome file with bowtie2, picard and "samtools faidx",
#plus generate the intervals lenghts file

import subprocess
import os
import common
from Bio import SeqIO
import gzip
import re

def validate(conf):
	'''validate incoming config parameters from .ini file and env variables'''
	#checking if files/paths exist
	if not os.path.exists(conf['genome_index']['reference_file']):
		msg = 'Reference file does not exist: ' + conf['genome_index']['reference_file']
		raise FileNotFoundError(msg)
	
	#checking if picard env variable is there
	if os.environ.get('PICARD') is None:
		msg = 'You need to set the environmental variable $PICARD to point to your picard.jar'
		raise EnvironmentError(msg)


def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	return(conf)

def genome_index(conf):

	#interface
	common.print_step_header('Indexing reference genome')
	
	#should we do something?
	RUN_THIS=conf['genome_index']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#----------- config
	REFERENCE_FILE=conf['genome_index']['reference_file']
	BOWTIE_INDEX=conf['genome_index']['bowtie_index']
	REGION_LENGTHS_FILE=conf['genome_index']['region_lengths_file']
		
	# ------------ bowtie
	#bowtie needs to run in the reference genome folder, so that
	#all the created files stay there
	cmd_str = 'cd ' + common.fn(os.path.dirname(REFERENCE_FILE)) + '; '
	cmd_str += 'bowtie2-build ' + common.fn(REFERENCE_FILE) + ' ' + common.fn(BOWTIE_INDEX)
	print(cmd_str)
	subprocess.run(cmd_str, shell=True)

	# ------------ samtools faidx
	cmd = ['samtools', 'faidx', REFERENCE_FILE]
	print(' '.join(cmd))
	subprocess.run(cmd, shell=False)
	
	# ------------ picard
	#substitute .fasta.gz, .fa.gz and all the upper/lower case variants with .dict
	picard_out = re.sub('(?i)\.fasta\.gz$', '.dict', REFERENCE_FILE)
	picard_out = re.sub('(?i)\.fa\.gz$', '.dict', picard_out)
	
	cmd = ['java', '-jar', os.environ.get('PICARD'), 'CreateSequenceDictionary']
	cmd += ['-R',  REFERENCE_FILE]
	cmd += ['-O',  picard_out]
	print(' '.join(cmd))
	subprocess.run(cmd, shell=False)

	# ------------ region lenghts
	print('Region lengths are stored in ' + REGION_LENGTHS_FILE)
	with open(REGION_LENGTHS_FILE, 'w') as fp_out:
		with gzip.open(REFERENCE_FILE, "rt") as fp_in:
			for seq_record in SeqIO.parse(fp_in, "fasta"):
				fp_out.write(seq_record.id + ':1-' + str(len(seq_record)) + '\n')	

