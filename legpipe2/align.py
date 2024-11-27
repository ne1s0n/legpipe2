#align reads to the indexed genome. Tools used: bowtie2, samtools, picard
#Notes:
# - from INFOLDER we take all forward and reverse reads (*_R1.fastq.gz and *_R2.fastq.gz)
# - the genome is expected to be indexed as per index step

#for each sample there's a list of commands:
	#bowtie2 -x ref_98 -1 sample1.1.fq -2 sample1.2.fq -S sample1.sam
	#samtools view -bS sample1.sam > sample1.bam
	#java -jar picard/2.4.1/ValidateSamFile.jar INPUT= sample1.bam
	#java -jar picard/2.4.1/AddOrReplaceReadGroups.jar I=sample1.bam O=sample1.Gr.bam LB=Whatever PL=Illumina PU=Whatever SM=sample1
	#samtools/1.3.1/samtools sort sample1.Gr.bam sample1.Gr.sorted.bam
	#samtools/1.3.1/samtools index sample1.Gr.sorted.bam

import subprocess
import glob
import os
import pandas as pd
import common
import sys

#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 

#----------- SUPPORT FUNCTIONS
def validate(conf):
	'''validate incoming config parameters from .ini file, plus environmental variables'''
	if os.environ.get('PICARD') is None:
		msg = 'You need to set the environmental variable $PICARD to point to your picard.jar'
		raise EnvironmentError(msg)
	if not os.path.exists(conf['align']['reference_file']):
		msg = 'Reference file does not exist: ' + conf['align']['reference_file']
		raise FileNotFoundError(msg)


def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#these values should be boolean
	conf['align']['skip_previously_completed'] = raw_conf['align'].getboolean('skip_previously_completed') 
	conf['align']['paired'] = raw_conf['align'].getboolean('paired') 

	#these values should be int
	conf['align']['cores'] = raw_conf['align'].getint('cores') 
	
	#if max_samples is zero it goes to +Infinity
	conf['align']['max_samples'] = raw_conf['align'].getint('max_samples')
	if conf['align']['max_samples'] == 0:
		conf['align']['max_samples'] = float('inf')

	return(conf)

def _do_align(infile_R1, infile_R2, outfolder, bowtie_index, paired):
	'''this function is designed to be executed in parallel, once per 
	input fastq R1/R2 files'''
	
	#because of how this function is invoked, all arguments are converted to string
	#let's compensate it
	paired = paired == 'True'
	
	#--------- filenames
	fn = _create_filenames(infile_R1, outfolder)
	
	#--------- bowtie2 align
	cmd = ['bowtie2', '-x', bowtie_index]
	if paired:
		cmd += ['-1', fn['infile_R1']]
		cmd += ['-2', fn['infile_R2']]
	else:
		cmd += ['-U', fn['infile_R1']]
		
	cmd += ['-S', fn['tmp_sam']]
	with open(fn['logfile'], "w") as fp:
		fp.write('\n\n--------- bowtie2 align\n')
		fp.write('Special bowtie call: -f for single end\n')
		fp.write('paired: "' + str(paired) + '"\n')
		fp.write('of type: "' + str(type(paired)) + '"\n')
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#--------- samtools for sam -> bam conversion
	cmd = ['samtools', 'view']
	cmd += ['-bS', fn['tmp_sam']]
	with open(fn['logfile'], "a") as fp:
		fp.write('\n\n--------- samtools for sam -> bam conversion\n')
		fp.write(' '.join(cmd) + '\n')
	with open(fn['tmp_bam'], "w") as fp:
		subprocess.run(cmd, shell=False, stdout=fp)
		
	#--------- picard, AddOrReplaceReadGroups
	cmd = ['java', '-jar', os.environ.get('PICARD'), 'AddOrReplaceReadGroups']
	cmd += ['-I',  fn['tmp_bam']]
	cmd += ['-O',  fn['tmp_bam_groups']]
	cmd += ['-LB', 'Whatever', '-PL', 'Illumina', '-PU', 'Whatever', '-SM', fn['core'], '-ID', fn['core']]
	with open(fn['logfile'], "a") as fp:
		fp.write('\n\n--------- picard, AddOrReplaceReadGroups\n')
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)

	#--------- picard, ValidateSamFile
	cmd = ['java', '-jar', os.environ.get('PICARD'), 'ValidateSamFile']
	cmd += ['-INPUT', fn['tmp_bam_groups']]
	with open(fn['logfile'], "a") as fp:
		fp.write('\n\n--------- picard, ValidateSamFile\n')
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#--------- samtools, sort
	cmd = ['samtools', 'sort', fn['tmp_bam_groups']]
	cmd += ['-o', fn['outfile']]
	with open(fn['logfile'], "a") as fp:
		fp.write('\n\n--------- samtools, sort\n')
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#--------- samtools, index
	cmd = ['samtools', 'index', fn['outfile']]
	with open(fn['logfile'], "a") as fp:
		fp.write('\n\n--------- samtools, index\n')
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#--------- cleanup of intermediate files
	subprocess.run(['rm', fn['tmp_sam']], shell=False)
	subprocess.run(['rm', fn['tmp_bam']], shell=False)
	subprocess.run(['rm', fn['tmp_bam_groups']], shell=False)

	#--------- done
	return(fn['core'])

def _create_filenames(infile_R1, outfolder):
	'''returns a dictionary with all the filenames derived from the input R1
	and the outfolder'''

	res = {'infile_R1' : infile_R1}
	
	#core sample name, without path and file extension
	res['core'] = os.path.basename(infile_R1).replace('_R1.fastq.gz', '')
	
	#the other input file
	res['infile_R2'] = infile_R1.replace('_R1', '_R2')

	#intermediate files, deletable at the end
	res['tmp_sam'] = outfolder + '/' + res['core'] + '.sam'
	res['tmp_bam'] = res['tmp_sam'].replace('.sam', '.bam')
	res['tmp_bam_groups'] = res['tmp_sam'].replace('.sam', '.gr.bam')
	
	#log files
	res['logfile']     = res['tmp_sam'].replace('.sam', '.align.log')
	
	#output file
	res['outfile'] = res['tmp_sam'].replace('.sam', '.gr.sorted.bam')
	
	#index for output file
	res['outfile_index'] = res['outfile'] + '.bai'
	
	return(res)

#----------- ALIGN (public) FUNCTION
def align(conf):
	#interface
	common.print_step_header('align')
	
	#should we do something?
	RUN_THIS=conf['align']['run_this']

	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	BOWTIE_INDEX=conf['align']['bowtie_index']
	INFOLDER=conf['align']['infolder']
	OUTFOLDER=conf['align']['outfolder']
	CORES=conf['align']['cores']
	MAX_SAMPLES=conf['align']['max_samples']
	SKIP_PREVIOUSLY_COMPLETED=conf['align']['skip_previously_completed']
	PAIRED=conf['align']['paired']
	
	#room for output
	cmd_str = "mkdir -p " + common.fn(OUTFOLDER)
	subprocess.run(cmd_str, shell=True)
	
	#collecting all the arguments for the parallel execution in a pandas df
	args = None
	skipped = 0
	for infile_R1 in glob.glob(INFOLDER + '/*_R1.fastq.gz'):
		#should we skip this file?
		fn = _create_filenames(infile_R1, OUTFOLDER)
		
		print('Aligning ' + fn['core'])
		if os.path.isfile(fn['outfile_index']) and SKIP_PREVIOUSLY_COMPLETED:
			skipped += 1
			print(' - skipping, previously aligned')
			continue
		
		#the arguments for the current file. Column order is important,
		#it should match the order for the parallel function, since 
		#the arguments are passed as positionals
		args_now = pd.DataFrame({
			'infile_R1' : [fn['infile_R1']], 
			'infile_R2' : [fn['infile_R2']], 
			'outfolder' : [OUTFOLDER], 
			'bowtie_index' : [BOWTIE_INDEX],
			'paired' : [PAIRED]
		})

		#storing in a single df
		args = pd.concat([args, args_now])

		#safeguard
		if len(args) >= MAX_SAMPLES:
			break

	#do we have something to execute?
	cnt = 0
	if args is None:
		print('All samples skipped, no operation required')
	else:
		#executing in parallel using multiprocessing module
		with ThreadPool(CORES) as pool:
			for result in pool.starmap(_do_align, args.itertuples(index = False)):
				#result variable contains the core of the processed sample, but
				#we don't want to flood the main screen (there's many log files)
				#so we just do nothing with it
				cnt += 1
				
	print('Total sample processed: ' + str(cnt))
	print('Skipped because previous runs: ' + str(skipped))

