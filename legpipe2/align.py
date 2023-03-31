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
#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 


#----------- SUPPORT FUNCTIONS
def validate(conf):
	'''validate incoming config parameters from .ini file'''

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#these values should be boolean
	conf['align']['skip_previously_completed'] = raw_conf['align'].getboolean('skip_previously_completed') 

	#these values should be int
	conf['align']['cores'] = raw_conf['align'].getint('cores') 
	
	#if max_samples is zero it goes to +Infinity
	conf['align']['max_samples'] = raw_conf['align'].getint('max_samples')
	if conf['align']['max_samples'] == 0:
		conf['align']['max_samples'] = float('inf')

	return(conf)

def _do_align(infile_R1, outfolder, bowtie_index):
	'''this function is designed to be executed in parallel, once per 
	input fastq R1/R2 files'''

	#--------- filenames
	fn = create_filenames(infile_R1, outfolder)
	
	#--------- bowtie2 align
	cmd = ['bowtie2', '-x', 'bowtie_index']
	cmd += ['-1', fn['infile_R1']]
	cmd += ['-2', fn['infile_R2']]
	cmd += ['-S', fn['tmp_sam']]
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(log_bowtie2_align, "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools for sam -> bam conversion
	cmd = ['samtools', 'view']
	cmd += ['-bS', fn['tmp_sam']]
	cmd += ['>', fn['tmp_bam']]
	subprocess.run(cmd, shell=True)
	
	#--------- picard, AddOrReplaceReadGroups
	cmd = ['java', '-jar', '/home/ubuntu/software/picard.jar', 'AddOrReplaceReadGroups']
	cmd += ['-I',  fn['tmp_bam']]
	cmd += ['-O',  fn['tmp_bam_groups']]
	cmd += ['-LB', 'Whatever', '-PL', 'Illumina', '-PU', 'Whatever', '-SM', fn['core']]
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(log_picard_readGroups, "w") as fp:
		fp.write(res.stdout)

	#--------- picard, ValidateSamFile
	cmd = ['java', '-jar', '/home/ubuntu/software/picard.jar', 'ValidateSamFile']
	cmd += ['-INPUT', fn['tmp_bam_groups']]
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(log_picard_validation, "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools, sort
	cmd = ['samtools', 'sort', fn['tmp_bam_groups']]
	cmd += ['-o', fn['outfile']]
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(log_samtools_sort, "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools, index
	cmd = ['samtools', 'index', fn['outfile']]
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(log_samtools_index, "w") as fp:
		fp.write(res.stdout)
	
	#--------- cleanup of intermediate files
	subprocess.run(['rm', fn['tmp_sam']], shell=True)
	subprocess.run(['rm', fn['tmp_bam']], shell=True)
	subprocess.run(['rm', fn['tmp_bam_groups']], shell=True)

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
	res['log_bowtie2_align']     = res['tmp_sam'].replace('.sam', '.bowtie2_align.log')
	res['log_picard_validation'] = res['tmp_sam'].replace('.sam', '.picard_validation.log')
	res['log_picard_readGroups'] = res['tmp_sam'].replace('.sam', '.picard_readGroups.log')
	res['log_samtools_sort']     = res['tmp_sam'].replace('.sam', '.samtools_sort.log')
	res['log_samtools_index']    = res['tmp_sam'].replace('.sam', '.samtools_index.log')
	
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
	
	#room for output
	cmd_str = "mkdir -p " + OUTFOLDER
	subprocess.run(cmd_str, shell=True)
	
	return(None)


	#collecting all the arguments for the parallel execution in a pandas df
	args = None
	skipped = 0
	for infile_R1 in glob.glob(INFOLDER + '/*_R1.fastq.gz'):
		#should we skip this file?
		fn = create_filenames(infile_R1, OUTFOLDER)
		if os.path.isfile(res['outfile_index']) and SKIP_ALIGNED:
			skipped += 1
			print('Skipping previously aligned sample ' + fn['core'])
			continue
		
		#the arguments for the current file. Column order is important,
		#it should match the order for the parallel function, since 
		#the arguments are passed as positionals.
		args_now = pd.DataFrame({
			'infile_R1' : [infile_R1], 
			'outfolder' : [OUTFOLDER], 
			'bowtie_index' : [BOWTIE_INDEX]
		})

		#storing in a single df
		args = pd.concat([args, args_now])

		#safeguard
		if len(args) >= MAX_SAMPLES:
			break

	#executing in parallel using multiprocessing module
	cnt = 0
	with ThreadPool(CORES) as pool:
		for result in pool.starmap(align, args.itertuples(index = False)):
			#result variable contains the core of the processed sample, but
			#we don't want to flood the main screen (there's many log files)
			#so we just do nothing with it
			cnt += 1
	print('Total sample processed: ' + str(cnt))
	print(' - of which, skipped because previous runs: ' + str(skipped))

