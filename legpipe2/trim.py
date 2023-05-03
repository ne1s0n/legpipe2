#trim reads to remove restriction sites and low quality bases using ‘fastp’ 
#(https://github.com/OpenGene/fastp). Put output in a new directory

#----------- IMPORT
import common
import subprocess
import glob
import os
import pandas as pd
import sys
import pickle
import shlex
#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 

#----------- SUPPORT FUNCTIONS
def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#the trim command, customized by user
	conf['trim']['cmd'] = shlex.split(conf['trim']['cmd'])
	
	#these values should be boolean
	conf['trim']['dry_run'] = raw_conf['trim'].getboolean('dry_run') 
	conf['trim']['skip_previously_completed'] = raw_conf['trim'].getboolean('skip_previously_completed') 

	#these values should be int
	conf['trim']['cores'] = raw_conf['trim'].getint('cores') 

	#if max_samples is zero it goes to +Infinity
	conf['trim']['max_samples'] = raw_conf['trim'].getint('max_samples')
	if conf['trim']['max_samples'] == 0:
		conf['trim']['max_samples'] = float('inf')
	
	return(conf)
	
def validate(conf):
	'''validate incoming config parameters from .ini file'''

def _create_filenames(infile_R1, outfolder):
	'''returns a dictionary with all the filenames derived from the input R1
	and the outfolder'''
	
	res = {'infile_R1' : infile_R1}
	
	#core sample name, without path and file extension
	res['core'] = os.path.basename(infile_R1).replace('_R1.fastq.gz', '')
	
	#the other input file
	res['infile_R2'] = infile_R1.replace('_R1', '_R2')

	#output files
	res['outfile_R1']   = outfolder + '/' + res['core'] + '_R1.fastq.gz'
	res['outfile_R2']   = outfolder + '/' + res['core'] + '_R2.fastq.gz'
	res['outfile_json'] = outfolder + '/' + res['core'] + '.json'
	res['outfile_html'] = outfolder + '/' + res['core'] + '.html'
	res['log']          = outfolder + '/' + res['core'] + '.fastp_trimming.log'
	
	return(res)

def _do_trim(infile_R1, outfolder, trim_cmd):
	'''executes the trimming for one R1/R2 pair'''
	#--------- filenames
	fn = _create_filenames(infile_R1, outfolder)
	
	#--------- fastp
	cmd = pickle.loads(trim_cmd)
	cmd += ['--in1' , fn['infile_R1'],  '--in2'  , fn['infile_R2']]
	cmd += ['--out1', fn['outfile_R1'], '--out2' , fn['outfile_R2']]
	cmd += ['-j' , fn['outfile_json']]
	cmd += ['-h' , fn['outfile_html']]
	
	print(cmd)
	
	res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log'], "w") as fp:
		fp.write(res.stdout)
	
	#--------- done
	return(fn['core'])

#----------- MAIN FUNCTION
def trim(conf):
	#interface
	common.print_step_header('trim')
	
	#should we do something?
	RUN_THIS=conf['trim']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#----------- config
	INFOLDER=conf['trim']['infolder']
	OUTFOLDER=conf['trim']['outfolder']
	#number of parallel threads
	CORES=conf['trim']['cores']
	#if True commands are just printed but not executed
	DRY_RUN=conf['trim']['dry_run']
	#maximum number of samples to process, for test purposes
	MAX_SAMPLES=conf['trim']['max_samples']
	#if TRUE samples already processed will be skipped
	SKIP_PREVIOUSLY_COMPLETED=conf['trim']['skip_previously_completed']
	#the actual trim command
	TRIM_CMD=conf['trim']['cmd']
	#room for output
	cmd_str = "mkdir -p " + OUTFOLDER
	subprocess.run(cmd_str, shell=True)

	#collecting all the arguments for the parallel execution in a pandas df
	args = None
	skipped = 0
	for infile_R1 in glob.glob(INFOLDER + '/*_R1.fastq.gz'):
		#should we skip this file?
		fn = _create_filenames(infile_R1, OUTFOLDER)
		if os.path.isfile(fn['outfile_R1']) and SKIP_PREVIOUSLY_COMPLETED:
			skipped += 1
			print('Skipping previously processed sample ' + fn['core'])
			continue
		
		#the arguments for the current file. Column order is important,
		#it should match the order for the parallel function, since 
		#the arguments are passed as positionals. Since TRIM_CMD is
		#itself a list we need to pickle it
		args_now = pd.DataFrame({
			'infile_R1' : [infile_R1], 
			'outfolder' : [OUTFOLDER],
			'trim_cmd'  : [pickle.dumps(TRIM_CMD)] 
		})

		#storing in a single df
		args = pd.concat([args, args_now])

		#safeguard
		if len(args) >= MAX_SAMPLES:
			break

	#do we have something to execute?
	if args is None:
		print('All samples skipped, no operation required')
		sys.exit(0)
		
	#executing in parallel using multiprocessing module
	cnt = 0
	with ThreadPool(CORES) as pool:
		for result in pool.starmap(_do_trim, args.itertuples(index = False)):
			#result variable contains the core of the processed sample, but
			#we don't want to flood the main screen (there's many log files)
			#so we just do nothing with it
			cnt += 1

	#closing interface
	print('Samples: ')
	print(' - processed: ' + str(cnt))
	print(' - skipped because of previous runs: ' + str(skipped))
