#demultiplexing from one batch of merged reads (all sharing the same
#barcodes) to a folder of separated reads 

import subprocess
import glob
import os
import common

def demultiplex_validate(conf):
	#checking if files/paths exist
	if not os.path.exists(conf['demultiplex']['infile_r1']):
		msg = 'Input file demultiplex/infile_r1 does not exist: ' + conf['demultiplex']['infile_r1']
		raise FileNotFoundError(msg)
	if not os.path.exists(conf['demultiplex']['infile_r2']):
		msg = 'Input file demultiplex/infile_r2 does not exist: ' + conf['demultiplex']['infile_r2']
		raise FileNotFoundError(msg)
	if not os.path.exists(conf['demultiplex']['barcodes']):
		msg = 'Input file demultiplex/barcodes does not exist: ' + conf['demultiplex']['barcodes']
		raise FileNotFoundError(msg)

def demultiplex(conf):
	#local copies of configuration variables, to have a leaner code
	RUN_THIS=conf['demultiplex']['run_this']
	INFILE_R1=conf['demultiplex']['infile_r1'] # multiplexed data, reads 1
	INFILE_R2=conf['demultiplex']['infile_r2'] # multiplexed data, reads 2
	BARCODES=conf['demultiplex']['barcodes']
	OUTFOLDER=conf['demultiplex']['outfolder']
	DEMUX_CMD=conf['demultiplex']['demux_cmd']
	DEMUX_PARAMS=conf['demultiplex']['demux_params']

	#derived variables
	STATS_FILE= OUTFOLDER + '/axe_stats.csv'
	LOG_FILE= OUTFOLDER + '/axe_log.txt'
	
	#interface
	common.print_step_header('demultiplex')
	
	#should we do something?
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#room for output
	cmd_str = "mkdir -p " + OUTFOLDER
	subprocess.run(cmd_str, shell=True)
	
	#a little interface
	print(' - infile R1  : ' + INFILE_R1)
	print(' - infile R2  : ' + INFILE_R2)
	print(' - barcodes   : ' + BARCODES)
	print(' - outfolder  : ' + OUTFOLDER)
	print(' - stats file : ' + STATS_FILE)
	print(' - log file   : ' + LOG_FILE)

	#paired end demux
	print('\nNow demultiplexing with axe demux, paired end')
	
	#building the command
	cmd = [DEMUX_CMD, DEMUX_PARAMS]
	cmd += ['-f' , INFILE_R1]
	cmd += ['-r' , INFILE_R2]
	cmd += ['-F' , OUTFOLDER + '/']
	cmd += ['-R' , OUTFOLDER + '/']
	cmd += ['-b' , BARCODES]
	cmd += ['-t' , STATS_FILE]
	
	print(cmd)

	res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(LOG_FILE, "w") as fp:
		fp.write(res.stdout)
