#demultiplexing from one batch of merged reads (all sharing the same
#barcodes) to a folder of separated reads 

import subprocess
import glob
import os
import common
import shlex
from common import fn

def validate(conf):
	'''validate incoming config parameters from .ini file'''

	#checking if files/paths exist, forward reads
	if not os.path.exists(conf['demultiplex']['infile_r1']):
		msg = 'Input file demultiplex/infile_r1 does not exist: ' + conf['demultiplex']['infile_r1']
		raise FileNotFoundError(msg)
	#if paired, let's check also the reverse
	if 	conf['demultiplex']['paired']:
		if not os.path.exists(conf['demultiplex']['infile_r2']):
			msg = 'Input file demultiplex/infile_r2 does not exist: ' + conf['demultiplex']['infile_r2']
			raise FileNotFoundError(msg)
	#check the barcode file
	if not os.path.exists(conf['demultiplex']['barcodes']):
		msg = 'Input file demultiplex/barcodes does not exist: ' + conf['demultiplex']['barcodes']
		raise FileNotFoundError(msg)

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#the demultiplex command, customized by user
	conf['demultiplex']['cmd'] = shlex.split(conf['demultiplex']['cmd'])
	return(conf)
	
def demultiplex(conf):
	#interface
	common.print_step_header('demultiplex')
	
	#should we do something?
	RUN_THIS=conf['demultiplex']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#local copies of configuration variables, to have a leaner code
	PAIRED=conf['demultiplex']['paired']
	INFILE_R1=conf['demultiplex']['infile_r1'] # multiplexed data, reads 1
	if PAIRED:
		INFILE_R2=conf['demultiplex']['infile_r2'] # multiplexed data, reads 2
	BARCODES=conf['demultiplex']['barcodes']
	OUTFOLDER=conf['demultiplex']['outfolder']
	DEMUX_CMD=conf['demultiplex']['cmd']

	#derived variables
	STATS_FILE= OUTFOLDER + '/axe_stats.csv'
	LOG_FILE= OUTFOLDER + '/axe_log.txt'

	#room for output
	cmd_str = 'mkdir -p ' + fn(OUTFOLDER)
	subprocess.run(cmd_str, shell=True)
	
	#a little interface
	print(' - infile R1  : ' + INFILE_R1)
	if PAIRED:
		print(' - infile R2  : ' + INFILE_R2)
	print(' - barcodes   : ' + BARCODES)
	print(' - outfolder  : ' + OUTFOLDER)
	print(' - stats file : ' + STATS_FILE)
	print(' - log file   : ' + LOG_FILE)

	if not PAIRED:
		print('\nNow demultiplexing with axe demux, single end')
	else:
		print('\nNow demultiplexing with axe demux, paired end')
	
	#building the command
	cmd = DEMUX_CMD
	cmd += ['-f', INFILE_R1]
	if PAIRED:
		cmd += ['-r', INFILE_R2]
	cmd += ['-F', OUTFOLDER + '/']
	cmd += ['-R', OUTFOLDER + '/']
	cmd += ['-b', BARCODES]
	cmd += ['-t', STATS_FILE]
	
	res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(LOG_FILE, "w") as fp:
		fp.write(res.stdout)
