#the .fastq.gz files so far are in two or mode folders (usually divided by plate).
#We want to create a single folder with all the files - not copies, just sybolic links

import subprocess
import glob
import os
import common
import copy

def validate(conf):
	'''validate incoming config parameters from .ini file'''

def interpolate(conf, raw_conf):
	'''two groups of values should become lists with the same internal order
		- INFOLDER_X (e.g. INFOLDER_1, INFOLDER_2, ...)
		- OUTFILE_PREFIX_X (e.g. OUTFILE_PREFIX_1, OUTFILE_PREFIX_2, ...) 
	'''
	#an editable copy of the dictionary, so that we can manipulate while iterating
	res = copy.deepcopy(conf)
	res['rename_reads']['infolders'] = []
	res['rename_reads']['outfile_prefixes'] = []
	
	for key in conf['rename_reads']:
		if not key.startswith('infolder_'):
			continue
		#found an infolder
		res['rename_reads']['infolders'].append(conf['rename_reads'][key])
		
		#looking for the corresponding outfile_prefix
		new_key = key.replace('infolder_', 'outfile_prefix_')
		res['rename_reads']['outfile_prefixes'].append(conf['rename_reads'][new_key])
		
		#we can now remove those entries
		res['rename_reads'].pop(key)
		res['rename_reads'].pop(new_key)
	
	return(res)

def rename_reads(conf):
	#interface
	common.print_step_header('rename_reads')
	
	#should we do something?
	RUN_THIS=conf['rename_reads']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#local copies of configuration variables, to have a leaner code
	INFOLDERS=conf['rename_reads']['infolders'] #where is the original data
	OUTFILE_PREFIXES=conf['rename_reads']['outfile_prefixes'] #prepend to each fastq
	OUTFOLDER=conf['rename_reads']['outfolder']
	
	#room for output
	cmd_str = "mkdir -p " + OUTFOLDER
	subprocess.run(cmd_str, shell=True)

	#for each infolder
	for i in range(len(INFOLDERS)):
		infolder = INFOLDERS[i]
		#a bit of interface
		print("Infolder:  " + infolder)
		
		#for each fastq file
		infiles = glob.glob(infolder + '/*.fastq.gz')
		for infile in infiles:
			#creating the links
			cmd_str = "ln --symbolic " + infile + ' ' + OUTFOLDER + '/' + OUTFILE_PREFIXES[i] + os.path.basename(infile) 
			subprocess.run(cmd_str, shell=True)
