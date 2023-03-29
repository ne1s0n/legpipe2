#the .fastq.gz files so far are in two or mode folders (usually divided by plate).
#We want to create a single folder with all the files - not copies, just sybolic links

import subprocess
import glob
import os
import common

def rename_reads(conf):
	#local copies of configuration variables, to have a leaner code
	INFOLDERS=conf['rename_reads']['infolders'] #where is the original data
	OUTFILE_PREFIXES=conf['rename_reads']['outfile_prefixes'] #prepend to each fastq
	OUTFOLDER=conf['rename_reads']['outfolder']
	RUN_THIS=conf['rename_reads']['run_this']
	
	#interface
	common.print_step_header('rename_reads')
	
	#should we do something?
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

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
