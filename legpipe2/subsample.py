#subsampling (either using seqtk or cat) all the fastq files
#named *_R1.fastq.gz and *_R2.fastq.gz present in the specified infolder. 
#Reading and compressing via gzip

import subprocess
import glob
import os
import common

def subsample_validate(conf):
	if conf['subsample']['tool'] not in ['cat', 'seqtk']:
		msg = 'Config parameter subsample/tool must be in [cat, seqtk], found : ' + conf['subsample']['tool']  
		raise ValueError(msg)
		
def subsample(conf):
	#basic folder
	INFOLDER=conf['subsample']['infolder']
	OUTFOLDER=conf['subsample']['outfolder']
	#random seed, same for everybody, so that to maintain paired ends reads
	SEED=conf['subsample']['seed']
	#how many reads are we going to keep
	READS=conf['subsample']['reads']
	#tool to be used
	TOOL=conf['subsample']['tool']
	RUN_THIS=conf['subsample']['run_this']

	#interface
	common.print_step_header('subsample')
	
	#should we do something?
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#room for output
	cmd_str = "mkdir -p " + OUTFOLDER
	subprocess.run(cmd_str, shell=True)

	#a bit of interface
	print("Infolder:  " + INFOLDER)
	print("Outfolder: " + OUTFOLDER)

	#retrieving the list of files
	infiles = glob.glob(INFOLDER + '/*_R1.fastq.gz') + glob.glob(INFOLDER + '/*_R2.fastq.gz')

	#for each infile, call either seqtk or just cat
	for infile in infiles:
		#build the outfile name
		outfile = OUTFOLDER + '/' + os.path.basename(infile) 
		print ('from ' + infile + ' to ' + outfile)
		
		#building the command
		if TOOL == 'seqtk':
			#example: seqtk sample -s 123 plate12_demuxed/1A10_R1.fastq.gz 1000          | gzip > foo.fq.gz
			cmd_str = "seqtk sample -s " + str(SEED) + ' ' + infile + ' ' + str(READS) + ' | gzip > ' + outfile
		else:		
			#example: zcat plate12_demuxed/1A10_R1.fastq.gz | head -n lines | gzip > foo.fq.gz
			cmd_str = 'zcat ' + infile + ' | head -n ' + str(READS) + ' | gzip > ' + outfile
		
		#runnning the command
		subprocess.run(cmd_str, shell=True)
