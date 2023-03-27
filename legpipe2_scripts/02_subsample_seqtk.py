#subsampling (using seqtk, which should be installed) all the fastq files
#named *_R1.fastq.gz and *_R2.fastq.gz present in the specified folder. 
#Reading and compressing via gzip (which should also be installed)

INFOLDER='/data/raw/all_fastqs'
SEED=123 #random seed, same for everybody, so that to maintain paired ends reads
READS=10000 #how many reads are we going to keep

import subprocess
import glob
import os

#room for output
outfolder = INFOLDER + '_subsampled'
cmd_str = "mkdir -p " + outfolder
subprocess.run(cmd_str, shell=True)

#a bit of interface
print("======================")
print("Infolder:  " + INFOLDER)
print("Outfolder: " + outfolder)

#retrieving the list of files
infiles = glob.glob(INFOLDER + '/*_R1.fastq.gz') + glob.glob(INFOLDER + '/*_R2.fastq.gz')

#for each infile, call seqtk
for infile in infiles:
	#build the outfile name
	outfile = outfolder + '/' + os.path.basename(infile) 
	print ('from ' + infile + ' to ' + outfile)
	
	#building the seqtk command
	#example: seqtk sample -s 123 plate12_demuxed/1A10_R1.fastq.gz 1000          | gzip > foo.fq.gz
	cmd_str = "seqtk sample -s " + str(SEED) + ' ' + infile + ' ' + str(READS) + ' | gzip > ' + outfile
	subprocess.run(cmd_str, shell=True)
