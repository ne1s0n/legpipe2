#script for collecting aligning statistics from bowtie2 logs

INFOLDER='/data/processed/aligned'
OUTFILE='/data/processed/aligned/aligning_stats.csv'

INFOLDER='/home/nelson/tmp/align'
OUTFILE='/home/nelson/tmp/align/aligning_stats.csv'

#----------- IMPORT
import glob
import pandas as pd

#----------- SUPPORT FUNCTIONS
def parse_bowtie2_log(logfile):
	#TODO write this function
	return(pd.DataFrame({'a' : [1], 'b' : ['fdasfdsa']}))

#----------- ACTUAL SCRIPT
#room for collecting result
res = None

#for each input logfile
for infile in glob.glob(INFOLDER + '/*.bowtie2_align.log'):
	#parse the file
	tmp = parse_bowtie2_log(infile)
	
	#collate the data
	res = pd.concat([res, tmp])

#save the results
res.to_csv(OUTFILE, index=False)
