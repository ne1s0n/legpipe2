import configparser
import copy
import sys

#legpipe2 modules
from config import read_config
from rename_reads import rename_reads
from subsample import subsample
from trim import trim
from demultiplex import demultiplex

def Legpipe2(infile):
	#reading the config in
	conf = read_config(infile)
	
	#the pipeline steps
	demultiplex(conf)
	rename_reads(conf)
	subsample(conf)
	trim(conf)
	
	#print(conf['rename_reads'])
	
	#done
	return(conf)

if __name__ == "__main__":
	print('Welcome to Legpipe2!')
	
	#very simple command line interface
	if len(sys.argv) != 2:
		print('Usage: python3 path/to/legpipe2.py config_file.ini')
		exit(0)
	
	#invoking legpipe
	Legpipe2(sys.argv[1])
