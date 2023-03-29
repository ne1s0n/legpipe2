import configparser
import copy

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
	
	print(conf['demultiplex'])
	
	#done
	return(conf)

