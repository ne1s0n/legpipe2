import configparser
import copy

from config import read_config
from rename_reads import rename_reads
from subsample import subsample

def Legpipe2(infile):
	#reading the config in
	conf = read_config(infile)
	
	#the pipeline steps
	rename_reads(conf)
	subsample(conf)
	
	#done
	return(conf)

