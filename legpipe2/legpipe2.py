import configparser
import copy

from config import read_config
from rename_reads import rename_reads


def Legpipe2(infile):
	conf = read_config(infile)
	return(conf)

