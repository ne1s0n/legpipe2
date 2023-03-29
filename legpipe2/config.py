#This module implements .ini config file parsing, plus several data transformation and
#validation. The core function is read_config(), which requires the .ini infile and
#returns a two (or more) levels dictionary (not a configparser). First level keys are 
#the .ini sections, second level are values.  
#We use the same notation as configparser objects, in particular we refer to "interpolate"
#as the operation of tweaking data, e.g. typing (so that values are actually booleans, int
#and so forth). The main extra functionality injected is the possibility of having lists,
#which are not supported by configparser. When possible we use however configparser
#native interpolation/validation mechanisms

import configparser
import subsample
import trim
import rename_reads
import demultiplex

def read_config(infile):
	'''
	Reads the config file in .ini format, parse some data so e.g. there's lists
	and not many keys, returns a dictionary
	'''
	
	config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
	config.read(infile)
	
	#loading everything
	res = {}
	for section in config.sections():
		res[section] = {}
		for key in config[section]:
			#this value is common to every section and we interpolate it on the fly
			if key == 'run_this':
				res[section][key] = config[section].getboolean(key)	
			else:
				#just copying the value
				res[section][key] = config[section][key]
	
	#interpolating some special cases (e.g. 
	#values that should actually be integer or lists)
	res = rename_reads.interpolate(res, config)
	res = subsample.interpolate(res, config)
	res = trim.interpolate(res, config)
	res = demultiplex.interpolate(res, config)
	
	#validating everything, look for invalid values
	subsample.validate(res)
	demultiplex.validate(res)
	
	return(res)

