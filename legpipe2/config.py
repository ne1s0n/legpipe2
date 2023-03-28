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
import copy

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
	
	#validating everything, look for invalid values
	_validate_subsample(res)
	
	#interpolating some special cases (e.g. 
	#values that should actually be integer or lists)
	res = _interpolate_rename_reads(res, config)
	res = _interpolate_subsample(res, config)
	res = _interpolate_trim(res, config)
	
	return(res)

def _interpolate_subsample(conf, raw_conf):
	conf['subsample']['seed'] = int(conf['subsample']['seed']) 
	conf['subsample']['reads'] = int(conf['subsample']['reads']) 
	return(conf)

def _interpolate_trim(conf, raw_conf):

	#these values should be boolean
	conf['trim']['dry_run'] = raw_conf['trim'].getboolean('dry_run') 
	conf['trim']['skip_previously_completed'] = raw_conf['trim'].getboolean('skip_previously_completed') 

	#these values should be int
	conf['trim']['cores'] = raw_conf['trim'].getint('cores') 

	#if max_samples is zero it goes to +Infinity
	conf['trim']['max_samples'] = raw_conf['trim'].getint('max_samples')
	if conf['trim']['max_samples'] == 0:
		conf['trim']['max_samples'] = float('inf')
	
	return(conf)
	
def _interpolate_rename_reads(conf, raw_conf):
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

def _validate_subsample(conf):
	if conf['subsample']['tool'] not in ['cat', 'seqtk']:
		msg = 'Config parameter subsample/tool must be in [cat, seqtk], found : ' + conf['subsample']['tool']  
		raise ValueError(msg)
