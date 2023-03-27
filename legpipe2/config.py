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
			#some values should be interpolated on the fly
			if key == 'run_this':
				res[section][key] = config[section].getboolean(key)	
			else:
				#just copying the value
				res[section][key] = config[section][key]
	
	#validating everything, look for invalid values
	_validate_subsample(res)
	
	#interpolating some special cases (e.g. 
	#values that should actually be lists)
	res = _interpolate_rename_reads(res)
	
	return(res)
	
def _interpolate_rename_reads(conf):
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
