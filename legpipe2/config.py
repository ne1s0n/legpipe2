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
import importlib

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
			#some on-the-fly interpolations for common values in all sections
			if key == 'run_this':
				res[section][key] = config[section].getboolean(key)	
			elif key == 'cmd':
				res[section][key] = config[section][key].split('\n')
			else:
				#just copying the value
				res[section][key] = config[section][key]
	
	#sections correspond to modules (e.g. operations)
	for m in config.sections():
		#should we import/interpolate/validate?
		if res[m]['run_this']:
			mymodule = importlib.import_module(m)
			res = mymodule.interpolate(res, config)
			mymodule.validate(res)
	
	return(res)

