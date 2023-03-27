import configparser
def read_config(infile):
	'''
	Reads the config file in .ini format, returns a dictionary
	'''
	
	config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
	config.read(infile)
	
	return(config)
