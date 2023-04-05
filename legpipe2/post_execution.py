#runs the specified script

import os
import common
import subprocess

def validate(conf):
	'''validate incoming config parameters from .ini file'''

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	return(conf)

def post_execution(conf):

	#interface
	common.print_step_header('Running post execution script')
	
	#should we do something?
	RUN_THIS=conf['post_execution']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)
	#config
	CMD=conf['post_execution']['cmd'][0]
	#execution
	if CMD == '':
		print('No script specified')
	else:
		print('Running script: ' + CMD)
		subprocess.run(CMD, shell=True)
	
