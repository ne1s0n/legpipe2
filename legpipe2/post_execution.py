#runs the specified script

import os
import common
import subprocess
import shlex

def validate(conf):
	'''validate incoming config parameters from .ini file'''
	#the command, customized by user
	conf['post_execution']['cmd'] = shlex.split(conf['post_execution']['cmd'])

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
	CMD=conf['post_execution']['cmd']
	#execution
	if len(CMD) == 0:
		print('No script specified')
	else:
		print('Running script: ' + ' '.join(CMD))
		subprocess.run(CMD, shell=False)
	
