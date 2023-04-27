#from the filtered vcf file, output data in several matrices
#fornat


#https://www.reneshbedre.com/blog/vcf-fields.html

#GT:AD:DP:GQ:PGT:PID:PL:PS       
#0/0:4,0:4:12:.:.:0,12,117       
#0/0:1,0:1:3:.:.:0,3,15  
#0/0:0,0:0:0:.:.:0,0,0   
#0/0:5,0:5:15:.:.:0,15,197       
#0/0:3,0:3:9:.:.:0,9,99  
#0/0:5,0:5:0:.:.:0,0,113 
#0/0:5,0:5:15:.:.:0,15,156

import common

def validate(conf):
	'''validate incoming config parameters from .ini file, plus environmental variables'''

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	return(conf)

def output_matrices(conf):
	'''main function for outputing matrices'''
	
	#interface
	common.print_step_header('Output matrices')
	
	#should we do something?
	RUN_THIS=conf['output_matrices']['run_this']

	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	print('WARNING: this is a stub function')
	return(None)
