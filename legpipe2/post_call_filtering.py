#runs the specified script

import os
import common
import subprocess

def validate(conf):
	'''validate incoming config parameters from .ini file'''
	#checking if files/paths exist
	if not os.path.exists(conf['post_call_filtering']['infile']):
		msg = 'Input file does not exist: ' + conf['post_call_filtering']['infile']
		raise FileNotFoundError(msg)
	if not os.path.exists(conf['post_call_filtering']['reference_file']):
		msg = 'Reference file does not exist: ' + conf['post_call_filtering']['reference_file']
		raise FileNotFoundError(msg)
		
def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#files
	conf['post_call_filtering']['infile'] = conf['post_call_filtering']['infolder'] + '/raw_SNPs_haplo.vcf.gz'
	conf['post_call_filtering']['outfile'] = conf['post_call_filtering']['infolder'] + '/filtered_SNPs_haplo.vcf.gz'
	
	#MAF
	conf['post_call_filtering']['min_maf'] = raw_conf['post_call_filtering'].getfloat('min_maf') 

	
	return(conf)

def post_call_filtering(conf):

	#interface
	common.print_step_header('Running post calling filtering')
	
	#should we do something?
	RUN_THIS=conf['post_call_filtering']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)
		
	#config
	INFILE=conf['post_call_filtering']['infile']
	OUTFILE=conf['post_call_filtering']['outfile']
	REFERENCE_FILE=conf['post_call_filtering']['reference_file']
	MIN_MAF=conf['post_call_filtering']['min_maf']
	
	#biallelic SNP
	print(' - select only biallelic SNPs and filter on allele frequency (MAF): ' + str(MIN_MAF))
	#https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
	#https://gatk.broadinstitute.org/hc/en-us/articles/360041850471-VariantFiltration
	#gatk --java-options "-Xmx4g" SelectVariants  \
	#   --reference REFERENCE_FILE \
	#   --variant INFILE \
	#   --output OUTFILE \
	#	--restrictAllelesTo BIALLELIC



#java ‐jar gatk/3.4.0/GenomeAnalysisTK.jar 
#‐T SelectVariants 
#‐R ../ASU/ref_98.fa ‐V biallelic_raw_SNPs_unified.vcf 
#‐‐filterExpression “AF <= 0.100 || AF >= 0.900” 
#‐‐filterName “AF_0.100_0.900” 
#‐o AF_biallelic_raw_SNPs_unified.vcf
