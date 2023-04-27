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
	conf['post_call_filtering']['tmpfile'] = conf['post_call_filtering']['infolder'] + '/tmp.vcf.gz'
	conf['post_call_filtering']['logfile'] = conf['post_call_filtering']['infolder'] + '/post_call_filterin.log'
	
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
	TMPFILE=conf['post_call_filtering']['tmpfile']
	LOGFILE=conf['post_call_filtering']['logfile']
	REFERENCE_FILE=conf['post_call_filtering']['reference_file']
	MIN_MAF=conf['post_call_filtering']['min_maf']
	
	#biallelic SNP
	print(' - select only biallelic SNPs and filter on allele frequency (MAF): ' + str(MIN_MAF))
	
	
	#FILTER ONE: biallelic
	#https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
	#https://gatk.broadinstitute.org/hc/en-us/articles/360041850471-VariantFiltration
	#gatk --java-options "-Xmx4g" SelectVariants  \
	#   --reference REFERENCE_FILE \
	#   --variant INFILE \
	#   --output OUTFILE \
	#	--restrictAllelesTo BIALLELIC
	cmd = ['gatk', '--java-options', '-Xmx4g', 'SelectVariants']
	cmd += ['--reference', REFERENCE_FILE]
	cmd += ['--variant', INFILE]
	cmd += ['--output', TMPFILE]
	cmd += ['--restrict-alleles-to', 'BIALLELIC']
	with open(LOGFILE, "w") as fp:
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#FILTER TWO: MAF
	#gatk --java-options "-Xmx4g" VariantFiltration  \
	#   --reference $REFERENCE_FILE \
	#   --variant $OUTFILE1 \
	#   --output $OUTFILE2 \
	#   --filter-expression "AF <= 0.050 || AF >= 0.950" \
	#   --filter-name "AF_0.05_0.95"
	cmd = ['gatk', '--java-options', '-Xmx4g', 'VariantFiltration']
	cmd += ['--reference', REFERENCE_FILE]
	cmd += ['--variant', TMPFILE]
	cmd += ['--output', OUTFILE]
	cmd += ['--filter-expression', 'AF <= ' + str(MIN_MAF) + ' || AF >= ' + str(1-MIN_MAF)]
	cmd += ['--filter-name', 'AF_'+ str(MIN_MAF) + '_' + str(1-MIN_MAF)]
	with open(LOGFILE, "a") as fp:
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)


	#cleanup
	subprocess.run(['rm', TMPFILE], shell=False)







