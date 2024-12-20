#runs the specified script

import os
import common
import subprocess

def validate(conf, runtime = False):
	'''validate incoming config parameters from .ini file. Since some 
	argument do exist only when the previous steps are executed, we have the runtime
	parameter which is true only when post_call_filtering() is called'''
	
	#checking if files/paths exist
	if runtime:
		if not os.path.exists(conf['post_call_filtering']['infile']):
			msg = 'Input file does not exist: ' + conf['post_call_filtering']['infile']
			raise FileNotFoundError(msg)
		if not os.path.exists(conf['post_call_filtering']['infile'] + '.tbi'):
			msg = 'Index file for input vcf does not exist: ' + conf['post_call_filtering']['infile'] + '.tbi'
			raise FileNotFoundError(msg)
		if not os.path.exists(conf['post_call_filtering']['reference_file']):
			msg = 'Reference file does not exist: ' + conf['post_call_filtering']['reference_file']
			raise FileNotFoundError(msg)
		
def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#files
	conf['post_call_filtering']['infile'] = conf['post_call_filtering']['infolder'] + '/raw_SNPs_haplo.vcf.gz'
	
	#MAF
	conf['post_call_filtering']['min_maf'] = raw_conf['post_call_filtering'].getfloat('min_maf') 
	
	#Mapping Quality MQ
	conf['post_call_filtering']['min_mq'] = raw_conf['post_call_filtering'].getfloat('min_mq') 

	
	return(conf)

def post_call_filtering(conf):

	#interface
	common.print_step_header('Post calling filtering')
	
	#should we do something?
	RUN_THIS=conf['post_call_filtering']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)
	
	#let's check if all required files are here
	validate(conf=conf, runtime=True) 
		
	#config
	INFILE=conf['post_call_filtering']['infile']
	OUTFOLDER=conf['post_call_filtering']['outfolder']
	REFERENCE_FILE=conf['post_call_filtering']['reference_file']
	MIN_MAF=conf['post_call_filtering']['min_maf']
	MIN_MQ=conf['post_call_filtering']['min_mq']
	
	#derived conf
	OUTFILE = OUTFOLDER + '/filtered_SNPs_haplo.vcf.gz'
	TMPFILE = OUTFOLDER + '/tmp.vcf.gz'
	LOGFILE = OUTFOLDER + '/post_call_filtering.log'
	
	#room for output
	cmd_str = "mkdir -p " + common.fn(OUTFOLDER)
	subprocess.run(cmd_str, shell=True)

	
	#FILTER ONE: biallelic
	#https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
	#https://gatk.broadinstitute.org/hc/en-us/articles/360041850471-VariantFiltration
	#gatk --java-options "-Xmx4g" SelectVariants  \
	#   --reference REFERENCE_FILE \
	#   --variant INFILE \
	#   --output OUTFILE \
	#	--restrictAllelesTo BIALLELIC
	print(' - select only biallelic SNPs')
	cmd = ['gatk', '--java-options', '-Xmx4g', 'SelectVariants']
	cmd += ['--reference', REFERENCE_FILE]
	cmd += ['--variant', INFILE]
	cmd += ['--output', TMPFILE]
	cmd += ['--restrict-alleles-to', 'BIALLELIC']
	with open(LOGFILE, "w") as fp:
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)
	
	#FILTER TWO: MAF, mapping quality
	#gatk --java-options "-Xmx4g" VariantFiltration  \
	#   --reference $REFERENCE_FILE \
	#   --variant $OUTFILE1 \
	#   --output $OUTFILE2 \
	#   --filter-expression "AF <= 0.050 || AF >= 0.950" \
	#   --filter-name "AF_0.05_0.95"
	#   --filter-expression "MQ <= 40 \
	#   --filter-name "MQ_40"
	print(' - filter on minor allele frequency (MAF): ' + str(MIN_MAF))
	print(' - filter on minimum mapping quality: ' + str(MIN_MQ))
	cmd = ['gatk', '--java-options', '-Xmx4g', 'VariantFiltration']
	cmd += ['--reference', REFERENCE_FILE]
	cmd += ['--variant', TMPFILE]
	cmd += ['--output', OUTFILE]
	cmd += ['--filter-expression', 'AF <= ' + str(MIN_MAF) + ' || AF >= ' + str(1-MIN_MAF)]
	cmd += ['--filter-name', 'AF_'+ str(MIN_MAF) + '_' + str(1-MIN_MAF)]
	cmd += ['--filter-expression', 'MQ < ' + str(MIN_MQ)]
	cmd += ['--filter-name', 'MQ_'+ str(MIN_MQ)]
	with open(LOGFILE, "a") as fp:
		fp.write(' '.join(cmd) + '\n')
		subprocess.run(cmd, shell=False, stdout=fp, stderr=subprocess.STDOUT, text=True)

	#cleanup
	subprocess.run(['rm', TMPFILE], shell=False)
