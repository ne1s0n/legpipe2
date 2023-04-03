#SNP calling with GATK Haplotype Caller module. following 
#the GVCF workflow 

#----------- IMPORTS
import subprocess
import glob
import os
import common

#----------- SUPPORT FUNCTIONS
def validate(conf):
	'''validate incoming config parameters from .ini file'''
	if not os.path.exists(conf['call']['reference_file']):
		msg = 'Reference file does not exist: ' + conf['genome_index']['reference_file']
		raise FileNotFoundError(msg)

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	#these values should be int
	conf['call']['ploidy'] = raw_conf['call'].getint('ploidy') 

	#if max_samples is zero it goes to +Infinity
	conf['call']['max_samples'] = raw_conf['call'].getint('max_samples')
	if conf['call']['max_samples'] == 0:
		conf['call']['max_samples'] = float('inf')
	return(conf)

#----------- CALL (public) FUNCTION
def call(conf):
	#interface
	common.print_step_header('call SNPs')
	
	#should we do something?
	RUN_THIS=conf['call']['run_this']
	if not RUN_THIS:
		print('SKIPPED')
		return(None)

	#currently used params
	INFOLDER=conf['call']['infolder']
	OUTFOLDER=conf['call']['outfolder']
	REFERENCE_FILE=conf['call']['reference_file']
	PLOIDY=conf['call']['ploidy']
	MAX_SAMPLES=conf['call']['max_samples']
	EXPERIMENT=conf['call']['experiment']
	TMP_FOLDER=conf['call']['tmp_folder']

	#a subfolder for GATK gvcf files
	OUTFOLDER_GVCF = OUTFOLDER + '/GVCF'

	#room for output, tmp
	cmd = "mkdir -p " + OUTFOLDER_GVCF
	subprocess.run(cmd, shell=True)
	cmd = "mkdir -p " + TMP_FOLDER
	subprocess.run(cmd, shell=True)

	#keeping track of the produced .g.vcf.gz files
	gvcf_list = []

	#------------ HaplotypeCaller
	#for each input bam
	for infile in glob.glob(INFOLDER + '/*.gr.sorted.bam'):
		#the produced gvcf file, log
		gvcf = OUTFOLDER_GVCF + '/' + os.path.basename(infile).replace('.gr.sorted.bam', '.g.vcf.gz')
		gvcf_list.append(gvcf)
		gvcf_log = OUTFOLDER_GVCF + '/' + os.path.basename(infile).replace('.gr.sorted.bam', '.log')
		
		#https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller
		#gatk --java-options "-Xmx4g" HaplotypeCaller  \
		#   -R Homo_sapiens_assembly38.fasta \
		#   -I input.bam \
		#   -O output.g.vcf.gz \
		#   -ERC GVCF \
		#   -ploidy PLOIDY \
		#   --min-pruning 1 #default is two \
		#   -stand-call-conf 30 #default 
		
		cmd = ['gatk', '--java-options', '"-Xmx4g"', 'HaplotypeCaller']
		cmd += ['-ERC', 'GVCF', '--min-pruning', '1', '-stand-call-conf', '30'] 
		cmd += ['-ploidy', str(PLOIDY)] 
		cmd += ['-R', REFERENCE_FILE]
		cmd += ['-I', infile]
		cmd += ['-O', gvcf]
		print(cmd)
		print(' '.join(cmd))
		res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
		with open(gvcf_log, "w") as fp:
			fp.write(res.stdout)
		
		if len(gvcf_list) >= MAX_SAMPLES:
			break

	#------------ GenomicsDBImport
	#import everything in a genomic db
	#https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	#gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	#      -V data/gvcfs/mother.g.vcf.gz \
	#      -V data/gvcfs/father.g.vcf.gz \
	#      -V data/gvcfs/son.g.vcf.gz \
	#      --genomicsdb-workspace-path my_database \
	#      --tmp-dir /path/to/large/tmp \
	cmd = ['gatk', '--java-options', '"-Xmx4g"', 'GenomicsDBImport']
	for g in gvcf_list:
		cmd += ['-V', g]
	cmd += ['--genomicsdb-workspace-path', EXPERIMENT]
	cmd += ['--tmp-dir', TMP_FOLDER]

	print(cmd)
	subprocess.run(cmd, shell=True)

	#------------ GenotypeGVCFs
	#joint variant calling
	#https://gatk.broadinstitute.org/hc/en-us/articles/9570489472411-GenotypeGVCFs
	# gatk --java-options "-Xmx4g" GenotypeGVCFs \
	#   -ploidy PLOIDY \
	#   -R Homo_sapiens_assembly38.fasta \
	#   -V gendb://my_database \
	#   -O output.vcf.gz \
	#   --tmp-dir /path/to/large/tmp
	cmd = ['gatk', '--java-options', '"-Xmx4g"', 'GenotypeGVCFs']
	cmd += ['-ploidy', str(PLOIDY)] 
	cmd += ['-R',  REFERENCE_FILE] 
	cmd += ['-V',  'gendb://' + EXPERIMENT]
	cmd += ['-O',  OUTFOLDER + '/raw_SNPs_haplo.vcf.gz']
	cmd += ['--tmp-dir' + TMP_FOLDER]
	print(cmd)
	subprocess.run(cmd, shell=True)
