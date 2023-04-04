#SNP calling with GATK Haplotype Caller module. following 
#the GVCF workflow 

#----------- IMPORTS
import subprocess
import glob
import os
import common
import pandas as pd
#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 


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
	conf['call']['cores'] = raw_conf['call'].getint('cores') 

	#these values should be boolean
	conf['call']['skip_previously_completed'] = raw_conf['call'].getboolean('skip_previously_completed') 
	conf['call']['dry_run'] = raw_conf['call'].getboolean('dry_run') 

	#if max_samples is zero it goes to +Infinity
	conf['call']['max_samples'] = raw_conf['call'].getint('max_samples')
	if conf['call']['max_samples'] == 0:
		conf['call']['max_samples'] = float('inf')
	return(conf)

def _create_filenames(sorted_bam, outfolder):
	'''returns a dictionary with all the derived filenames/folders'''
	
	#storing the input data, just for the record
	res = {'sorted_bam' : sorted_bam}
	res = {'outfolder' : outfolder}
	
	#core sample name, without path and file extension
	res['core'] = os.path.basename(sorted_bam).replace('.gr.sorted.bam', '')
	
	#output files
	res['gvcf']     = outfolder + '/' + os.path.basename(sorted_bam).replace('.gr.sorted.bam', '.g.vcf.gz')
	res['gvcf_log'] = outfolder + '/' + os.path.basename(sorted_bam).replace('.gr.sorted.bam', '.log')

	return(res)



def _haplotypecaller(ploidy, reference_file, infile, outfile, logfile, dry_run):
	#https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller
	#gatk --java-options "-Xmx4g" HaplotypeCaller  \
	#   -R Homo_sapiens_assembly38.fasta \
	#   -I input.bam \
	#   -O output.g.vcf.gz \
	#   -ERC GVCF \
	#   -ploidy PLOIDY \
	#   --min-pruning 1 #default is two \
	#   -stand-call-conf 30 #default 
	cmd = ['gatk', '--java-options', '-Xmx4g', 'HaplotypeCaller']
	cmd += ['-ERC', 'GVCF', '--min-pruning', '1', '-stand-call-conf', '30'] 
	cmd += ['-ploidy', ploidy] 
	cmd += ['-R', reference_file]
	cmd += ['-I', infile]
	cmd += ['-O', outfile]
	if dry_run:
		cmd += ['--dry-run']
	res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(logfile, "w") as fp:
		fp.write(res.stdout)

#----------- CALL (public) FUNCTION
def call(conf):
	
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
	CORES=conf['call']['cores']
	MAX_SAMPLES=conf['call']['max_samples']
	EXPERIMENT=conf['call']['experiment']
	TMP_FOLDER=conf['call']['tmp_folder']
	SKIP_PREVIOUSLY_COMPLETED=conf['call']['skip_previously_completed']
	DRY_RUN=conf['call']['dry_run']

	#a subfolder for GATK gvcf files
	OUTFOLDER_GVCF = OUTFOLDER + '/GVCF'

	#room for output, tmp
	cmd = "mkdir -p " + OUTFOLDER_GVCF
	subprocess.run(cmd, shell=True)
	cmd = "mkdir -p " + TMP_FOLDER
	subprocess.run(cmd, shell=True)

	#list of all .g.vcf.gz files produced/available
	gvcf_list = []
	#list of the arguments for parallel execution
	args = None
	#number of skipped samples (because already called in previous runs) 
	skipped = 0

	#------------ HaplotypeCaller
	#interface
	common.print_step_header('call SNPs - HaplotypeCaller')
	#for each input bam
	for infile in glob.glob(INFOLDER + '/*.gr.sorted.bam'):
		#the produced gvcf file, log
		fn = _create_filenames(infile, OUTFOLDER)
		gvcf_list.append(fn['gvcf'])
		
		#should we skip this file?
		if os.path.isfile(fn['gvcf']) and SKIP_PREVIOUSLY_COMPLETED:
			skipped += 1
			print('Skipping previously processed sample ' + fn['gvcf'])
			continue
		
		#the arguments for the current file. Column order is important,
		#it should match the order for the parallel function, since 
		#the arguments are passed as positionals
		args_now = pd.DataFrame({
			'ploidy' : [str(PLOIDY)], 
			'reference_file' : [REFERENCE_FILE],
			'infile'  : [infile],
			'outfile' : fn['gvcf'],
			'logfile' : fn['gvcf_log'],
			'dry_run' : [DRY_RUN]
		})
		
		#storing in a single df
		args = pd.concat([args, args_now])

		#safeguard
		if len(args) >= MAX_SAMPLES:
			break

	#do we have something to execute?
	if args is None:
		print('All samples skipped/inherited from previous runs')
	else:		
		#executing in parallel using multiprocessing module
		with ThreadPool(CORES) as pool:
			for result in pool.starmap(_haplotypecaller, args.itertuples(index = False)):
				#nothing to do here (result is None)
				pass

		#closing interface
		print('Samples: ')
		print(' - samples called: ' + str(len(args)))
		print(' - inherited from previous runs: ' + str(skipped))
	
	#------------ GenomicsDBImport
	#interface
	common.print_step_header('call SNPs - GenomicsDBImport')
	#import everything in a genomic db
	#https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	#gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	#      -V data/gvcfs/mother.g.vcf.gz \
	#      -V data/gvcfs/father.g.vcf.gz \
	#      -V data/gvcfs/son.g.vcf.gz \
	#      --genomicsdb-workspace-path my_database \
	#      --tmp-dir /path/to/large/tmp \
	cmd = ['gatk', '--java-options', '-Xmx4g', 'GenomicsDBImport']
	for g in gvcf_list:
		cmd += ['-V', g]
	cmd += ['--genomicsdb-workspace-path', EXPERIMENT]
	cmd += ['--tmp-dir', TMP_FOLDER]
	if DRY_RUN:
		cmd += ['--dry-run']
	subprocess.run(cmd, text=True)

	#------------ GenotypeGVCFs
	#interface
	common.print_step_header('call SNPs - GenotypeGVCFs')
	#joint variant calling
	#https://gatk.broadinstitute.org/hc/en-us/articles/9570489472411-GenotypeGVCFs
	# gatk --java-options "-Xmx4g" GenotypeGVCFs \
	#   -ploidy PLOIDY \
	#   -R Homo_sapiens_assembly38.fasta \
	#   -V gendb://my_database \
	#   -O output.vcf.gz \
	#   --tmp-dir /path/to/large/tmp
	cmd = ['gatk', '--java-options', '-Xmx4g', 'GenotypeGVCFs']
	cmd += ['-ploidy', str(PLOIDY)] 
	cmd += ['-R',  REFERENCE_FILE] 
	cmd += ['-V',  'gendb://' + EXPERIMENT]
	cmd += ['-O',  OUTFOLDER + '/raw_SNPs_haplo.vcf.gz']
	cmd += ['--tmp-dir' + TMP_FOLDER]
	if DRY_RUN:
		cmd += ['--dry-run']
	subprocess.run(cmd, text=True)


