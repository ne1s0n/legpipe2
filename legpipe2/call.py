#SNP calling with GATK Haplotype Caller module. following 
#the GVCF workflow 

#----------- IMPORTS
import subprocess
import glob
import os
import common
import pandas as pd
import copy
#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 


#----------- SUPPORT FUNCTIONS
def validate(conf):
	'''validate incoming config parameters from .ini file'''
	if not os.path.exists(conf['call']['reference_file']):
		msg = 'Reference file does not exist: ' + conf['call']['reference_file']
		raise FileNotFoundError(msg)
		
	if not os.path.exists(conf['call']['region_lengths_file']):
		msg = 'Regions file does not exist: ' + conf['call']['region_lengths_file']
		raise FileNotFoundError(msg)
		
def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	
	#an editable copy of the dictionary, so that we can manipulate while iterating
	res = copy.deepcopy(conf)

	#these values should be int
	res['call']['ploidy'] = raw_conf['call'].getint('ploidy') 
	res['call']['cores'] = raw_conf['call'].getint('cores') 

	#these values should be boolean
	res['call']['skip_previously_completed'] = raw_conf['call'].getboolean('skip_previously_completed') 
	res['call']['dry_run'] = raw_conf['call'].getboolean('dry_run') 

	#if max_samples is zero it goes to +Infinity
	res['call']['max_samples'] = raw_conf['call'].getint('max_samples')
	if res['call']['max_samples'] == 0:
		res['call']['max_samples'] = float('inf')
	
	return(res)

def _create_filenames(sorted_bam, outfolder):
	'''returns a dictionary with all the derived filenames/folders'''
	
	#storing the input data, just for the record
	res = {'sorted_bam' : sorted_bam}
	res = {'outfolder' : outfolder}
	res = {'outfolder_GVCF' : outfolder + '/GVCF'}
	
	#core sample name, without path and file extension
	res['core'] = os.path.basename(sorted_bam).replace('.gr.sorted.bam', '')
	
	#output files
	res['gvcf']     = res['outfolder_GVCF'] + '/' + os.path.basename(sorted_bam).replace('.gr.sorted.bam', '.g.vcf.gz')
	res['gvcf_log'] = res['outfolder_GVCF'] + '/' + os.path.basename(sorted_bam).replace('.gr.sorted.bam', '.log')
	res['GenomicsDBImport_log'] = outfolder + '/GenomicsDBImport.log'
	res['GenotypeGVCFs_log'] = outfolder + '/GenotypeGVCFs.log'
	
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
	REGION_LENGTHS_FILE = conf['call']['region_lengths_file']

	#tmp folder
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
		#the produced gvcf files, logs
		fn = _create_filenames(infile, OUTFOLDER)
		gvcf_list.append(fn['gvcf'])
		
		#room for result
		cmd = "mkdir -p " + fn['outfolder_GVCF']
		subprocess.run(cmd, shell=True)
		
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
	
	#a mock names dictionary, so that we have the common reference to log files
	fn = _create_filenames('/path/to/fakesample.bam', OUTFOLDER)
	
	#import everything in a genomic db, for each region in the region file
	#https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	#gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	#      -L path/to/intervals.list \
	#      -V data/gvcfs/mother.g.vcf.gz \
	#      -V data/gvcfs/father.g.vcf.gz \
	#      -V data/gvcfs/son.g.vcf.gz \
	#      --tmp-dir /path/to/large/tmp \
	#      --genomicsdb-workspace-path my_database \
	#      --overwrite-existing-genomicsdb-workspace 
	cmd = ['gatk', '--java-options', '-Xmx4g', 'GenomicsDBImport']
	cmd += ['-L', REGION_LENGTHS_FILE]
	cmd += ['--tmp-dir', TMP_FOLDER]
	cmd += ['--genomicsdb-workspace-path', EXPERIMENT]
	cmd += ['--overwrite-existing-genomicsdb-workspace', 'true']
	#calling on all the samples
	for g in gvcf_list:
		cmd += ['-V', g]
	
	#ready to run
	res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['GenomicsDBImport_log'], 'w') as fp:
		fp.write(res.stdout)
		fp.write('------------------------------\n')
	

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
	#res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	#with open(fn['GenotypeGVCFs_log'], "w") as fp:
	#	fp.write(res.stdout)



