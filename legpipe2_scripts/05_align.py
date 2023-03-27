#align reads to the indexed genome. Tools used: bowtie2, samtools, picard
#Notes:
# - from INFOLDER we take all forward and reverse reads (.fastq.gz)
# - the genome is expected to be indexed as per 05_genome_index.py script

#for each sample there's a list of commands:
	#bowtie2 -x ref_98 -1 sample1.1.fq -2 sample1.2.fq -S sample1.sam
	#samtools view -bS sample1.sam > sample1.bam
	#java -jar picard/2.4.1/ValidateSamFile.jar INPUT= sample1.bam
	#java -jar picard/2.4.1/AddOrReplaceReadGroups.jar I=sample1.bam O=sample1.Gr.bam LB=Whatever PL=Illumina PU=Whatever SM=sample1
	#samtools/1.3.1/samtools sort sample1.Gr.bam sample1.Gr.sorted.bam
	#samtools/1.3.1/samtools index sample1.Gr.sorted.bam


#----------- CONFIG 
REFERENCE_FILE='/data/reference/zm-4.genome.largest_chromosomes.fasta.gz'
PICARD_COMMAND='java -jar /home/ubuntu/software/picard.jar'
INFOLDER='/data/processed/trimmed_fastp'
OUTFOLDER='/data/processed/aligned_fastp'
CORES=7 #number of parallel threads, put to zero to run as many as possible
MAX_SAMPLES=float("inf") #maximum number of samples to process, for test purposes
SKIP_PREVIOUSLY_COMPLETED=True #if TRUE samples already processed (i.e. with a .bam.bai) will be skipped

#----------- IMPORT
import subprocess
import glob
import os
import pandas as pd
#instead of basic Pool, for complicated reasons linked to shared memory
#that would prevent pandas to be pickable, we use ThreadPool 
from multiprocessing.pool import ThreadPool 

#----------- SUPPORT FUNCTIONS
#this function is designed to be executed in parallel, once per 
#input fastq R1/R2 files
def align(infile_R1, outfolder, reference_file, picard_command):
	#--------- filenames
	fn = create_filenames(infile_R1, outfolder, reference_file)
	
	#--------- bowtie2 align
	cmd = 'bowtie2 -x ' + fn['bowtie_index_base']
	cmd += ' -1 ' +  fn['infile_R1']
	cmd += ' -2 ' +  fn['infile_R2']
	cmd += ' -S ' +  fn['tmp_sam']
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log_bowtie2_align'], "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools for sam -> bam conversion
	cmd = 'samtools view'
	cmd += ' -bS ' +  fn['tmp_sam']
	cmd += ' > ' +  fn['tmp_bam']
	subprocess.run(cmd, shell=True)
	
	#--------- picard, AddOrReplaceReadGroups
	cmd = picard_command + ' AddOrReplaceReadGroups'
	cmd += ' -I ' +  fn['tmp_bam']
	cmd += ' -O ' +  fn['tmp_bam_groups']
	cmd += ' -LB Whatever -PL Illumina -PU Whatever -SM ' + fn['core']
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log_picard_readGroups'], "w") as fp:
		fp.write(res.stdout)

	#--------- picard, ValidateSamFile
	cmd = picard_command + ' ValidateSamFile'
	cmd += ' -INPUT ' +  fn['tmp_bam_groups']
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log_picard_validation'], "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools, sort
	cmd = 'samtools sort ' + fn['tmp_bam_groups']
	cmd += ' -o ' +  fn['outfile']
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log_samtools_sort'], "w") as fp:
		fp.write(res.stdout)
	
	#--------- samtools, index
	cmd = 'samtools index ' + fn['outfile']
	res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	with open(fn['log_samtools_index'], "w") as fp:
		fp.write(res.stdout)
	
	#--------- cleanup of intermediate files
	subprocess.run('rm ' + fn['tmp_sam'], shell=True)
	subprocess.run('rm ' + fn['tmp_bam'], shell=True)
	subprocess.run('rm ' + fn['tmp_bam_groups'], shell=True)

	#--------- done
	return(fn['core'])

#returns a dictionary with all the filenames derived from the input R1
#and the outfolder
def create_filenames(infile_R1, outfolder, reference_file):
	res = {'infile_R1' : infile_R1}
	
	#core sample name, without path and file extension
	res['core'] = os.path.basename(infile_R1).replace('_R1.fastq.gz', '')
	
	#the other input file
	res['infile_R2'] = infile_R1.replace('_R1', '_R2')

	#intermediate files, deletable at the end
	res['tmp_sam'] = outfolder + '/' + res['core'] + '.sam'
	res['tmp_bam'] = res['tmp_sam'].replace('.sam', '.bam')
	res['tmp_bam_groups'] = res['tmp_sam'].replace('.sam', '.gr.bam')
	
	#log files
	res['log_bowtie2_align']     = res['tmp_sam'].replace('.sam', '.bowtie2_align.log')
	res['log_picard_validation'] = res['tmp_sam'].replace('.sam', '.picard_validation.log')
	res['log_picard_readGroups'] = res['tmp_sam'].replace('.sam', '.picard_readGroups.log')
	res['log_samtools_sort']     = res['tmp_sam'].replace('.sam', '.samtools_sort.log')
	res['log_samtools_index']    = res['tmp_sam'].replace('.sam', '.samtools_index.log')
	
	#bowtie2 reference genome name
	res['bowtie_index_base'] = reference_file.replace('.fasta.gz', '')

	#output file
	res['outfile'] = res['tmp_sam'].replace('.sam', '.gr.sorted.bam')
	
	#index for output file
	res['outfile_index'] = res['outfile'] + '.bai'
	
	return(res)
	
#----------- ACTUAL SCRIPT
#room for output
cmd_str = "mkdir -p " + OUTFOLDER
subprocess.run(cmd_str, shell=True)

#collecting all the arguments for the parallel execution in a pandas df
args = None
skipped = 0
for infile_R1 in glob.glob(INFOLDER + '/*_R1.fastq.gz'):
	#should we skip this file?
	fn = create_filenames(infile_R1, OUTFOLDER, REFERENCE_FILE)
	if os.path.isfile(fn['outfile_index']) and SKIP_PREVIOUSLY_COMPLETED:
		skipped += 1
		print('Skipping previously processed sample ' + fn['core'])
		continue
	
	#the arguments for the current file. Column order is important,
	#it should match the order for the parallel function, since 
	#the arguments are passed as positionals
	args_now = pd.DataFrame({
		'infile_R1' : [infile_R1], 
		'outfolder' : [OUTFOLDER], 
		'reference_file' : [REFERENCE_FILE], 
		'picard_command' : [PICARD_COMMAND]
	})
	
	#storing in a single df
	args = pd.concat([args, args_now])

	#safeguard
	if len(args) >= MAX_SAMPLES:
		break

#do we have something to execute?
if args is None:
	print('All samples skipped, no operation required')
	sys.exit(0)
	
#executing in parallel using multiprocessing module
cnt = 0
with ThreadPool(CORES) as pool:
	for result in pool.starmap(align, args.itertuples(index = False)):
		#result variable contains the core of the processed sample, but
		#we don't want to flood the main screen (there's many log files)
		#so we just do nothing with it
		cnt += 1

#closing interface
print('Samples: ')
print(' - processed: ' + str(cnt))
print(' - skipped because of previous runs: ' + str(skipped))
