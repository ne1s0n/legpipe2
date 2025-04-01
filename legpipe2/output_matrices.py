#from the filtered vcf file, output data in several matrices formats
#
#Output:
# - matrix_GT.csv.gz : most likely genotypes. Alleles are separated 
#                      by / or |, following the scheme (for diploid):
#                      0/0 -> the sample is a homozygous reference
#                      0/1 -> the sample is heterozygous (carries both reference and alternate alleles)
#                      1/1 -> the sample is a homozygous alternate
#                      ./. -> No genotype called or missing genotype
#                      In case of a polyploid call the scheme is the same,
#                      but with more digits. E.g. a tetraploid will
#                      have 0/0/0/0 or ./././. or other combinations
# - matrix_DP.csv.gz : overall read depth from all target samples
#                      supporting the genotype call.
# - matrix_AD.csv.gz : refers to the allele depth. AD reports the 
#                      informative reads supporting each allele. AD may 
#                      not always sum to DP.
# - matrix_countRef.csv.gz : from AD, counting the supporting reads 
#                            for the SECOND allele (reference). Meaningful
#                            only for diploid calls
# - matrix_countAlt.csv.gz : from AD, counting the supporting reads 
#                            for the SECOND allele (alternative). Meaningful
#                            only for diploid calls
#
#example of vcf file:
# GT:AD:DP:GQ:PGT:PID:PL:PS       
# 0/0:4,0:4:12:.:.:0,12,117       
# 0/0:1,0:1:3:.:.:0,3,15  
# 0/0:0,0:0:0:.:.:0,0,0   
# 0/0:5,0:5:15:.:.:0,15,197       
# 0/0:3,0:3:9:.:.:0,9,99  
# 0/0:5,0:5:0:.:.:0,0,113 
# 0/0:5,0:5:15:.:.:0,15,156

import common
import subprocess
import os
import gzip

def validate(conf, runtime = False):
	'''validate incoming config parameters from .ini file, plus environmental variables.
	Validation at runtime requires all files to be ready'''
	if runtime:
		if not os.path.exists(conf['output_matrices']['infile']):
			msg = 'Input file does not exist: ' + conf['output_matrices']['infile']
			raise FileNotFoundError(msg)

def interpolate(conf, raw_conf):
	'''transform incoming config parameters from .ini file'''
	conf['output_matrices']['infile'] = conf['output_matrices']['infolder'] + '/filtered_SNPs_haplo.vcf.gz'
	
	return(conf)

def _split_AD(infile):
	'''splits AD file, creates two more matrices with counted reference and alternative alleles'''
	outfileRef = infile.replace('_AD.csv.gz', '_countRef.csv.gz')
	outfileAlt = infile.replace('_AD.csv.gz', '_countAlt.csv.gz')
	
	with gzip.open(infile, 'rt') as fp_in:
		with gzip.open(outfileRef, 'wt') as fp_ref:
			with gzip.open(outfileAlt, 'wt') as fp_alt:
				for line in fp_in:
					#splitting the various AD values, expected as cntRef,cntAlt
					pieces = line.strip().split(' ')
					
					#splitting again
					for p in pieces:
						cnt = p.split(',')
						if len(cnt) == 2:
							fp_ref.write(cnt[0] + ' ')
							fp_alt.write(cnt[1] + ' ')
						else:
							#missing data or something similarly unparsable
							fp_ref.write(' ')
							fp_alt.write(' ')
					fp_ref.write('\n')
					fp_alt.write('\n')

def output_matrices(conf):
	'''main function for outputing matrices'''
	
	#interface
	common.print_step_header('Output matrices')
	
	#should we do something?
	RUN_THIS=conf['output_matrices']['run_this']

	if not RUN_THIS:
		print('SKIPPED')
		return(None)
	
	#let's check if all required files are here
	validate(conf=conf, runtime=True) 
	#core config
	INFILE=conf['output_matrices']['infile']
	OUTFOLDER=conf['output_matrices']['outfolder']
	
	#derived file names
	LOGFILE = OUTFOLDER + '/output_matrices.log'
	OUTFILE_SAMPLES   = OUTFOLDER + '/sample_list.txt'
	OUTFILE_SNP_INFO  = OUTFOLDER + '/SNP_info.csv'
	OUTFILE_SNP_BED   = OUTFOLDER + '/SNP_info.bed'
	OUTFILE_MATRIX_AD = OUTFOLDER + '/matrix_AD.csv.gz'
	OUTFILE_MATRIX_DP = OUTFOLDER + '/matrix_DP.csv.gz'
	OUTFILE_MATRIX_GT = OUTFOLDER + '/matrix_GT.csv.gz'
	
	#room for output
	cmd_str = "mkdir -p " + common.fn(OUTFOLDER)
	subprocess.run(cmd_str, shell=True)
	
	#reference links, for future
	#https://www.reneshbedre.com/blog/vcf-fields.html
	#https://samtools.github.io/bcftools/bcftools.html#query

	#sample list, one per row
	#bcftools query -l filtered_SNPs_haplo.vcf.gz > sample_list.txt
	cmd = ['bcftools', 'query', '-l', INFILE]
	with open(LOGFILE, 'w') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd) + ' > ' + OUTFILE_SAMPLES + '\n')
		log_fp.flush()
		with open(OUTFILE_SAMPLES, "w") as out_fp:
			subprocess.run(cmd, shell=False, stdout=out_fp, stderr=log_fp)

	#SNP info: chromosome, position, ref allele and the first alternate allele
	#bcftools view --apply-filters PASS filtered_SNPs_haplo.vcf.gz | bcftools query -f '%CHROM,%POS,%REF,%ALT{0}\n' filtered_SNPs_haplo.vcf.gz > SNP_info.csv
	cmd1 = ['bcftools', 'view', '--apply-filters', 'PASS', INFILE]
	cmd2 = ['bcftools', 'query', '-f', '%CHROM,%POS,%REF,%ALT{0}\n']
	with open(LOGFILE, 'a') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd1) + ' | ' + ' '.join(cmd2) + ' > ' + OUTFILE_SNP_INFO + '\n')
		log_fp.flush()
		with open(OUTFILE_SNP_INFO, "w") as out_fp:
			out_fp.write('Chromosome,position,reference_allele,alternative_allele\n')
			out_fp.flush()
			p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE, stderr = log_fp)
			subprocess.run(cmd2, shell=False, stdin = p1.stdout, stdout=out_fp, stderr=log_fp)

	#make a BED file: chr, pos (0-based), end pos (1-based), id
	#bcftools view --apply-filters PASS filtered_SNPs_haplo.vcf.gz | bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' > SNP_info.bed
	cmd1 = ['bcftools', 'view', '--apply-filters', 'PASS', INFILE]
	cmd2 = ['bcftools', 'query', '-f', '%CHROM\t%POS0\t%END\t%ID\n']
	with open(LOGFILE, 'a') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd1) + ' | ' + ' '.join(cmd2) + ' > ' + OUTFILE_SNP_BED + '\n')
		log_fp.flush()
		with open(OUTFILE_SNP_BED, "w") as out_fp:
			p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE, stderr = log_fp)
			subprocess.run(cmd2, shell=False, stdin = p1.stdout, stdout=out_fp, stderr=log_fp)

	#allele count matrices: AD, DP
	#general structure for the commands
	#bcftools view --apply-filters PASS infile.vcf.gz | bcftools query -f <some other filter> | bgzip -c > output.csv.gz

	#bcftools view --apply-filters PASS filtered_SNPs_haplo.vcf.gz | bcftools query -f '[%AD]\n' | bgzip -c > output_AD.csv.gz
	cmd1 = ['bcftools', 'view', '--apply-filters', 'PASS', INFILE]
	cmd2 = ['bcftools', 'query', '-f', '[%AD ]\n']
	cmd3 = ['bgzip', '-c']
	
	with open(LOGFILE, 'a') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd1) + ' | ')
		log_fp.write(' '.join(cmd2) + ' | ')
		log_fp.write(' '.join(cmd3) + ' > ' + OUTFILE_MATRIX_AD + '\n')
		log_fp.flush()
		with open(OUTFILE_MATRIX_AD, "w") as out_fp:
			p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE, stderr = log_fp)
			p2 = subprocess.Popen(cmd2, stdin = p1.stdout, stderr = log_fp, stdout = subprocess.PIPE)
			p3 = subprocess.run(cmd3, stdin = p2.stdout, stderr = log_fp, stdout = out_fp)
	
	#same thing, but for DP matrix
	cmd2 = ['bcftools', 'query', '-f', '[%DP ]\n']
	with open(LOGFILE, 'a') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd1) + ' | ')
		log_fp.write(' '.join(cmd2) + ' | ')
		log_fp.write(' '.join(cmd3) + ' > ' + OUTFILE_MATRIX_DP + '\n')
		log_fp.flush()
		with open(OUTFILE_MATRIX_DP, "w") as out_fp:
			p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE, stderr = log_fp)
			p2 = subprocess.Popen(cmd2, stdin = p1.stdout, stderr = log_fp, stdout = subprocess.PIPE)
			p3 = subprocess.run(cmd3, stdin = p2.stdout, stderr = log_fp, stdout = out_fp)

	#same thing, but for GT matrix
	cmd2 = ['bcftools', 'query', '-f', '[%GT ]\n']
	with open(LOGFILE, 'a') as log_fp:
		log_fp.write('------------------------\n')
		log_fp.write(' '.join(cmd1) + ' | ')
		log_fp.write(' '.join(cmd2) + ' | ')
		log_fp.write(' '.join(cmd3) + ' > ' + OUTFILE_MATRIX_GT + '\n')
		log_fp.flush()
		with open(OUTFILE_MATRIX_GT, "w") as out_fp:
			p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE, stderr = log_fp)
			p2 = subprocess.Popen(cmd2, stdin = p1.stdout, stderr = log_fp, stdout = subprocess.PIPE)
			p3 = subprocess.run(cmd3, stdin = p2.stdout, stderr = log_fp, stdout = out_fp)

	#splitting the AD field in reference and alternative
	_split_AD(OUTFILE_MATRIX_AD)
	
	return(None)
