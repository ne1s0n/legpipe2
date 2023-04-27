#from the filtered vcf file, output data in several matrices
#fornat



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
	
	#reference, for future
	#https://www.reneshbedre.com/blog/vcf-fields.html
	#https://samtools.github.io/bcftools/bcftools.html#query

	#sample list, one per row
	#bcftools query -l filtered_SNPs_haplo.vcf.gz > sample_list.txt

	#SNP info: chromosome, position, ref allele and the first alternate allele
	#bcftools query -f '%CHROM,%POS,%REF,%ALT{0}\n' filtered_SNPs_haplo.vcf.gz > SNP_info.txt

	#make a BED file: chr, pos (0-based), end pos (1-based), id
	#bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' file.bcf > SNP_info.bed

	#general structure for output
	#bcftools view --apply-filters PASS infile.vcf.gz | bcftools query -f <some other filter> | bgzip -c > output.csv.gz

	#allele count matrices: AD, DP
	#bcftools view --apply-filters PASS filtered_SNPs_haplo.vcf.gz | bcftools query -f '[%AD]\n' | bgzip -c > output_AD.csv.gz
	#bcftools view --apply-filters PASS filtered_SNPs_haplo.vcf.gz | bcftools query -f '[%DP]\n' | bgzip -c > output_DP.csv.gz
	
	print('WARNING: this is a stub function')
	return(None)
