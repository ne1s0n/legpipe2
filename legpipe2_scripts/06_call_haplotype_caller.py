#SNP calling with GATK Haplotype Caller module. following 
#the GVCF workflow 



INFOLDER='/data/processed/aligned'
OUTFOLDER='/data/processed/called'
REFERENCE_FILE='/data/reference/zm-4.genome.largest_chromosomes.fasta.gz'
PLOIDY=4
EXPERIMENT='alfalfa_30006_tetra'
TMP_FOLDER='/data/tmp'

INFOLDER='/home/nelson/tmp/foo'
OUTFOLDER='/home/nelson/tmp/foo'

import subprocess
import glob
import os


#a subfolder for GATK gvcf files
OUTFOLDER_GVCF = OUTFOLDER + '/GVCF'

#room for output, tmp
cmd_str = "mkdir -p " + OUTFOLDER_GVCF
#subprocess.run(cmd_str, shell=True)
cmd_str = "mkdir -p " + TMP_FOLDER
#subprocess.run(cmd_str, shell=True)

#keeping track of the produced g.vcf.gz files
gvcf_list = []

#------------ HaplotypeCaller
#for each input bam
for infile in glob.glob(INFOLDER + '/*.gr.sorted.bam'):
	#the produced gvcf file
	gvcf = OUTFOLDER_GVCF + '/' + os.path.basename(infile).replace('.gr.sorted.bam', '.g.vcf.gz')
	gvcf_list.append(gvcf)
	
	#https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller
	#gatk --java-options "-Xmx4g" HaplotypeCaller  \
	#   -R Homo_sapiens_assembly38.fasta \
	#   -I input.bam \
	#   -O output.g.vcf.gz \
	#   -ERC GVCF \
	#   -ploidy PLOIDY \
	#   --min-pruning 1 #default is two \
	#   -stand-call-conf 30 #default 
	
	cmd_str = 'gatk --java-options "-Xmx4g" HaplotypeCaller'
	cmd_str += ' -ERC GVCF --min-pruning 1 -stand-call-conf 30' 
	cmd_str += ' -ploidy ' +  str(PLOIDY) 
	cmd_str += ' -R ' +  REFERENCE_FILE 
	cmd_str += ' -I ' +  infile
	cmd_str += ' -O ' +  gvcf
	print(cmd_str)
	subprocess.run(cmd_str, shell=True)


#------------ GenomicsDBImport
#import everything in a genomic db
#https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
#gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
#      -V data/gvcfs/mother.g.vcf.gz \
#      -V data/gvcfs/father.g.vcf.gz \
#      -V data/gvcfs/son.g.vcf.gz \
#      --genomicsdb-workspace-path my_database \
#      --tmp-dir /path/to/large/tmp \
cmd_str = 'gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport '
for g in gvcf_list:
	cmd_str += ' -V ' + g
cmd_str += ' --genomicsdb-workspace-path ' +  EXPERIMENT
cmd_str += ' --tmp-dir ' + TMP_FOLDER

print(cmd_str)
subprocess.run(cmd_str, shell=True)


#------------ GenotypeGVCFs
#joint variant calling
#https://gatk.broadinstitute.org/hc/en-us/articles/9570489472411-GenotypeGVCFs
# gatk --java-options "-Xmx4g" GenotypeGVCFs \
#   -ploidy PLOIDY \
#   -R Homo_sapiens_assembly38.fasta \
#   -V gendb://my_database \
#   -O output.vcf.gz \
#   --tmp-dir /path/to/large/tmp
cmd_str = 'gatk --java-options "-Xmx4g" GenotypeGVCFs'
cmd_str += ' -ploidy ' +  str(PLOIDY) 
cmd_str += ' -R ' +  REFERENCE_FILE 
cmd_str += ' -V ' +  'gendb://' + EXPERIMENT
cmd_str += ' -O ' +  OUTFOLDER + '/raw_SNPs_haplo.vcf.gz'
cmd_str += ' --tmp-dir ' + TMP_FOLDER
print(cmd_str)
subprocess.run(cmd_str, shell=True)
