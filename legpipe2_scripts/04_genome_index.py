#prepare the reference genome (index with bowtie, use picard to create sequence dictionary)

REFERENCE_FILE='/data/reference/zm-4.genome.largest_chromosomes.fasta.gz'
PICARD_COMMAND='java -jar /home/ubuntu/software/picard.jar'

#example commands:
#bowtie2/2.2.9/bin/bowtie2-build ref_98.fa ref_98
#java -jar <picard.jar> CreateSequenceDictionary -R=ref_98.fa O=ref_98.dict
#samtools faidx ref_98.fa

import subprocess
import glob
import os

# ------------ bowtie
print('===== bowtie =====')
#bowtie needs to run in the reference genome folder, so that
#all the created files stay there
cmd_str = 'cd ' + os.path.dirname(REFERENCE_FILE) + '; '
bowtie_index_base = os.path.basename(REFERENCE_FILE).replace('.fasta.gz', '')
cmd_str += 'bowtie2-build ' + REFERENCE_FILE + ' ' + bowtie_index_base
print(cmd_str)
subprocess.run(cmd_str, shell=True)

# ------------ picard
print('===== picard =====')
picard_out = REFERENCE_FILE.replace('.fasta.gz', '.dict')
cmd_str = PICARD_COMMAND + ' CreateSequenceDictionary'
cmd_str += ' -R ' + REFERENCE_FILE
cmd_str += ' -O ' + picard_out
print(cmd_str)
subprocess.run(cmd_str, shell=True)

# ------------ samtools
print('===== samtools =====')
cmd_str = 'samtools faidx ' + REFERENCE_FILE
print(cmd_str)
subprocess.run(cmd_str, shell=True)


