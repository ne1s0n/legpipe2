# legpipe2: a SNP calling pipeline

The SNP calling pipeline used at CREA-ZA (Lodi, Italy). Main application is legumes (diploid and tetraploid) and restricted sequencing (mainly Genotyping By Sequencing). You can use it for everything else, though.

## Status

Work in progress & not ready to run. Stay tuned.

## Main inspirations

- [dDocent](https://github.com/jpuritz/dDocent)
- [UGbS-Flex](https://github.com/madgenetics/UGbS-Flex)

## TODO - Big chunks

Modules still missing that need to be implemented

- write installation tutorial

## TODO - Functionalities

Checks and features that ensure the pipeline fails gracefully
when something is wrong

- support both .fa.gz and .fasta.gz genome files (see in particular genome indexing module)
- the picard jar should be read in the config file like all other configs, it does not make sense to use an env variable only of this)
	- creating maybe a common section in the .ini
- check if the reference genome is gzipped and not bgzipped
- indexing module needs a better management of the logs

## WISHLIST

Stuff that it would be nice to have, once the above blocks are empty

- export to conf file all the commands, especially from align/filter
- we are using both samtools and bcftools. Is it really necessary?
- call module does a lot of work. it could probably be split in three
  modules (HaplotypeCaller, GenomicsDBImport and GenotypeGVCFs). On the
  other hand, I'm not sure how common is to do each step separatedly
- implement test for required software (maybe module based? like what
  is done for `interpolate()` and `validate()` functions). Also check for python 3.6+
- post calling filtering from UGbS-Flex:
	- [optional] Remove adjacent SNPs
	- [optional] Consolidate SNPs
	- [optional] Select SNPs based on parental scores (only applies to some mapping populations)
	- [optional] Remove cosegregating SNP markers
- remove any subprocess.run(..., shell=True)
