# legpipe2: a SNP calling pipeline

The SNP calling pipeline used at CREA-ZA (Lodi, Italy). Main application is legumes (diploid and tetraploid) and restricted sequencing (mainly Genotyping By Sequencing). You can use it for everything else, though.

## Status

Work in progress & not ready to run. Stay tuned.

## Main inspirations

- [dDocent](https://github.com/jpuritz/dDocent)
- [UGbS-Flex](https://github.com/madgenetics/UGbS-Flex)

## TODO

- test align module
- implement post-execution script (for synch, shutdown...)
- implement post calling filtering:
	- select only biallelic SNPs and filter on allele frequency (MAF)
	- [optional?] Remove adjacent SNPs
	- [optional?] Consolidate SNPs
	- [optional?] Select SNPs based on parental scores (only applies to some mapping populations)
	- [optional?] Remove cosegregating SNP markers

## WISHLIST

- call module does a lot of work. it could probably be split in three
  modules (HaplotypeCaller, GenomicsDBImport and GenotypeGVCFs). On the
  other hand, I'm not sure how common is to do each step separatedly
- implement test for required software (maybe module based? like what
  is done for `interpolate()` and `validate()` functions)

