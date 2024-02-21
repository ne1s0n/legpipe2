# legpipe2: a SNP calling pipeline

The SNP calling pipeline used at CREA-ZA (Lodi, Italy). Main application is legumes (diploid and tetraploid) and restricted sequencing (mainly Genotyping By Sequencing). You can use it for everything else, though.

Main features:

- clean separation in modules (alignment, trimmin, SNP calling, ...)
- each module can be (re)executed easily, carries its own configuration
- easy to hack. You don't like a module and want to improve it? Just go to the corresponding file
- internally use of [GVCF workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format) from GATK/HaplotypeCaller suite, allowing for low-memory fast calling

## Status

Usable.

## Dependencies

Legpipe2 internally uses many other tools. The script [setup.sh](setup.sh) allows to install everything
starting from a clean ubuntu machine (here "ubuntu" means: apt package manager, sudo for installation privileges, basing Unix shell).
Even if you are not using [setup.sh](setup.sh) you should check it and use it as a guideline for setting up
your machine.

## Using Legpipe2

Since it's python3 all the way down, you just have to clone the repo. Then do a:

`python3 /path/to/legpipe2.py your_config_file.ini`

The configuration file is the core of the pipeline, since it dictates what steps are going to be executed and on what data. A lot of interesting details are provided in the [sample_config.ini](sample_config.ini) example file.

## Main inspirations

- [dDocent](https://github.com/jpuritz/dDocent)
- [UGbS-Flex](https://github.com/madgenetics/UGbS-Flex)

## TODO

### Big chunks

Nothing :)

### Functionalities

Checks and features that ensure the pipeline fails gracefully
when something is wrong:

- support both .fa.gz and .fasta.gz genome files (see in particular genome indexing module)
- check if the reference genome is gzipped and not bgzipped
- indexing module needs a better management of the logs

### Bugs

None known :)

### Wishlist

Stuff that it would be nice to have, once the above blocks are empty:

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
