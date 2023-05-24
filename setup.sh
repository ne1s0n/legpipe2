#!/bin/bash

#software base folder, where to save stuff
LEGPIPE_ROOT=~/software
mkdir -p $LEGPIPE_ROOT

#we add software paths to ~/.bashrc
#change it if you have a different configuration
PATHFILE=~/.bashrc
echo export LEGPIPE_ROOT=~/software >> $PATHFILE

#minimal ubuntu setup
sudo apt-get update
sudo apt-get upgrade
sudo apt install -y parallel
sudo apt install -y fastp

#python3 and modules
#TODO how to check if python is at least version 3.6?
sudo apt install -y python3
sudo apt install -y python3-pip
pip install pandas
pip install Bio

#tabix for bgzip
sudo apt install -y tabix

#fastx_trimmer
#check last version on http://hannonlab.cshl.edu/fastx_toolkit/download.html 
cd $LEGPIPE_ROOT
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
mv bin fastx_folder
echo PATH=$LEGPIPE_ROOT/fastx_folder:$PATH >> $PATHFILE

#Install Bowtie2
sudo apt install -y bowtie2

#Install samtools, here for version 1.3.1, update accordingly
sudo apt install -y make
sudo apt install -y gcc
sudo apt install -y libncurses-dev
sudo apt install -y libz-dev
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.3.1
make
sudo make prefix=/usr/local/bin install
echo PATH=$LEGPIPE_ROOT/samtools-1.3.1:$PATH >> $PATHFILE

#Install picard
cd $LEGPIPE_ROOT
wget https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar
#add thus to ~/.bashrc
echo export PICARD=$LEGPIPE_ROOT/picard.jar >> $PATHFILE

#Install GATK
#GATK asks for java 17, so...
sudo apt install -y openjdk-17-jre-headless

#No installation require for GATK, download version 4.3 and just change the python/python3 shebang
cd $LEGPIPE_ROOT
sudo apt install -y unzip
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip 
unzip gatk-4.3.0.0.zip
#edit the very first line of ./gatk from "#!/usr/bin/env python" to "#!/usr/bin/python3"
sed -i 's|#!/usr/bin/env python|#!/usr/bin/python3|' gatk-4.3.0.0/gatk
echo PATH=$LEGPIPE_ROOT/gatk-4.3.0.0:$PATH >> $PATHFILE

#install bcftools
sudo apt install -y bcftools

#done
echo Several software have been installed, either using local apt manager or manually installing in $LEGPIPE_ROOT
echo Paths have been added to $PATHFILE, please reload \(or just open a new terminal\)
