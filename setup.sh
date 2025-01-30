#!/bin/bash

#This is an installation script for all software required by legpipe2
#Its intended scope is to bring a clean slate, ubuntu 64 (aws based)
#machine from mint conditions to ready-to-run. Your mileage may vary.
#Worst case: use the script as a shopping list for everything you need
#to install

#software base folder, where to save stuff that is not installed via apt
LEGPIPE_ROOT="~/software"
mkdir -p "$LEGPIPE_ROOT"

#we add software paths to ~/.bashrc
#change it if you have a different configuration
PATHFILE='~/.bashrc'
echo export LEGPIPE_ROOT=~/software >> "$PATHFILE"

#minimal ubuntu setup
echo "------------------ UBUNTU base updates"
sudo apt-get update
sudo apt-get upgrade
sudo apt install -y parallel
sudo apt install -y fastp

#python3 and modules
echo "------------------ python install"
sudo apt install -y python3
sudo apt install -y python3-pip
sudo apt install python3.12-venv
python3 -m venv "$LEGPIPE_ROOT"/legpipe_venv
source "$LEGPIPE_ROOT"/legpipe_venv/bin/activate
pip install pandas
pip install Bio
deactivate

#tabix for bgzip
echo "------------------ tabix"
sudo apt install -y tabix

#fastx_trimmer
echo "------------------ fastx_trimmer"
#check last version on http://hannonlab.cshl.edu/fastx_toolkit/download.html 
cd $LEGPIPE_ROOT
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
mv bin fastx_folder
echo PATH=$LEGPIPE_ROOT/fastx_folder:$PATH >> $PATHFILE

#bowtie2
echo "------------------ bowtie2"
sudo apt install -y bowtie2

#samtools, here for version 1.3.1, update accordingly
echo "------------------ samtools"
sudo install samtools

#picard
echo "------------------ picard"
cd $LEGPIPE_ROOT
wget https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar
#add this to the path file, so that to allow commands like "java -jar $PICARD ..."
echo export PICARD=$LEGPIPE_ROOT/picard.jar >> $PATHFILE

#GATK
echo "------------------ gatk"
sudo apt install -y openjdk-17-jre-headless
sudo apt install -y unzip
#No installation require for GATK, download version 4.3 and just change the python/python3 shebang
cd $LEGPIPE_ROOT
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip 
unzip gatk-4.3.0.0.zip
#edit the very first line of ./gatk from "#!/usr/bin/env python" to "#!/usr/bin/python3"
sed -i 's|#!/usr/bin/env python|#!/usr/bin/python3|' gatk-4.3.0.0/gatk
echo PATH=$LEGPIPE_ROOT/gatk-4.3.0.0:$PATH >> $PATHFILE

#bcftools
echo "------------------ bcftools"
sudo apt install -y bcftools

#done
echo Several software have been installed, either using local apt manager or manually installing in $LEGPIPE_ROOT
echo Paths have been added to $PATHFILE, please reload \(or just open a new terminal\)
echo Remember to activate the python virtual environment via: python -m venv "$LEGPIPE_ROOT"/legpipe_venv
