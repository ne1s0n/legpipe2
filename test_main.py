#incantation to import the legpipe2 module wherever you execute this script
import os 
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/legpipe2')

#ready to import
import legpipe2

conf = legpipe2.Legpipe2('sample_config.ini')

#print(conf)
