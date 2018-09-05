# -*-coding:utf8 -*

######################################################
#
# Script to calculate binned replication landscape
# Requires origin position, mapped read (as bed files) and a "region" files -- usually the mappable portion of the
# studied genomes (this file is called karyo file and should be a bed file).
# To generate figure 2 of the paper, sort the output file by 4th column and do a cumulative plot
#
#
##########################################################


import os
import sys
import pickle
#import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
################################



parser = OptionParser()
parser.add_option("-k","--karyo", dest="karyo_file",action="store",
                  help="file containing chr length -- or mappable regions of the genome")
parser.add_option("-r","--reads", dest="reads_file",action="store",
                  help="file containing SNS reads")
parser.add_option("-s", "--species",dest="sp",action="store",
                  help="species studied")
parser.add_option("-w", "--wind_length",dest="wind_length",type="int",action="store",default=1000,
                  help="window size to calculate density. default is 1000")
parser.add_option("-o", "--output",dest="output_prefix",type="str",action="store",default="./",
                  help="prefix for output file")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
sp=options.sp
karyo_file=options.karyo_file
wind_length=options.wind_length
reads_file=options.reads_file
output_prefix=options.output_prefix


f=open("temp.bed","w")
with open(karyo_file,'r') as karyo:
	for line in karyo.readlines():
		line=line.rstrip()
		row=line.split("\t")
		if (row[0] != "chrY"):
			for i in range(int(row[1]),int(row[2]),int(wind_length)):
				f.write(row[0]+"\t"+str(i)+"\t"+str(i+wind_length)+"\n")

f.close()		
comm="sort -k1,1 -k2,2n temp.bed >temp2.bed"
os.system(comm)

comm2="intersectBed -a temp2.bed -b "+reads_file+" -c -sorted >"+output_prefix+"_density_wind"+str(wind_length)+"_"+sp+"_w_ori_strength.bed"
os.system(comm2)


