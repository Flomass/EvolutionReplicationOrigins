# -*-coding:utf8 -*

################################################
##
## Script to count the average number of bp of regions 
## of the same size. Typically, 1000bp around the SNS peaks
## or 1000 bp around the TSSs
## 
##
##################################################

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
parser.add_option("--ori", dest="ori_file",action="store",
                  help="file containing oris")
parser.add_option("-s",dest="size",type="int",action="store",default=50,
                  help="length of the analyzed sequences")
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
ori_file=options.ori_file
size=int(options.size)

#ori=defaultdict(list)

GC_cont=[0]*(size)
G_cont=[0]*(size)
T_cont=[0]*(size)
A_cont=[0]*(size)
C_cont=[0]*(size)



count=[0]*(size)

with open (ori_file,'r') as ori:
	for line in ori.readlines():
		if re.search("^[^>]",line):
			seq=line.rstrip()
			subseq=seq[i]
			subseq=re.sub("[Nn]",'',subseq)
			if len(subseq)>0:
				GC_cont[i]=len(re.findall("[GCgc]",subseq))+GC_cont[i]
				G_cont[i]=len(re.findall("[Gg]",subseq))+G_cont[i]
				C_cont[i]=len(re.findall("[Cc]",subseq))+C_cont[i]
				A_cont[i]=len(re.findall("[Aa]",subseq))+A_cont[i]
				T_cont[i]=len(re.findall("[Tt]",subseq))+T_cont[i]
				count[i]=count[i]+1
			if int(count[int(size/2)]/1000) == count[int(size/2)]/1000:
				eprint(str(count[int(size/2)])+" lines have been done")


print("pos\tGC_cont\tcount\tA_cont\tC_cont\tG_cont\tT_cont")

for i in range(0,size):
	print(str(i),str(GC_cont[i]),str(count[i]),str(A_cont[i]),str(C_cont[i]),str(G_cont[i]),str(T_cont[i]),sep='\t')
	if (A_cont[i]+C_cont[i]+G_cont[i]+T_cont[i]!=count[i]):
		print ("# smth wrong here"+str(i))
	


