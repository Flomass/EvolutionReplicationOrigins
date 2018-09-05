# -*-coding:utf8 -*

import os
import sys
import pickle
#import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
################################
import datetime


# Program that compute the average SNP density on regions of 1500 bp (default).
# Requires tabix and "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz" file from the 1000 genome
# project phase 3. 
# One also needs to have a directory called ./temp/ to run this script


parser = OptionParser()
parser.add_option("-i","--Bed", dest="Bed_file",action="store",
                  help="file containing Phast cons table")
parser.add_option("-p","--para", dest="para",action="store",default=0,
                  help="in order to run several times the program at the same time")
parser.add_option("-k","--keep_size",dest="keep_size",action="store",default=0,
		help="if we wnat to run this prog on regions of different sizes")
parser.add_option("-w","--wind_size",dest="wind_size",action="store",default=750,
		help="wind size around Ori peak where we should look for SNPs")


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


(options, args) = parser.parse_args()
patt=re.compile(",")
para=options.para
Bed_file=options.Bed_file
keep_size=int(options.keep_size)
wind_size=int(options.wind_size)

count_SV=0

wind_size=wind_size*2+2

mean_all=[0]*wind_size
count=0


mean_all_INDELs=[0]*wind_size
mean_rare_INDELs=[0]*wind_size
mean_common_INDELs=[0]*wind_size

SNP = defaultdict(dict)

for i in ["Cto","Gto","Ato","Tto"]:
	for j in ["T","A","C","G"]:
		mutTyp=i+j
		mutTyp2=mutTyp+"_Common"
		mutTyp3=mutTyp+"_Rare"
		SNP[mutTyp]=[0]*wind_size
		SNP[mutTyp2]=[0]*wind_size
		SNP[mutTyp3]=[0]*wind_size

SNP["INDEL_Rare"]=[0]*wind_size
SNP["INDEL"]=[0]*wind_size
SNP["INDEL_Common"]=[0]*wind_size

with open (Bed_file,'r') as Bed:
	for line in Bed.readlines():
		line=line.rstrip()
		ori_coord=line.split('\t')
		ori_coord[0]=ori_coord[0].replace("chr","")
		if keep_size==0:
			start=int((int(ori_coord[1])+int(ori_coord[2]))/2-((wind_size-2)/2))
			stop=int((int(ori_coord[1])+int(ori_coord[2]))/2+((wind_size-2)/2))
		else:
			start=int(ori_coord[1])
			stop=int(ori_coord[2])


		comm=str("tabix ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz "+ori_coord[0]+":"+str(start)+"-"+str(stop)+" >temp/tabix"+str(para))
		os.system(comm)

		tabix_file=open("temp/tabix"+str(para),'r')
		count=count+1
		if int(count/1000)==count/1000:
			date = datetime.datetime.now()
			eprint(str(count)+" lines have been processed;"+str(date)+"; ori: "+ori_coord[0])

		for tabix_line in tabix_file.readlines():
			tabix_row=tabix_line.split('\t')
			pos=int(tabix_row[1])-start

			if pos < wind_size:

				infos=tabix_row[7]
				A_ref=tabix_row[3]
				A_alt=tabix_row[4]

				AF=re.search(";AF=[0-9\.,]+",infos)
				AF=AF.group(0)
				AF=AF.replace(";AF=","")

				AFs=AF.split(",")

				AF=re.search("AF=[0-9\.,]+",infos)
				AF=AF.group(0)
				AF=AF.replace("AF=","")

				if re.search(',',AF):
					AFs=AF.split(",")
					AF=0
					for i in range(0,len(AFs)):
						AF=AF+float(AFs[i])

				AF=float(AF)
				mut_type=re.search("VT=[A-z]+",infos)
				mut_type=mut_type.group(0)
				mut_type=mut_type.replace("VT=","")

				# Excluding SVs

				if mut_type == "SV":
					count_SV=count_SV+1
				
				else:			

				### All alleles
					mean_all[pos]=mean_all[pos]+1
				
					if mut_type == "SNP":

						# SNP polarization
						AA=re.search(";AA=[A-z]",infos)
						if AA:
							AA=AA.group(0)
							AA=AA.replace(";AA=","")
							AA=AA.upper()
							A_ref=A_ref.upper()
							
							A_alt=A_alt.split(',')
							
							for var in range(0,len(A_alt)):		
								
								mutTyp=""
								mutTyp2=""
								if AA==A_ref:
									A_mut=A_alt[var].upper()
									this_AF=float(AFs[var])
								else:						
									A_mut=A_ref
									this_AF=1-float(AFs[var])

								if re.search("[Cc]",AA):
									mutTyp="Cto"
								elif re.search("[Gg]",AA):
									mutTyp="Gto"
								elif re.search("[Tt]",AA):
									mutTyp="Tto"
								elif re.search("[Aa]",AA):
									mutTyp="Ato"
								
								if mutTyp:
									if re.search("[Cc]",A_mut):
										mutTyp=mutTyp+"C"
									elif re.search("[Gg]",A_mut):
										mutTyp=mutTyp+"G"
									elif re.search("[Aa]",A_mut):
										mutTyp=mutTyp+"A"
									elif re.search("[Tt]",A_mut):
										mutTyp=mutTyp+"T"

									if this_AF>0.1:
										mutTyp2=mutTyp+"_Common"
										SNP[mutTyp2][pos]=SNP[mutTyp2][pos]+1
									elif this_AF<0.01:
										mutTyp2=mutTyp+"_Rare" 		
										SNP[mutTyp2][pos]=SNP[mutTyp2][pos]+1

									SNP[mutTyp][pos]=SNP[mutTyp][pos]+1	
							
					else:
						SNP["INDEL_Common"][pos]=SNP["INDEL_Common"][pos]+1

					## High freq alleles
						if AF>0.1 :
							if mut_type != "SNP":
								SNP["INDEL_Common"][pos]=SNP["INDEL_Common"][pos]+1

					##### Low freq alleles

						if AF<0.01 :
							if mut_type != "SNP":
								SNP["INDEL_Rare"][pos]=SNP["INDEL_Rare"][pos]+1

				


print("#count\tcount_SVs")
print("#"+str(count)+"\t"+str(count_SV))

print("pos",end="")
for cle in SNP:
	print("\t"+cle,end="")

print("")

for pos in range(0,(wind_size-2)):
	print(str(pos),end='')
	for cle in  SNP:
		print("\t"+str(SNP[cle][pos]/count),end='')
	print("")

