#####################################################################################################
#
# A script to compute the mean phastCons score of Origins or TSS (or any other bed type data).
#
# Requires a bedfile as an entry, as well as the position of TE on the genome (bed file format as well)
# and the position of G4s along the genome (-- bed again) 
# This last part is not essential for the analysis of the paper, but the script won't work in its current state if you
# don't provide these.
#
#  you'll also need the phastcons data (replace lines
#  "/home/florian/work_ori/PhastCons/all_data_human_Primates/"+chrom+".phastCons46way.primates.wigFix" by the path to
#  you directory containing PhastCons files)
#
# These data can be downloaded easily on the UCSC website
#
#####################################################################################################


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
parser.add_option("-p", "--peak",
                  action="store_true", dest="peak", default=False,
                  help="set to True if you provide a peak in the ori_file (peak should be the 5th row)")
parser.add_option("--peak_file",
                  action="store", dest="peak_file", default='',
                  help="when peak is true, is this the file storing the peaks")
parser.add_option("--out", dest="out_file",action="store",
                  help="outfiles will be <>_chr1:21 and <>_mean")
parser.add_option("-a", "--alignment",dest="ali",action="store",
                  help="which set of alignments do you want to use? possibilities are placental or all or primates")
parser.add_option("-w", "--window_size",dest="wind",type="int",action="store",default=50,
                  help="size of the window. To avoid trouble, give an even number. default is 50")
parser.add_option("-s", "--species",dest="sp",action="store",
                  help="species studied")
parser.add_option("-e", "--extension",dest="ext",type="int",action="store",default=500,
                  help="size of region to consider around each side of the ori. default is 500")
parser.add_option("-l", "--ori_length",dest="ori_length",type="int",action="store",default=500,
                  help="size of region to consider around the peak (ori size). default is 500. Better be even")
parser.add_option("-c", "--col_number",dest="col_number",type="int",action="store",default=3,
                  help="Number of column of the input bedfile that should appear in the output. deffault is 3 but if you
want to keep additional information in your bed file in your output you can")
parser.add_option("-r", "--random",dest="rando",type="int",action="store",default=0,
                   help="if different from 1, randomly move the center of the input peak (or mean position of the feature) of +/- the value given for r. This allows to consider windows of 1500 bp with a TSS positionned -- to imitate what we do for oris")
parser.add_option("--accept_TE",dest="accept_TE",action="store_true",default=False,
                   help="if set to true, keep the phastcons score of TE regions")



(options, args) = parser.parse_args()


patt=re.compile("start=([0-9]+)")

sp=options.sp
wind=options.wind
peak=options.peak
ali=options.ali
out_file=options.out_file
ori_file=options.ori_file
ext=options.ext
ori_length=options.ori_length
peak_file=options.peak_file
col_number=options.col_number
rando=options.rando
accept_TE=options.accept_TE

ch_patt=re.compile("chr(.+)")
inv_patt=re.compile("strand|inv")

TE_temp=[]

chrom=''
TE_row=["chr0",0,0]

out_file_mean=out_file+"_mean"
f_mean=open(out_file_mean,'w')
if (peak):
	peak=open(peak_file,'r')

with open (ori_file,'r') as ori:
	legend=ori.readline()

## Check whether there is an inv column
	if (inv_patt.search(legend)):

		legend=legend.rstrip()
		legend_row=legend.split('\t')

		f_mean.write("ori\tlength\tmean\tdownstream_mean\tupstream_mean\tmax\tpos_max\tmin\tpos_min")
		for n in range(3,col_number):
			f_mean.write('\t'+legend_row[n])
		f_mean.write("\n")
	else:
		sys.exit("there has to be an 'inv' column in the ori file!\n You can generate it with ~/work_ori/conserved_oris/scripts/orientating_ori_according_to_GGG.sh for instance\n")

	for line in ori.readlines():
		line=line.rstrip()
		row=line.split('\t')

		if (peak):
			this_peak=int(peak.readline())
		else:
			this_peak=int((int(row[1])+int(row[2]))/2)

		if (rando>0):
			this_rando=int(random.random()*2*rando)-rando
			this_peak=this_peak+this_rando

		row[1]=int(this_peak-ori_length/2)
		row[2]=int(this_peak+ori_length/2)

		if (chrom!=row[0]):			
			phastPos=0
			chrom=row[0]
			if sp=="human":
				if ali=="primates":
					file_phast=open("/home/florian/work_ori/PhastCons/all_data_human_Primates/"+chrom+".phastCons46way.primates.wigFix","r")
				elif ali=="placentals":
					file_phast=open("/home/florian/work_ori/PhastCons/all_data_human_Placentals/"+chrom+".phastCons46way.placental.wigFix","r")

				elif ali=="all":
					file_phast=open("/home/florian/work_ori/PhastCons/all_data_human_46_sp/"+chrom+".phastCons46way.wigFix","r")
				else:
					print ("wrong ali\n")
					sys.exit(1)

				TEcomm=str("grep '"+chrom+"\s' /home/florian/work_ori/conserved_oris/data/human/TE/TE_hg19_no_simple_and_low_complexity.bed |sort -k1,1 -k2,2n |mergeBed >~/work_ori/PhastCons/temp/TE_human")
				os.system(TEcomm)
				TE_file=open("/home/florian/work_ori/PhastCons/temp/TE_human",'r')
				TE=TE_file.readline()
				TE_row=["chr0",0,0]
				m=ch_patt.search(chrom)
				short_chr=m.group(1)
				G4_file=open("/home/florian/work_ori/conserved_oris/data/human/G4/Homo_sapiens.GRCh37.75.dna.chromosome."+short_chr+".fa.G4.bed",'r')
				G4_row=["chr0",0,0]

				triplet_file=open("/home/florian/work_ori/conserved_oris/data/human/G4/Homo_sapiens.GRCh37.75.dna.chromosome."+short_chr+".fa.G4_triplet_pos",'r')
				triplet_row=["chr0",0,0]

				f = open(str(out_file+"_"+chrom), 'w')
				f_G4=open(str(out_file+"_"+chrom+"_G4"),'w')
				f_triplet=open(str(out_file+"_"+chrom+"_G4_triplet"),'w')
				print("I opened "+chrom)
			
				f.write("ori")
				f_G4.write("ori")
				f_triplet.write("ori")
				for n in range(1,2*ext+ori_length+1):
					f.write("\tpos"+str(n))
					f_G4.write("\tpos"+str(n))
					f_triplet.write("\tpos"+str(n))
				for n in range(3,col_number):
					f.write('\t'+legend_row[n])
					f_G4.write('\t'+legend_row[n])
					f_triplet.write('\t'+legend_row[n])

				f.write("\n")
				f_G4.write("\n")
				f_triplet.write("\n")	
				TEs=[]	
				for line in TE_file.readlines():
					line=line.rstrip()
					TEs.append(line.rsplit('\t'))

				last_TE="chr0\t450000000\t450000002"
				TEs.append(last_TE.rsplit('\t'))
				TE_pos=0
				previous_TE_pos=0

			elif sp=="mouse":
				if ali=="placental":
					file_phast=open("/home/florian/work_ori/PhastCons/all_data_mouse_Placental/"+chrom+".phastCons60way.placental.wigFix","r")
				elif ali=="all":
					file_phast=open("/home/florian/work_ori/PhastCons/all_data_mouse_60_sp/"+chrom+".phastCons60way.wigFix","r")
				else:
					print ("wrong ali\n")
					sys.exit(1)

				TEcomm=str("grep '"+chrom+"\s' /home/florian/work_ori/conserved_oris/data/mouse/TE/TE_mm10_no_simple_and_low_complexity.bed |sort -k1,1 -k2,2n |mergeBed >~/work_ori/PhastCons/temp/TE_mouse")
				os.system(TEcomm)
				TE_file=open("/home/florian/work_ori/PhastCons/temp/TE_mouse",'r')
				TE=TE_file.readline()
				TE_row=["chr0",0,0]
				m=ch_patt.search(chrom)
				short_chr=m.group(1)
				G4_file=open("/home/florian/work_ori/conserved_oris/data/mouse/G4/Mus_musculus.GRCm38.dna.chromosome."+short_chr+".fa.G4.bed",'r')
				G4_row=["chr0",0,0]
				f = open(str(out_file+"_"+chrom), 'w')
				f_G4=open(str(out_file+"_"+chrom+"_G4"),'w')
				f.write("ori")
				f_G4.write("ori")
				for n in range(1,2*ext+ori_length+1):
					f.write("\tpos"+str(n))
					f_G4.write("\tpos"+str(n))
				for n in range(3,col_number):
					f.write('\t'+legend_row[n])
					f_G4.write('\t'+legend_row[n])

				f.write("\n")
				f_G4.write("\n")
				print("I opened "+chrom)
				TEs=[]	
				for line in TE_file.readlines():
					line=line.rstrip()
					TEs.append(line.rsplit('\t'))

				last_TE="chr0\t450000000\t450000002"
				TEs.append(last_TE.rsplit('\t'))
				TE_pos=0
				previous_TE_pos=0

			else:
				print(sp+" is not a valid specie name!\n")
				sys.exit(1)
	
			

		if int(row[1])-ext-int(wind/2+1)<phastPos:
			phastPos=previous_ori_phastPos
			file_phast.seek(previous_ori_file)
		if   int(row[1])-ext-int(wind/2+1)<int(TEs[TE_pos][2]):
			TE_pos=previous_TE_pos

		previous_ori_file=file_phast.tell()
		previous_ori_phastPos=phastPos
		previous_TE_pos=TE_pos

		while phastPos<int(row[1])-ext-int(wind/2+1):
			line_phast=file_phast.readline()			

			line_phast=line_phast.rstrip()
			if line_phast=='':
				phastPos=450000000
				line_phast="0"
			if patt.search(line_phast):
				m = patt.search(line_phast)
				phastPos=int(m.group(1))
			else:
				phastPos+=1

		this_ori_phastCons=defaultdict(list)
		while phastPos<int(row[2])+ext+int(wind/2+1):
			line_phast=file_phast.readline()			
			if line_phast=='':
				phastPos=450000001
				line_phast="0"


			if patt.search(line_phast):
				m = patt.search(line_phast)
				phastPos=int(m.group(1))
			else:
				while phastPos>int(TEs[TE_pos][2]):
					TE_pos+=1
					
				if accept_TE:
					this_ori_phastCons[phastPos]=float(line_phast)
				else:
					if phastPos not in range(int(TEs[TE_pos][1]),int(TEs[TE_pos][2])+1):
						this_ori_phastCons[phastPos]=float(line_phast)
					else:
						this_ori_phastCons[phastPos]=-1

				phastPos+=1


		nb_for_this_mean=0
		nb_for_this_up=0
		nb_for_this_down=0
		this_mean=0
		this_downstream_mean=0
		this_upstream_mean=0
		this_max=-1000
		this_max_pos=-1000-ext
		this_min=1000
		this_min_pos=-1000-ext

		row[2]=str(row[2])
		row[1]=str(row[1])
		f.write(str(row[0])+"-"+str(row[1])+"-"+str(row[2]))
		f_G4.write(str(row[0])+"-"+str(row[1])+"-"+str(row[2]))
		f_triplet.write(str(row[0])+"-"+str(row[1])+"-"+str(row[2]))
		for i in range(int(row[1])-ext,int(row[2])+ext):
			this_phast=0
			nb_for_this_wind=0
			for j in range(int(i-wind/2),int(i+wind/2)):
				if j in this_ori_phastCons:
					this_phast+=this_ori_phastCons[j]
					nb_for_this_wind+=1
			if nb_for_this_wind>0:#int(wind/2):
				this_phast=this_phast/nb_for_this_wind
				if this_phast>this_max:
					this_max=this_phast
					this_max_pos=j
				if this_phast<this_min:
					this_min=this_phast
					this_min_pos=j

				f.write("\t"+str(this_phast))
			else:
				f.write("\tNA")
				this_phast="NA"

			while i>int(G4_row[2]):
				G4=G4_file.readline()
				if G4=='':
					G4="chr0\t500000000\t500000000\n"

				G4=G4.rstrip()
				G4_row=G4.split('\t')
			
			if i in range(int(G4_row[1]),int(G4_row[2])):
				f_G4.write("\t"+str(this_phast))
			else:
				f_G4.write("\tNA")



			while i>int(triplet_row[2]):
				triplet=triplet_file.readline()
				if triplet=='':
					triplet="chr0\t500000000\t500000000\n"

				triplet=triplet.rstrip()
				triplet_row=triplet.split('\t')
			
			if i in range(int(triplet_row[1]),int(triplet_row[2])):
				f_triplet.write("\t"+str(this_phast))
			else:
				f_triplet.write("\tNA")

				
			if i in this_ori_phastCons:
				if this_ori_phastCons[i]>=0:
					if i<int(row[1]):
						this_downstream_mean+=this_ori_phastCons[i]
						nb_for_this_down+=1
					elif i in range(int(row[1]),int(row[2])):
						this_mean+=this_ori_phastCons[i]
						nb_for_this_mean+=1
					else:
						this_upstream_mean+=this_ori_phastCons[i]
						nb_for_this_up+=1

		for n in range(3,col_number):
			f.write("\t"+row[n])
			f_G4.write("\t"+row[n])
			f_triplet.write("\t"+row[n])
		f.write("\n")
		f_G4.write("\n")
		f_triplet.write("\n")


		if nb_for_this_mean>(int(row[2])-int(row[1]))/2:
			this_mean=this_mean/nb_for_this_mean
		else:
			this_mean="NA"
		if nb_for_this_up>ext/2:
			this_upstream_mean=this_upstream_mean/nb_for_this_up
		else:
			this_upstream_mean="NA"
		if nb_for_this_down>ext/2:
			this_downstream_mean=this_downstream_mean/nb_for_this_down
		else:
			this_downstream_mean="NA"
		if this_max_pos==-1000-ext:
			this_max="NA"
			this_max_pos="NA"
			this_min="NA"
			this_min_pos="NA"

		f_mean.write(row[0]+"-"+row[1]+"-"+row[2]+"\t")
		f_mean.write(str(int(row[2])-int(row[1]))+"\t")
		f_mean.write(str(this_mean)+"\t")
		f_mean.write(str(this_downstream_mean)+"\t")
		f_mean.write(str(this_upstream_mean)+"\t")
		f_mean.write(str(this_max)+"\t")
		f_mean.write(str(this_max_pos)+"\t")
		f_mean.write(str(this_min)+"\t")
		f_mean.write(str(this_min_pos))
		for n in range(3,col_number):
			f_mean.write("\t"+row[n])
		f_mean.write("\n")



