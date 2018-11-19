# -*-coding:utf8 -*


####
# Python3 script to generate random position resembling oris:
# The randomized regions have the same chromosomal distribution than the original set (same number of random regions per
# chromosomes than in the input set) and the same length distribution.
# The script also requires a target region as input. The randomized segments will have to fall into this target regions
# (the mappable part of the genomes in our study).
# Appart from that, the random positions are randomly and uniformly drawn in the target regions.
#####


import sys
import pickle
import random


def uniq(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked
################################################
def column(matrix,i):
	return [row[i] for row in matrix ]

###################################################"#

###number of random pos:

ori_file=sys.argv[1]  # a bed file containing the list of origins to randomize
Region_file=sys.argv[2] 
# a bed file containing the the regions in which the randomized oris are allowed to fall (in our
# case, the mappable regions of the human genome

if len(sys.argv)<3:
	print ("NOT ENOUGH ARGUMENTS!! \n\n")
	sys.exit(1)

with open(ori_file,'r') as oris:
	## Reading the replication origin bed file
	data_matrix=[]
	for line in oris.readlines(): 
		data_matrix.append(line.split('\t'))
		
	distrib_matrix={}
	distrib_random_pos={}
	i=0

	ori_chr={}
	ori_chr[data_matrix[0][0]]=1
	old_chro=data_matrix[0][0]

	compt=0
	start=0
	for chro in column(data_matrix,0):
		if chro != old_chro:
			ori_chr[old_chro]=[start,compt-1]
			start=compt
			old_chro=chro

		compt+=1

	ori_chr[old_chro]=[start,compt]

	#Compute and store the number of oris on each chromosome	
	for chro in uniq(column(data_matrix,0)):
		distrib_matrix[chro]=ori_chr[chro][1]-ori_chr[chro][0]
		distrib_random_pos[chro]=0


	distrib_random_pos=distrib_matrix

	for chro in distrib_random_pos:
		distrib_random_pos[chro]+=1


	with open(Region_file,'r') as target_pos:
		
		target_matrix=[]
		for line in target_pos.readlines(): 
			target_matrix.append(line.split('\t'))
		
		for mychr in distrib_random_pos.keys():
			#pick a chromosome
			if int(distrib_random_pos[mychr])>0:
				# Retrieve the number of oris on that chromosome
				compt=distrib_random_pos[mychr]		
				# Compute the size distribution of oris on that chromosome
				distrib_ori_length=[]
				for i in range(ori_chr[mychr][0]-1,ori_chr[mychr][1]):
					distrib_ori_length.append(int(data_matrix[i][2])-int(data_matrix[i][1]))
				# Retrieve the target regions for this chromosome
				this_chr_target=[]
				for row in target_matrix:
					if row[0]==mychr:
						this_chr_target.append(row)
		
##########################################Choose a Position

				#Draw random position until we've reached the number of oris on that chromosome
				while compt>0:
					# draw a random position 						
					this_pos=int(this_chr_target[len(this_chr_target)-1][2])*random.random()
					this_pos=int(this_pos)
					# draw randomly a size for this random region from the size distribution of
					# origins on this chromosome
					this_length=distrib_ori_length[int(len(distrib_ori_length)*random.random())]

					# check that in falls into a targeted regions. If not, draw a new one.
					i=0
					while this_pos>int(this_chr_target[i][2]):
						i+=1
				
					while this_pos<int(this_chr_target[i][1]) or this_pos+this_length>int(this_chr_target[i][2]):
						this_pos=int(this_chr_target[len(this_chr_target)-1][2])*random.random()
						this_length=distrib_ori_length[int(len(distrib_ori_length)*random.random())]
						this_pos=int(this_pos)	
						i=0
						while this_pos>int(this_chr_target[i][2]):
							i+=1

					#print the region position										
					this_stop=this_pos+this_length
					print (mychr,this_pos,this_stop,sep='\t')
					compt-=1
		

