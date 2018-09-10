# Evolution of Replication Origins in Vertebrates

You will find in this repository the data and code described in the article "Evolution of Replication Origins in
Vertebrate Genomes: Rapid Turnover Despite Selective Constraints" which is available as a preprint on the biorXiv
platform [here](https://www.biorxiv.org/content/early/2018/09/10/411470)



# Data

The data folder contains the position of replication origins (bed files).
Some dataset comes with the SNS peak position (i.e. the loci maximizing the read accumulation, see Methods of the paper
for more details) and the strength of the origin (defined as the raw count of reads mapping in the origin -- it might be
a good idea to renormalize by origin size before interpreting this column).

## Human data

3 subfolder in this directory. Raw data come from (Besnard et al 2012 nat. str. and mol. biol.) and (picard et al. 2014
plos. biol.)

* Oris_no_cov_correction  
Contains one set of Origins per technical replicate for each of the 5 available cell line

* Oris_seq_depth_32M 
Contains one set of Origin per human cell line. In these sets, mapped read where first combined and
subsampled in order to generate 5 datasets (one per cell line) with equal sequencing depth (32M reads per set).
Data in this folders should be used to compare origins of different cell types. 

* Oris_for_vertebrate_comparison_66M_reads
Contains the H9 origin data set that was used to compare to the mouse and chicken datasets.
Detection was performed on a resampled read set so that experiment in the 3 species have similar coverage (i.e between
24 and 28 reads per kb)


## Mouse data

2 subfolder in this directory:
Raw data come from (almeida et al 2018, nat.comm.) or (Cayrou et al. 2015 Gen. Res.)

* Oris_no_cov_correction  
Contains one set of Origins per technical replicate for the 2 cell lines available (mESC and MEF) in mouse.

* Oris_for_vertebrate_comparison_50M_reads
Contains the mESC origin (original data in Cayrou et al, 2015) data set that was used to compare to the mouse and
chicken datasets. Detection was performed on a resampled read set so that experiment in the 3 species have similar
coverage (i.e between 24 and 28 reads per kb)

## Chicken Data

*Data will be made available upon publication*

2 subfolder in this directory. Raw data come from (massip et al., soon to be submitted).
* Oris_no_cov_correction  
Contains one set of Origins as there is only one technical replicate available for chicken (DT40 cell line).

* Oris_for_vertebrate_comparison_50M_reads
Contains the DT40 origin data set that was used to compare to the mouse and human datasets. Detection was performed on a
resampled read set so that experiment in the 3 species have similar coverage (i.e between 24 and 28 reads per kb).
Raw reads were subsampled 3 times independently hence leading to 3 different origin sets.


# Scripts

**Warning:  These scripts are intended to give an idea/ help reproduce the analysis described in the mentionned paper,
but are most of the time not worknig as it is / will require a little work to be adapted**

Each script file has a header quickly explaining what the script does and what are the required input files.

Five files are provided to:
* Compute the position of the SNS peak

* Compute the average SNP / INDEL density around the SNS peak

* Compute the PhastCons scores around the SNS peak

* Compute the mean nucleotide content of fixed size regions

* Compute the replication landscape (fig2)

What I do not provide here (but is available on demand):

* Scripts to assess evolutionary conservation of origins. This mainly require instersecting origins w/ UCSC pairwise
alignment using bedtools and is relatively easy. 

* Scripts to plot the actual data -- this is available on demand.

