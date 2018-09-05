#####################################################################################################################
#
# R script that find the position of the SNS peak for origins that have been called already 
# As an input, this scripts require a file w/ the following structure:
#
#  CHR-start-stop\tread_start\tread_end  => One such line for each read mapping in (or close by -- I use 500bp limit) to
#  an origin.
#
# You can generate such a file with the following command (bedtools required, as well a karyotype file for the
# genome that contain the length of each chr) 
#
#  <input.SNS> is a bed file of all mapped read (sorted by chr and pos w/ sort -k1,1 -k2n) 
#  
#    bedtools slop -b 500 -g ../data/genomic_data/human/hg19_karyotype.txt \
#   -i <Ori_bedfile> |sort -k1,1 -k2,2n >temp/ori
#
#   intersectBed -a <input.SNS> -b temp/ori  -sorted -wa -wb  \
#   | awk '{{ print $4 "-" $5+500 "-" $6-500 "\t" $2 "\t" $1 }}' \
#   |sort -k1,1 -k2,2n -t '-' >inFile"""
#
#
# For efficiency, I do this chr by chr (one input file per chr).
# 
#####################################################################################################################

args = commandArgs(trailingOnly=TRUE) 


inFile=args[1]
cellType=args[2]
this_chr=args[3]
outFile=args[4]

temp_read_file=inFile

comm3=paste("cut -f 1 ",temp_read_file, "|uniq",sep='')
oris=read.table(pipe(comm3),sep='-')

comm4=paste("cut -f 1 ",temp_read_file,"|uniq -c",sep='')
keys=read.table(pipe(comm4))


oris$peaks=rep(NA,length(oris$V1))
dist=c()
len=c()
len2=c()
cum_keys=c(cumsum(keys$V1))
inv_cum_keys=rev(cumsum(rev(keys$V1)))

for (i in 1:length(keys$V1)){
  HH=cum_keys[i]
  TT=keys$V1[i]
  ori=keys$V2[i]

  if ( cum_keys[i]<cum_keys[floor(length(cum_keys)/2)] ){
    comm5=paste("head -n",HH,temp_read_file,"|tail -n",TT)
  }   else{
  TT=inv_cum_keys[i]  
  HH=keys$V1[i]
  comm5=paste("tail -n",TT,temp_read_file," 2>Err_files/err_R |head -n",HH)
}
a=read.table(pipe(comm5))
start=oris$V2[i]
stop=oris$V3[i]
xlim1=c(start-500,stop+500)
dens=density(a$V2,bw=sqrt(500))
  
dens$y=dens$y[which(dens$x>start & dens$x<stop )]
dens$x=dens$x[which(dens$x>start & dens$x<stop )]
store=dens$x[which.max(dens$y)]
oris$peaks[i]=round(store)

write.table(file = outFile,quote=FALSE,sep='\t',x = oris[i,],row.name=FALSE,col.name=FALSE,append = T)
  
}
